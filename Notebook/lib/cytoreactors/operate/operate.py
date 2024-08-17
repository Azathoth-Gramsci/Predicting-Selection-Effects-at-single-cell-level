from cytoreactors.design.reactors import reactors
from cytoreactors.design.program import TurbidostatProgram, GrowDiluteProgram, Preculture
from cytoreactors.operate.api.manager import Manager
from time import time, sleep
from threading import Thread
import logging
import pandas as pd
import shutil
import guava2data
from datetime import datetime, date
import os

# just to forbid creation of two sessions
session_created = False

### the main class
class Session():

    def __init__(self, virtual_mode=False):
        global session_created
        if session_created:
            raise Exception('Session already running, stop it before starting a new one.')
        self.status = 'created'
        session_created = True
        self.virtual_mode = virtual_mode
        self.manager = Manager(virtual_mode=virtual_mode)
        if self.virtual_mode:
            self.loop_duration_s = 120. / 1000.
        else:
            self.loop_duration_s = 120.
        self.programs = {} # dict indexed by reactor id
        self.loop_counter = 0
        self.sampling_schedule = []
        self.samplings = []
        self.ot2_remaining_tips_by_8 = 24 # by 8, so 12 per rack
        self.ot2_state_sampling_plates = {1:[13], 2:[13]} # store cols used in list for each plate. we put a 13 ebcause we go downwards and because it helps my code (ugly hack)
        self.ot2_state_reservoir = {c:0 for c in range(1,13)} # store number of times reservoir col used for wash
        self.additional_cytoplate_dilution = False # to dilute 20X more (not needed if OD below 1.0)
        self.guava_num_events = 5000
        self.should_wash = False
        self.waiting_for_guava_data = False
        self.run_thread = None
        self.reactor_ids_to_clean = []

    @property
    def programs_with_cytometry(self):
        return {rid:program for rid,program in self.programs.items() if program.active_cytometry}

    # from excel file
    def load_programs(self):
        df_programs = pd.read_excel('description-programs.xlsx')
        for _,df_prog in df_programs.iterrows():
            #
            prec = Preculture(strain_id=df_prog['preculture_strain_id'],
                            strain_name=df_prog['preculture_strain_name'],
                            media=df_prog['preculure_media'])
            # decide on program type based on OD_low and OD_high
            if df_prog['OD_low'] == df_prog['OD_low']:
                prog = TurbidostatProgram(user=df_prog['user'],
                                          campaign=df_prog['campaign'],
                                          short_name=df_prog['short_name'],
                                          description=df_prog['description'],
                                          reactor_id=df_prog['reactor_id'],
                                          preculture=preculture,
                                          media=df_prog['media'],
                                          OD_setpoint=df_prog['OD_low'])
            else:
                prog = GrowDiluteProgram(user=df_prog['user'],
                                         campaign=df_prog['campaign'],
                                         short_name=df_prog['short_name'],
                                         description=df_prog['description'],
                                         reactor_id=df_prog['reactor_id'],
                                         preculture=preculture,
                                         media=df_prog['media'],
                                         OD_low=df_prog['OD_low'],
                                         OD_high=df_prog['OD_high'])

    def add_program(self, program):
        if program.reactor_id in self.programs:
            raise Exception('Reactor %i already used by a program !')
        self.programs[program.reactor_id] = program
        program._manager = self.manager

    def remove_program(self, rid):
        del self.programs[rid]

    def schedule_sampling(self, time_sampling):
        if time_sampling > time() - 1:
            self.sampling_schedule.append(time_sampling)
        else:
            print('Cannot schedule the sampling: invalid sampling time (in the past)')

    def start(self):
        if self.status != 'running':
            print('Starting !')
            self.run_thread = SessionRun(self)
            self.run_thread.start()

    def pause(self):
        logging.info('T= {}| Pausing the session'.format(time()))
        print('Pausing session')
        if self.status == 'running':
            self.status = 'pausing'
            if self.run_thread is None:
                logging.error('T= {}| Trying to pause a session with no more run thread.'.format(time()))
            else:
                if not self.run_thread.is_alive():
                    logging.error('T= {}| Trying to pause a session with an inactive run thread ??'.format(time()))
            while self.status != 'paused':
                sleep(0.01)
        else:
            logging.info('T= {}| Trying to pause a session that is not running (status is {})'.format(time(),self.status))
            print('Not running ? status is {}.'.format(self.status))

    def stop(self):
        print('Stopping session')
        if self.status == 'running':
            self.status = 'stopping'
            while self.status != 'stopped':
                sleep(0.01)
        else:
            self.status = 'stopping'
            self.manager.shut_off_all_pumps_and_valves()
            self.status = 'stopped'
        # last data writing
        for program in self.programs.values():
            dfs = program.all_data_to_dfs()
            for data_type,df in dfs.items():
                if df is not None:
                    df.to_csv(program.output_path + '/' + data_type + '.csv', index=False)
                    # atlas
                    try:
                        df.to_csv(program.atlas_path + '/' + data_type + '.csv', index=False)
                    except BaseException as e:
                        logging.error('Something went wrong in atlas backup (error: {})'.format(e))

    def reset_all_sampling_state(self):
        self.ot2_remaining_tips_by_8 = 24
        self.ot2_state_sampling_plates = {1:[13], 2:[13]}
        self.ot2_state_reservoir = {c:0 for c in range(1,13)}
        self.manager.send_ot2_request('tips/reset')
        return 'Done sampling state reset... Means should be two clean sampling plates, full reservoir, two full boxes of tips.'

    def next_sampling_plate_and_col(self):
        # find plate first
        i_plate = 1
        # we should distinguish wether we can use col 1 or not
        if min(self.programs_with_cytometry.keys()) > 8:
            while min(self.ot2_state_sampling_plates[i_plate]) == 2: # because cannot use col 1
                i_plate += 1
        else:
            while min(self.ot2_state_sampling_plates[i_plate]) == 1: # means we used the whole plate
                i_plate += 1
        # find col
        col = min(self.ot2_state_sampling_plates[i_plate]) - 1
        # for now,  if only 9-16 reactors, we cannot use the col 1...
        if col == 1 and min(self.programs_with_cytometry.keys()) > 8:
            col = 12
            i_plate += 1
        return i_plate,col

    def loop(self):
        try:
            ### beginning of the loop !
            loop_start_time = time()

            ## blink once
            try:
                self.manager.send_ot2_request('blink/1')
            except:
                pass

            ## OD measurements
            # first stop the bubbling because bubbles can affect OD
            bubbling_set_ids = set([reactors[reactor_id].bubbling_set_id for reactor_id in self.programs])
            for bubbling_set_id in bubbling_set_ids:
                self.manager.set_valve_state(bubbling_set_id, True)
            if self.status != 'running':
                return
            # OD measurements per se
            reactor_set_ids = set([reactors[reactor_id].reactor_set_id for reactor_id in self.programs])
            OD_measuring_thread = ODMeasurement(self.manager, reactor_set_ids)
            OD_measuring_thread.start()
            while OD_measuring_thread.is_alive():
                if self.status != 'running':
                    sleep(0.1)
                    return
            ODs = OD_measuring_thread.ODs
            # provide the ODs to the programs
            for reactor_id in self.programs:
                self.programs[reactor_id].receive_OD_reading(ODs[reactor_id])
            # put back bubbling
            for bubbling_set_id in bubbling_set_ids:
                self.manager.set_valve_state(bubbling_set_id, False)

            ## drain via output
            # set output valve to the rest position, for which output pump takes from drain tube (should be useless because already off)
            for reactor_id,program in self.programs.items():
                self.manager.set_valve_state(reactors[reactor_id].out_valve_id, False)
            # position the robot to the bin so that we can collect the excess culture volume in the bioharzard tank
            self.manager.send_ot2_request('gotobin')
            # open the output pumps for the duration
            drain_threads = []
            for reactor_id,program in self.programs.items():
                if program.drain_out_pump_duration_s > 0:
                    drain_threads.append(self.manager.open_pump_for_duration(reactors[reactor_id].pump_out_slot_id, program.drain_out_pump_duration_s))
            if not self.virtual_mode:
                for drain_thread in drain_threads:
                    drain_thread.join()

            ## blink twice
            try:
                self.manager.send_ot2_request('blink/2')
            except:
                pass

            ## media addition
            dilution_threads = []
            for reactor_id,program in self.programs.items():
                dilution_duration = program.compute_dilution_duration()
                if dilution_duration:
                    program.drain_out_pump_duration_s = 2. * dilution_duration
                    dilution_threads.append(
                        self.manager.open_pump_for_duration(reactors[reactor_id].pump_in_slot_id, dilution_duration))
                else:
                    program.drain_out_pump_duration_s = 0.
            for dilution_thread in dilution_threads:
                dilution_thread.join()

            ## blink thrice
            try:
                self.manager.send_ot2_request('blink/3')
            except:
                pass

            ## are we waiting for data ? if arrives, we should wash
            if self.waiting_for_guava_data:
                # gather the fcs files on atlas
                fcs3_files = guava2data.gather_fcs3_files_on_atlas('Z:')
                fcs3_files.sort(key=lambda x: x[-1])
                if self.samplings and self.samplings[-1] > 0:
                    sampling_datetime = datetime.fromtimestamp(self.samplings[-1])
                    fcs3_files_after_sampling = [(root,f,d) for root,f,d in fcs3_files if d > sampling_datetime]
                    if fcs3_files_after_sampling: # there is data !
                            root,f,d = fcs3_files_after_sampling[0]
                            if d < datetime.today(): # to be sure
                                self.waiting_for_guava_data = False
                                self.should_wash = True
                                file_to_fetch = os.path.join(root, f)
                                sleep(5) # in case the file is being written ?
                                logging.info('T= {}| Found the fcs3 file {} to be the one corresponding to last timepoint'.format(time(), file_to_fetch))
                                fcs3 = guava2data.GuavaFCS3(file_to_fetch, only_meta=True)
                                for reactor_id,program in self.programs_with_cytometry.items():
                                    if program.samplings['guava_fcs_file']:
                                        program.samplings['guava_fcs_file'][-1] = file_to_fetch
                                        guava_well = reactors[reactor_id].guava_row + reactors[reactor_id].guava_col
                                        if guava_well not in fcs3.samples:
                                            logging.info('T= {}| cannot find data for well {} ?!'.format(time(),guava_well))
                                            continue
                                        logging.info('T= {}| loading data for well {}'.format(time(),guava_well))
                                        sample = fcs3.get_sample(guava_well)
                                        sample_df = sample.to_df(min_n_return_none=2)
                                        logging.info('T= {}| finished loading data for well {}'.format(time(),guava_well))
                                        # only take sample if expected number of events
                                        if sample_df is not None and len(sample_df) == self.guava_num_events:
                                            logging.info('T= {}| data for well {} has good number of events'.format(time(),guava_well))
                                            # add the sampling_time_s as a columns
                                            sample_df['sampling_time_s'] = program.samplings['time_s'][-1]
                                            if program.cells is None:
                                                program.cells = sample_df
                                            else:
                                                program.cells = pd.concat([program.cells, sample_df], ignore_index=True)
                                            program.cells.to_csv(f'{program.output_path}/cells.csv', index=False)
                                            program.cells.to_csv(f'{program.atlas_path}/cells.csv', index=False)
                                        else:
                                            if sample_df is not None:
                                                logging.info('T= {}| not enough events ! only {}'.format(time(),len(sample_df)))
                                            else:
                                                logging.info('T= {}| less than 2 events !!!!!!!'.format(time()))
                            else:
                                logging.info('T= {}| Weird error / problem: cytometry data from the future ? I assume its not the good one'.format(time()))

            ## should we wash the guava plate ?
            if self.should_wash:
                self.manager.send_guava_request('toggle_tray')
                # find wich cyto col should be washed
                cyto_cols_to_wash = []
                if min(self.programs_with_cytometry.keys()) < 9:
                    cyto_cols_to_wash.append(2)
                if max(self.programs_with_cytometry.keys()) > 8:
                    cyto_cols_to_wash.append(1)
                # iterate on cyto cols to wash
                for cyto_col_to_wash in cyto_cols_to_wash:
                    # find which col of the reservoir to use
                    reservoir_col = 12
                    while self.ot2_state_reservoir[reservoir_col] >= 2:
                        reservoir_col -= 1
                    if reservoir_col < 1:
                        raise Exception('cannot wash anymore ?!!?')
                    # do the wash
                    for i in range(3):
                        self.manager.send_ot2_request(f'wash_from_trough/{reservoir_col}/{cyto_col_to_wash}')
                    # update reservoir data
                    self.ot2_state_reservoir[reservoir_col] += 1
                self.manager.send_ot2_request('tips/drop')
                self.ot2_remaining_tips_by_8 -= 1
                self.manager.send_guava_request('toggle_tray')
                self.should_wash = False

            ## should we sample to the guava and start acquire ?
            if not self.waiting_for_guava_data and self.ot2_remaining_tips_by_8 > 0:
                # check if we passed the time of the next sampling
                if self.sampling_schedule:
                    next_sampling_time = sorted(self.sampling_schedule)[0]
                    if next_sampling_time < time():
                        logging.info('T= {}| We passed time of a scheduled sampling'.format(time()))
                        # check that we can do the sampling
                        i_plate,ref_sampling_col = self.next_sampling_plate_and_col()
                        if i_plate not in [1,2]:
                            logging.error('T= {}| Asked for a cyto sampling on a wrong plate ! (#{}). Forgot to reset num of samplings ?'.format(time(), i_plate))
                        else:
                            # decide programs that we can sample based on volume estimate
                            programs_to_sample = {}
                            for reactor_id, program in self.programs_with_cytometry.items():
                                vol_estim = program.volume_estimates
                                if True:#vol_estim.empty or vol_estim['volume_estimation_mL'].iloc[-1] > 15:
                                    programs_to_sample[reactor_id] = program
                            # start the sampling !
                            real_sampling_time = time()
                            self.sampling_schedule.remove(next_sampling_time)
                            self.samplings.append(real_sampling_time)
                            # record the sampling into each program so that can be output
                            for reactor_id, program in programs_to_sample.items():
                                program.samplings['time_s'].append(real_sampling_time)
                                program.samplings['guava_well'].append(reactors[reactor_id].guava_row + reactors[reactor_id].guava_col)
                                program.samplings['guava_fcs_file'].append('not_fetched_yet')
                                if max(self.programs_with_cytometry.keys()) < 9 or min(self.programs_with_cytometry.keys()) > 8:
                                    program.samplings['sampling_plate_and_col'].append(f'plate-{i_plate}_col-{ref_sampling_col}')
                                else:
                                    if reactor_id < 9:
                                        program.samplings['sampling_plate_and_col'].append(f'plate-{i_plate}_col-{ref_sampling_col-1}')
                                    else:
                                        program.samplings['sampling_plate_and_col'].append(f'plate-{i_plate}_col-{ref_sampling_col}')
                            # move the sampling output above the thrash
                            self.manager.send_ot2_request('gotobin')
                            # set the out valve to sampling
                            for reactor_id,program in programs_to_sample.items():
                                self.manager.set_valve_state(reactors[reactor_id].out_valve_id, True)
                            # gather the out slots and dead volumes
                            pump_out_slot_ids = [reactors[reactor_id].pump_out_slot_id for reactor_id in programs_to_sample]
                            dead_volumes = [program.dead_volume_sampling_line_mL for program in programs_to_sample.values()]
                            # dead volume first
                            pump_threads = []
                            for pump_out_slot_id,dead_volume in zip(pump_out_slot_ids, dead_volumes):
                                pump_threads.append(
                                    self.manager.open_pump_for_volume(pump_out_slot_id, dead_volume))
                            if not self.virtual_mode:
                                for pump_thread in pump_threads:
                                    pump_thread.join()
                            # do volume renewal if activated
                            pump_in_slot_ids = [reactors[reactor_id].pump_in_slot_id for reactor_id,program in programs_to_sample.items() if program.renew_sampled_volume]
                            dead_volumes = [program.dead_volume_sampling_line_mL for program in programs_to_sample.values() if program.renew_sampled_volume]
                            pump_threads = []
                            for pump_in_slot_id,dead_volume in zip(pump_in_slot_ids, dead_volumes):
                                pump_threads.append(
                                    self.manager.open_pump_for_volume(pump_in_slot_id, dead_volume))
                            # declare the dilutions to the program
                            for reactor_id,program in programs_to_sample.items():
                                if program.renew_sampled_volume:
                                    program.dilutions['time_s'].append(time())
                                    program.dilutions['duration_s'].append(program.dead_volume_sampling_line_mL / self.manager.input_pump_flow_rate_mL_min(reactor_id) * 60.)
                                    program.dilutions['est_flow_rate_uL_per_s'].append(self.manager.input_pump_flow_rate_mL_min(reactor_id)/60.*1000.)
                            if not self.virtual_mode:
                                for pump_thread in pump_threads:
                                    pump_thread.join()
                            # shake the ot2 arm to remove ready to fall drops
                            self.manager.send_ot2_request('shakeit')
                            # compute what move we should do for sampling positioning
                            plate_labware_name = {1:'sampling_metal', 2:'sampling_metal_2'}[i_plate]
                            # if not only 1-8 we have the positioning is not the sampling
                            col_sampling_positioning = ref_sampling_col
                            if max(self.programs_with_cytometry.keys()) > 8:
                                col_sampling_positioning -= 1
                            if col_sampling_positioning > 12 or col_sampling_positioning < 1:
                                logging.error('T= {}| Invalid col sampling pos ? ! ({}). Forgot to reset num of samplings ?'.format(time(),col_sampling_positioning))
                            else:
                                # move into sampling position (first cytoplate for path optim reason)
                                self.manager.send_ot2_request(f'move_to/cytoplate/1')
                                self.manager.send_ot2_request(f'move_to/{plate_labware_name}/{col_sampling_positioning}')
                                # sampling !
                                pump_threads = []
                                for pump_out_slot_id in pump_out_slot_ids:
                                    pump_threads.append(
                                        self.manager.open_pump_for_duration(pump_out_slot_id, 0.75))
                                if not self.virtual_mode:
                                    for pump_thread in pump_threads:
                                        pump_thread.join()
                                # update the state of sampling plate
                                self.ot2_state_sampling_plates[i_plate].append(ref_sampling_col)
                                if min(self.programs_with_cytometry.keys()) <= 8 and max(self.programs_with_cytometry.keys()) > 8:
                                    self.ot2_state_sampling_plates[i_plate].append(ref_sampling_col-1)
                                # prepare the guava (will make the worklist, and open the tray)
                                guava_wells = [reactors[reactor_id].guava_row + reactors[reactor_id].guava_col for reactor_id in programs_to_sample]
                                guava_wells_str = '_'.join(guava_wells)
                                self.manager.send_guava_request(f'prepare/{guava_wells_str}/{self.guava_num_events}')
                                # move to cytoplate plate then to bin to minimize cross-cont while moving
                                self.manager.send_ot2_request(f'move_to/cytoplate/1')
                                self.manager.send_ot2_request('gotobin')
                                # fill guava plate
                                reservoir_col = 12
                                while self.ot2_state_reservoir[reservoir_col] >= 2:
                                    reservoir_col -= 1
                                if reservoir_col < 1:
                                    raise Exception('no more reservoir slots ?!!?')
                                plate_labware_name = {1:'sampling', 2:'sampling_2'}[i_plate]
                                if max(self.programs_with_cytometry.keys()) < 9: # 1-8
                                    self.manager.send_ot2_request(f'dilute/10/{plate_labware_name}/{reservoir_col}/{ref_sampling_col}/2')
                                    if self.additional_cytoplate_dilution:
                                        self.manager.send_ot2_request(f'dilute_cytoplate/20/{reservoir_col}/2')
                                elif min(self.programs_with_cytometry.keys()) > 8: # 9-16
                                    self.manager.send_ot2_request(f'dilute/10/{plate_labware_name}/{reservoir_col}/{ref_sampling_col}/1')
                                    if self.additional_cytoplate_dilution:
                                        self.manager.send_ot2_request(f'dilute_cytoplate/20/{reservoir_col}/1')
                                else: # both
                                    self.manager.send_ot2_request(f'dilute/10/{plate_labware_name}/{reservoir_col}/{ref_sampling_col-1}/2')
                                    if self.additional_cytoplate_dilution:
                                        self.manager.send_ot2_request(f'dilute_cytoplate/20/{reservoir_col}/2')
                                    self.manager.send_ot2_request(f'dilute/10/{plate_labware_name}/{reservoir_col}/{ref_sampling_col}/1')
                                    if self.additional_cytoplate_dilution:
                                        self.manager.send_ot2_request(f'dilute_cytoplate/10/{reservoir_col}/1')
                                # start guava acquisition (careful !!!!!!!! on the background, we know its finished only when fcs data appears on atlas)
                                self.manager.send_guava_request('acquire')
                                # store that the sampling has been done, so one should wait for guava data
                                self.waiting_for_guava_data = True

            ## output data
            for program in self.programs.values():
                dfs = program.all_data_to_dfs()
                for data_type,df in dfs.items():
                    if df is not None:
                        df.to_csv(program.output_path + '/' + data_type + '.csv', index=False)
                        # also write data to ATLAS
                        try:
                            df.to_csv(program.atlas_path + '/' + data_type + '.csv', index=False)
                        except BaseException as e:
                            logging.error('Something went wrong in atlas backup (error: {})'.format(e))

            ## update LED change data
            led_histories = self.manager.get_all_leds_history()
            for reactor_id,program in self.programs.items():
                program.update_LED_data(led_histories[reactor_id])

            ## check and apply events
            logging.info('T= {}| Starting to check events'.format(time()))
            for program in self.programs.values():
                logging.info('T= {}| Starting to check events of reactor id {}'.format(time(), program.reactor_id))
                for event in program.events:
                    event.check_and_apply(program)

            ## wait until next loop
            self.loop_counter += 1
            while time() - loop_start_time < self.loop_duration_s:
                if self.status != 'running':
                    return
                sleep(0.1)

        except BaseException as e:
            message = f'Something went wrong (error: {e}). Main loop continues, BE CAREFUL / CHECK'
            logging.error(message)
            self.manager.send_discord(message)
            while time() - loop_start_time < self.loop_duration_s:
                sleep(0.1)

    ## useful functions to use outside of main loop

    def prime_input_pumps(self,open_duration_s=6):
        try:
            reactor_ids = self.programs.keys()
            logging.info('T= %f| Started an input pumps priming sequence' % time())
            while input('Priming of pumps with liquid. Make sure pump output is setup as desired. Enter y to flow for %f seconds, otherwise ok\n'%open_duration_s) == 'y':
                pump_threads = []
                for reactor_id in reactor_ids:
                    pump_threads.append(
                        self.manager.open_pump_for_duration(reactors[reactor_id].pump_in_slot_id, open_duration_s))
                for pump_thread in pump_threads:
                    pump_thread.join()
        except BaseException as e:
            logging.error('Something went wrong, shutting off all coils. (error: %s)' % str(e))
            self.manager.shut_off_all_pumps_and_valves()

    def pinch_drain_valves(self, reactor_ids=None):
        try:
            if reactor_ids is None:
                reactor_ids = self.programs.keys()
            logging.info('T= %f| Pinching front valves' % time())
            for reactor_id in reactor_ids:
                    self.manager.set_valve_state(reactors[reactor_id].out_valve_id, True)
        except BaseException as e:
            logging.error('Something went wrong, shutting off all coils. (error: %s)' % str(e))
            self.manager.shut_off_all_pumps_and_valves()

    def unpinch_drain_valves(self, reactor_ids=None):
        try:
            if reactor_ids is None:
                reactor_ids = self.programs.keys()
            logging.info('T= %f| Pinching front valves' % time())
            for reactor_id in reactor_ids:
                    self.manager.set_valve_state(reactors[reactor_id].out_valve_id, False)
        except BaseException as e:
            logging.error('Something went wrong, shutting off all coils. (error: %s)' % str(e))
            self.manager.shut_off_all_pumps_and_valves()

    def add_media(self,volume_mL, cleaning=False):
        try:
            if cleaning:
                reactor_ids = self.reactor_ids_to_clean
            else:
                reactor_ids = self.programs.keys()
            logging.info('T= %f| Started an input pumps sequence' % time())
            pump_threads = []
            for reactor_id in reactor_ids:
                pump_threads.append(
                    self.manager.open_pump_for_volume(reactors[reactor_id].pump_in_slot_id, volume_mL))
            for pump_thread in pump_threads:
                pump_thread.join()
        except BaseException as e:
            logging.error('Something went wrong, shutting off all coils. (error: %s)' % str(e))
            self.manager.shut_off_all_pumps_and_valves()

    def remove_volume(self, volume_mL, cleaning=False, pinch_drain_valves=True):
        try:
            self.manager.send_ot2_request('gotobin')
            if cleaning:
                reactor_ids = self.reactor_ids_to_clean
            else:
                reactor_ids = self.programs.keys()
            if pinch_drain_valves:
                self.pinch_drain_valves(reactor_ids)
            pump_threads = []
            for reactor_id in reactor_ids:
                pump_threads.append(
                    self.manager.open_pump_for_volume(reactors[reactor_id].pump_out_slot_id, volume_mL))
            for pump_thread in pump_threads:
                pump_thread.join()
            self.unpinch_drain_valves(reactor_ids)
        except BaseException as e:
            logging.error('Something went wrong, shutting off all coils. (error: %s)' % str(e))
            self.manager.shut_off_all_pumps_and_valves()

    # for cleaning
    def register_reactors_for_cleaning(self, reactor_ids):
        for reactor_id in reactor_ids:
            if reactor_id in self.programs:
                del self.programs[reactor_id]
        self.reactor_ids_to_clean = reactor_ids

    def start_input_flow(self, reactor_ids=None):
        if reactor_ids is None:
            reactor_ids = self.programs.keys()
        for reactor in [reactors[id] for id in reactor_ids]:
            self.manager.set_pump_state(reactor.pump_in_slot_id, True)

    def stop_input_flow(self, reactor_ids=None):
        if reactor_ids is None:
            reactor_ids = self.programs.keys()
        for reactor in [reactors[id] for id in reactor_ids]:
            self.manager.set_pump_state(reactor.pump_in_slot_id, False)

    def start_output_flow(self, reactor_ids=None):
        self.manager.send_ot2_request('gotobin')
        if reactor_ids is None:
            reactor_ids = self.programs.keys()
        for reactor in [reactors[id] for id in reactor_ids]:
            self.manager.set_pump_state(reactor.pump_out_slot_id, True)

    def stop_output_flow(self, reactor_ids=None):
        if reactor_ids is None:
            reactor_ids = self.programs.keys()
        for reactor in [reactors[id] for id in reactor_ids]:
            self.manager.set_pump_state(reactor.pump_out_slot_id, False)

    def start_output_pump_calibration(self):
        duration_seconds = {}
        for rid in self.reactor_ids_to_clean:
            self.manager.set_pump_state(f'H{rid}', True)
            tstart = time()
            input()
            tend = time()
            self.manager.set_pump_state(f'H{rid}', False)
            duration_seconds[rid] = tend - tstart
            sleep(3)
        pump_flow_rates = {rid:9.5/T*60 for rid,T in duration_seconds.items()}
        for rid,flow_rate in pump_flow_rates.items():
            pump_id = self.manager.output_pumps[rid]
            old_flow_rate = self.manager.pump_calib_pars.loc[pump_id, 'flow_rate_mL_min']
            print(f'For pump {pump_id} (installed as output for reactor {rid}), new flow rate: {flow_rate} mL/min (was {old_flow_rate})')
            logging.info(f'For pump {pump_id} (installed as output for reactor {rid}), new flow rate: {flow_rate} mL/min (was {old_flow_rate})')
            self.manager.pump_calib_pars.loc[pump_id, 'flow_rate_mL_min'] = flow_rate
            self.manager.pump_calib_pars.loc[pump_id, 'date'] = date.today()
        self.manager.pump_calib_pars.to_csv('cytoreactors/calibration/pump-calibration-parameters.csv')

    def start_input_pump_calibration(self):
        duration_seconds = {}
        for rid in self.reactor_ids_to_clean:
            self.manager.set_pump_state(f'L{rid}', True)
            tstart = time()
            input()
            tend = time()
            self.manager.set_pump_state(f'L{rid}', False)
            duration_seconds[rid] = tend - tstart
            sleep(3)
        pump_flow_rates = {rid:9.5/T*60 for rid,T in duration_seconds.items()}
        for rid,flow_rate in pump_flow_rates.items():
            pump_id = self.manager.input_pumps[rid]
            old_flow_rate = self.manager.pump_calib_pars.loc[pump_id, 'flow_rate_mL_min']
            print(f'For pump {pump_id} (installed as input for reactor {rid}), new flow rate: {flow_rate} mL/min (was {old_flow_rate})')
            logging.info(f'For pump {pump_id} (installed as input for reactor {rid}), new flow rate: {flow_rate} mL/min (was {old_flow_rate})')
            self.manager.pump_calib_pars.loc[pump_id, 'flow_rate_mL_min'] = flow_rate
            self.manager.pump_calib_pars.loc[pump_id, 'date'] = date.today()
        self.manager.pump_calib_pars.to_csv('cytoreactors/calibration/pump-calibration-parameters.csv')



    ## plotting utilities

    def plot_reactor_group(self, reactor_ids=None, xlims=None, show_fl=False, ch_fl='GRN-B-HLin'):
        if reactor_ids is None:
            reactor_ids = list(self.programs.keys())
        if show_fl:
            ax_OD, ax_gr, ax_LEDs, ax_fl = self.programs[reactor_ids[0]].plot(show_fl=show_fl, ch_fl=ch_fl)
        else:
            ax_OD, ax_gr, ax_LEDs = self.programs[reactor_ids[0]].plot(show_fl=show_fl, ch_fl=ch_fl)
        for rid in reactor_ids[1:]:
            if show_fl:
                self.programs[rid].plot(axs=(ax_OD, ax_gr, ax_LEDs, ax_fl), show_fl=show_fl, ch_fl=ch_fl)
            else:
                self.programs[rid].plot(axs=(ax_OD, ax_gr, ax_LEDs), show_fl=show_fl, ch_fl=ch_fl)
        if show_fl:
            for ax in [ax_OD,ax_gr,ax_LEDs,ax_fl]:
                ax.grid()
                if xlims:
                    ax.set_xlim(xlims)
            return ax_OD, ax_gr, ax_LEDs, ax_fl
        for ax in [ax_OD, ax_gr, ax_LEDs]:
            ax.grid()
            if xlims:
                ax.set_xlim(xlims)
        return ax_OD, ax_gr, ax_LEDs


# thread that terminates when session is not running
class SessionRun(Thread):
    def __init__(self, session):
        super().__init__()
        self.session = session
    def run(self):
        if self.session.status != 'running':
            self.session.status = 'running'
            while self.session.status == 'running':
                self.session.loop()
        if self.session.status == 'pausing':
            self.session.manager.shut_off_all_pumps_and_valves()
            self.session.status = 'paused'
            while self.session.status == 'paused':
                sleep(0.1)
        elif self.session.status == 'stopping':
            self.session.manager.shut_off_all_pumps_and_valves()
            self.session.status = 'stopped'
        else:
            logging.error('T= {}| Unknown session status {}'.format(time(),self.session.status))
            print('unknown session status {}'.format(self.session.status))

# thread for OD measurements
class ODMeasurement(Thread):
    def __init__(self, manager, reactor_set_ids):
        super().__init__()
        self.manager = manager
        self.reactor_set_ids = reactor_set_ids
        self.ODs = None
    def run(self):
        self.ODs = self.manager.get_turbidity_readings(self.reactor_set_ids)
