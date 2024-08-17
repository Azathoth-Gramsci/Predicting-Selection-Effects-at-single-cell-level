clear
close all


%% DATA
% 2mi_pulses2 = 5 min pulse every 4 hours and a total of 5 pulses
% 2mi_pulses1 = 5 min pulse every hour and a total of 3 pulses
% 2mi_fulllight = continuous light
% 2mi_reppulses2 = 2 min pulse every 4 hours 
% 2mi_reppulses1 = 1 min pulses every hour 
% replace these names in the lines below to change the experiment
load('..\..\Matlab_Data\2mi_reppulses1_data') % here

timeScale = 60;
X = table2array(readtable('..\..\Matlab_Data\2mi_reppulses1_LightProfile.csv')); % and here
u = zeros(1,floor(X(end,2)*timeScale)+1);
%note: rounding and discretizing time in minutes
for i = 1:size(X,1)
    indStart = round(X(i,1)*timeScale)+1;
    indEnd = round(X(i,2)*timeScale);
    u(indStart:indEnd) = X(i,3);
end
Tmax = length(u);
u(Tmax+1) = u(Tmax);


DataAna

fracRec = numCells1./(numCells0+numCells1);

%% MODEL
tDiv = 103;
lambda = log(2)/tDiv;

%FIX From integrated
meanX = 10;
scale = 5.0617;
CVTF = 0.3525; 

%zero plasmid death rate
pKill = lambda/2.2;

%plasmid paras
bP = lambda;
aP = 0.85*bP;

%transcription factor
b = lambda;

%bursts
shapeTF = 1/CVTF^2;
scaleTF = meanX/shapeTF;
burst = scaleTF;
a = shapeTF*b;

%parameters of Hill function that models cell differentiation
rH = 1;
Kh = 132/scale;
nH = 4;
uIn = 0.1*rH;

%combine
paras = [a,b,burst,aP,bP,pKill,rH,Kh,nH];
parasDiffFunction = [rH,Kh,nH];

% %finite state projection size (computationally expensive)
% truncX = 250;
% truncP = 30;

% to test if the code runs
truncX = 60;
truncP = 20;

% % to test if the code runs
% truncX = 60;
% truncP = 10;

%with zero light input
[Au0] = MatrixBuilder(truncX,truncP,paras,0); 
[Au0NoDiff] = MatrixBuilderNoDiff(truncX,truncP,paras); %can calculate with a model without differentiation to determine the initial distribution

%initial distribution
[V,D] = eig(Au0NoDiff);
indMax = find(abs(diag(D))<=min(abs(diag(D))),1,'first');
pQSD = V(:,indMax)./sum(V(:,indMax));
PZQSD = reshape(pQSD,truncX+1,truncP+1);
PXQSD = sum(PZQSD,2);
PPQSD = sum(PZQSD,1)';
meanPQSD = [0:1:truncP]*PPQSD;

initdiff = 0.06; %initial fraction of differentiated cells
p0 = [pQSD*(1-initdiff);pQSD*initdiff]; %setting initial conditions for both diff and undiff cells to quasi stationary

%with active light input
[Au1] = MatrixBuilder(truncX,truncP,paras,uIn);

%solve master equation using sparse matrix techniques
ts = 1;
Au1s = sparse(Au1);
Au0s = sparse(Au0);
pOut = zeros(size(Au1,1),Tmax+1);
pOutNorm = zeros(size(Au1,1),Tmax+1);
pt1 = p0;
pOut(:,1) = p0;
pOutNorm(:,1) = p0;
for i = 1:ts:Tmax
    if u(i) == 0
        Aus = Au0s;
    else
        Aus = Au1s;
    end
    [pt1, err, hump] = expv( ts, Aus, pt1);
    pOut(:,i+1) = pt1;
    pOutNorm(:,i+1) = pOut(:,i+1)./sum(pOut(:,i+1));
end
Yout_FSP = pOutNorm;


%% analyze and plot

steps = 1;
MeanX = zeros(1,Tmax/steps);
MeanP = zeros(1,Tmax/steps);
MeanX2 = zeros(1,Tmax/steps);
MeanP2 = zeros(1,Tmax/steps);
MedianX = zeros(1,Tmax/steps);
pGrowVec = zeros(1,Tmax/steps);
pDiffVec = zeros(1,Tmax/steps);

k = 1;
for i = 1:steps:Tmax+1
    pGrow = sum(Yout_FSP(1:(truncX+1)*(truncP+1),i));
    pDiff = 1-sum(Yout_FSP(1:(truncX+1)*(truncP+1),i));
    pGrowVec(k) = pGrow;
    pDiffVec(k) = pDiff;

    PZ = reshape(Yout_FSP(1:(truncX+1)*(truncP+1),i)/pGrow,truncX+1,truncP+1);
    PX = sum(PZ,2);
    PP = sum(PZ,1)';
    meanX = sum([0:1:truncX]'.*PX);
    meanP = sum([0:1:truncP]'.*PP);
    MeanX(k) = meanX;
    MeanP(k) = meanP;
    
    FX = cumsum(PX);
    indX = find(FX > 0.5, 1, 'first');
    %continuous median
    medianX = indX-1;
    diff1 = FX(indX)-0.5;
    if indX > 1
        diff2 = 0.5-FX(indX-1);
        weight1 = diff2/(diff1+diff2);
        weight2 = diff1/(diff1+diff2);
        medianX2 = medianX*weight1 + (medianX-1)*weight2;
    else
        medianX2 = indX-1;
    end
    MedianX(k) = medianX2;

    %differentiated
    PZ2 = reshape(Yout_FSP((truncX+1)*(truncP+1)+1:end,i)/pDiff,truncX+1,truncP+1);
    PX2 = sum(PZ2,2);
    PP2 = sum(PZ2,1)';
    meanX2 = sum([0:1:truncX]'.*PX2);
    meanP2 = sum([0:1:truncP]'.*PP2);
    MeanX2(k) = meanX2;
    MeanP2(k) = meanP2;

    k = k + 1;
end


%% PLOTTING

color1 = 'b';
color2 = 'r';
color3 = 'g';
color4 = 'k';

figure(1)
hold on
plot([0:steps:Tmax]./timeScale,pDiffVec./(pGrowVec+pDiffVec),color2,'Linewidth',2)
hold on
plot(tvec,fracRec,'bx','Linewidth',1.5)
xlabel('Time')
ylabel('Differentiated Fraction')
axis([0 65 0 1])
hold on
plot([0:1:Tmax]./timeScale,u)


figure(2)
hold on
plot([0:steps:Tmax]./timeScale,MeanX*scale,color2,'Linewidth',2)
hold on
plot([0:steps:Tmax]./timeScale,MeanX2*scale,color4,'Linewidth',2)
hold on
plot(tvec,meanTF0,color2,'Linewidth',1)
hold on
plot(tvec,meanTF1,color4,'Linewidth',1)
xlabel('Time')
ylabel('Mean X')
%axis([0 50 95 420])

figure(3)
hold on
plot([0:steps:Tmax]./timeScale,MedianX*scale,color1,'Linewidth',2)
hold on
plot(tvec,medTF0,color2,'Linewidth',2)
xlabel('Time')
ylabel('Median X')
axis([0 65 0 200])

figure(4)
hold on
plot([0:steps:Tmax]./timeScale,MeanP,color1,'Linewidth',2)
hold on
plot([0:steps:Tmax]./timeScale,MeanP2,color2,'Linewidth',2)
xlabel('Time')
ylabel('Mean P')
axis([0 65 1.5 6.5])

%normal AF for convolution with error
VarAF = 1.5;
ymin = -100;
ymax = 100;
ysteps = 1;
y = ymin:ysteps:ymax;
muY = 0;
sigmaY = sqrt(VarAF);
distY = normpdf(y,muY,sigmaY)*ysteps;

PXQSDconv = conv(PXQSD,distY);
z = 0+ymin:1:truncX+ymax;

m = 1;
for i = 1:1:length(tvec)-1
    distTF0 = histc(cellsTF0{i},edgesTF)./length(cellsTF0{i});
    distTF1 = histc(cellsTF1{i},edgesTF)./length(cellsTF1{i});

    testTime = max(round(tvec(i)*timeScale),0)+1;
    pGrow = sum(Yout_FSP(1:(truncX+1)*(truncP+1),testTime));
    pDiff = 1-sum(Yout_FSP(1:(truncX+1)*(truncP+1),testTime)); 
    PZ = reshape(Yout_FSP(1:(truncX+1)*(truncP+1),testTime)/pGrow,truncX+1,truncP+1);
    PX = sum(PZ,2);
    PXconv = conv(PX,distY);
    PP = sum(PZ,1)';
    PZ2 = reshape(Yout_FSP((truncX+1)*(truncP+1)+1:end,testTime)/pDiff,truncX+1,truncP+1);
    PX2 = sum(PZ2,2);
    PX2conv = conv(PX2,distY);
    PP2 = sum(PZ2,1)';

    titleI = num2str(tvec(i));
    
    %plotting EL222 distribution at all time points
    if m < 73
        figure(10)
        hold on
        subplot(8,9,m)
        title(titleI)
        hold on
        plot(edgesTF,distTF0/stepsTF,color2,'Linewidth',2)
        hold on
        plot(edgesTF,distTF1/stepsTF,color3,'Linewidth',2)
        hold on
        plot(z*scale,PXconv/scale,color1,'Linewidth',2)
        hold on
        plot(z*scale,PX2conv/scale,color4,'Linewidth',2)
        hold on
        plot(z*scale,PXQSDconv/scale,'m','Linewidth',1)
        axis([-20 400 0 0.01])
    end
    
    %plotting EL222 distribution at chosen time points
    if i == 6
        figure(11)
        title(titleI)
    elseif i == 13
        figure(12)
        title(titleI)
    elseif i == 19
        figure(13)
        title(titleI)
 	elseif i == 41
        figure(14)
        title(titleI)
    elseif i == 30
        figure(15)
        title(titleI)
    elseif i == 32
        figure(16)
        title(titleI)
    elseif i == 38
        figure(17)
        title(titleI)
    elseif i == 52
        figure(18)
        title(titleI)
    end
    if i == 6 || i == 13 || i == 19 || i == 41 || i == 30 || i == 32 ||  i == 38 || i == 52
        hold on
        plot(edgesTF,distTF0/stepsTF,color2,'Linewidth',2)
        hold on
        plot(edgesTF,distTF1/stepsTF,color3,'Linewidth',2)
        hold on
        plot(z*scale,PXconv/scale,color1,'Linewidth',2)
        hold on   
        plot(z*scale,PX2conv/scale,color4,'Linewidth',2)
        hold on
        plot(z*scale,PXQSDconv/scale,'m','Linewidth',1)
        axis([-20 700 0 0.02])
    end
    
    %plotting plasmid distribution at chosen time points
    if i == 16
        figure(21)
        title(titleI)
    elseif i == 17
        figure(22)
        title(titleI)
    elseif i == 23
        figure(23)
        title(titleI)
    elseif i == 24
        figure(24)
        title(titleI)
    elseif i == 29
        figure(25)
        title(titleI)
    elseif i == 31
        figure(26)
        title(titleI)
    elseif i == 37
        figure(27)
        title(titleI)
    elseif i == 51
        figure(28)
        title(titleI)
    end
    if i == 16 || i == 17 || i == 23 || i == 24 || i == 29 || i == 31 ||  i == 37 || i == 51
        hold on
        plot([0:1:truncP],PP(1:truncP+1),color1,'Linewidth',2)
        hold on
        plot([0:1:truncP],PP2(1:truncP+1),color2,'Linewidth',2)
        hold on
        plot([0:1:truncP],PPQSD(1:truncP+1),'m','Linewidth',2)
        axis([0 20 0 0.6])
    end
    
    m = m+1;
    
end
