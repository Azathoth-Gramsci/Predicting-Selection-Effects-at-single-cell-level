clear
close all

%% DATA

% int_fulllight = continuous light
% int_pulses = 1h pulse every h for 5 pulses followed by 10h of rest then repeat the first 5 pulses

% replace these names in the lines below to change the experiment 
load('..\..\Matlab_Data\int_pulses_data') % here

timeScale = 60;
X = table2array(readtable('..\..\Matlab_Data\int_pulses_LightProfile.csv')); % and here

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
paras = [a,b,burst,rH,Kh,nH];
parasDiffFunction = [rH,Kh,nH];

%finite state projection size
% truncX = 150;

%for testing
truncX = 80;

%initial TF distribution
VarTF = meanX^2*CVTF^2;
pNB = 1 - meanX/(VarTF);
rNB = meanX^2/(VarTF-meanX);
x = 0:1:truncX;
p0 = nbinpdf(x,rNB,1-pNB);
diffFrac0 = 0;
p0u = p0*(1-diffFrac0);
p0d = p0*diffFrac0;


%with zero light inpus
[Au0] = MatrixBuilder(truncX,paras,0); 

%with active light input
[Au1] = MatrixBuilder(truncX,paras,uIn);

%solve master equation using sparse matrix techniques
ts = 1;
Au1s = sparse(Au1);
Au0s = sparse(Au0);
pOut = zeros(size(Au1,1),Tmax+1);
pOutNorm = zeros(size(Au1,1),Tmax+1);
pt1 = [p0u,p0d]';
pOut(:,1) = [p0u,p0d]';
pOutNorm(:,1) = [p0u,p0d]';
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

%analyze and plot
steps = 1;
MeanX = zeros(1,Tmax/steps);
MeanX2 = zeros(1,Tmax/steps);
MedianX = zeros(1,Tmax/steps);
pGrowVec = zeros(1,Tmax/steps);
pDiffVec = zeros(1,Tmax/steps);

k = 1;
for i = 1:steps:Tmax+1
    pGrow = sum(Yout_FSP(1:(truncX+1),i));
    pDiff = 1-sum(Yout_FSP(1:(truncX+1),i)); 
    pGrowVec(k) = pGrow;
    pDiffVec(k) = pDiff;
    PX = Yout_FSP(1:(truncX+1),i)/pGrow;
    if pDiff < 10^(-5)
        PX2 = PX;
    else
        PX2 = Yout_FSP(truncX+2:end,i)/pDiff;
    end
    meanX = sum([0:1:truncX]'.*PX);
    meanX2 = sum([0:1:truncX]'.*PX2);
    MeanX(k) = meanX;
    MeanX2(k) = meanX2;
    FX = cumsum(PX);
    indX = find(FX > 0.5, 1, 'first');
    medianX = indX-1;
    %continuous median
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
    
    k = k + 1;
end


%% PLOTTING

color1 = 'b';
color2 = 'r';
color3 = 'g';
color4 = 'k';

figure(1)
plot([0:steps:Tmax]./timeScale,pDiffVec./(pGrowVec+pDiffVec),color2,'Linewidth',2)
hold on
plot(tvec,fracRec,'bx','Linewidth',1.5)
xlabel('Time')
ylabel('Differentiated Frac')
axis([0 65 0 1])
hold on
plot([0:1:Tmax]./60,u)


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

figure(3)
hold on
plot([0:steps:Tmax]./timeScale,MedianX*scale,color1,'Linewidth',2)
hold on
plot(tvec,medTF0,color2,'Linewidth',2)
xlabel('Time')
ylabel('Median X')

%normal AF for convolution with error
VarAF = 1.5;
ymin = -100;
ymax = 100;
ysteps = 1;
y = ymin:ysteps:ymax;
muY = 0;
sigmaY = sqrt(VarAF);
distY = normpdf(y,muY,sigmaY)*ysteps;

PX0conv = conv(p0,distY);
z = 0+ymin:1:truncX+ymax;

m = 1;
for i = 1:1:length(tvec)-1
    distTF0 = histc(cellsTF0{i},edgesTF)./length(cellsTF0{i});
    distTF1 = histc(cellsTF1{i},edgesTF)./length(cellsTF1{i});
    
    testTime = max(round((tvec(i))*timeScale),0)+1;
    pGrow = sum(Yout_FSP(1:(truncX+1),testTime));
    pDiff = 1-sum(Yout_FSP(1:(truncX+1),testTime)); 
    PX = Yout_FSP(1:(truncX+1),testTime)/pGrow;
    PXconv = conv(PX,distY);
    PX2 = Yout_FSP(truncX+2:end,testTime)/pDiff;
    PX2conv = conv(PX2,distY);
 
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
        plot(z*scale,PX0conv/scale,'m','Linewidth',1)
        axis([0 120 0 0.025])
    end
    
    if i == 7
        figure(71)
        title(titleI)
    elseif i == 11
        figure(72)
        title(titleI)
    elseif i == 18
        figure(73)
        title(titleI)
 	elseif i == 22
        figure(74)
        title(titleI)
    elseif i == 20
        figure(75)
        title(titleI)
    elseif i == 16
        figure(76)
        title(titleI)
    elseif i == 32
        figure(77)
        title(titleI)
    elseif i == 35
        figure(78)
        title(titleI)
    end
    if i == 7 || i == 11 || i == 18 || i == 22 || i == 20 || i == 16 ||  i == 32 || i == 35
        hold on
        plot(edgesTF,distTF0/stepsTF,color2,'Linewidth',2)
        hold on
        plot(edgesTF,distTF1/stepsTF,color3,'Linewidth',2)
        hold on
        plot(z*scale,PXconv/scale,color1,'Linewidth',2)
        hold on   
        plot(z*scale,PX2conv/scale,color4,'Linewidth',2)
        hold on
        plot(z*scale,PX0conv/scale,'m','Linewidth',1)
        axis([0 120 0 0.03])
    end
    
    m = m+1;
    
end
