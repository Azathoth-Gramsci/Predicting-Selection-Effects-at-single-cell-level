cellsTF0 = cell(1,length(tvec));
cellsTF1 = cell(1,length(tvec));
cellsTF = cell(1,length(tvec));
cellsRec0 = cell(1,length(tvec));
cellsRec1 = cell(1,length(tvec));
cellsRec = cell(1,length(tvec));
numCells0 = zeros(1,length(tvec));
numCells1 = zeros(1,length(tvec));

cellsCerulean0 = cell(1,length(tvec));
cellsCerulean1 = cell(1,length(tvec));
cellsCeruleanNeon = cell(1,length(tvec));

meanTF0 = zeros(1,length(tvec));
meanTF1 = zeros(1,length(tvec));
meanRec0 = zeros(1,length(tvec));
meanRec1 = zeros(1,length(tvec));

medTF0 = zeros(1,length(tvec));
medTF1 = zeros(1,length(tvec));
medRec0 = zeros(1,length(tvec));
medRec1 = zeros(1,length(tvec));

varTF0 = zeros(1,length(tvec));
varTF1 = zeros(1,length(tvec));
varRec0 = zeros(1,length(tvec));
varRec1 = zeros(1,length(tvec));

cellsTFPool = [];
cellsRecPool = [];

for i = 1:length(tvec)
    
    %gating
    indclear = find(data_matrix(:,15,i)<0.2  | data_matrix(:,16,i)>0.5);
    data_matrix(indclear,22,i) = 2;
    
    indR = find(data_matrix(:,22,i) < 2);
    dataValidRaw = data_matrix(indR,:,i);

    %AF clean
    CeruleanAF = 10;
    NeonAF = 300;
    indAF = find(dataValidRaw(:,17) > CeruleanAF | dataValidRaw(:,18) > NeonAF );
    dataValid = dataValidRaw(indAF,:);
      
    ind0 = find(dataValid(:,22) == 0);
    ind1 = find(dataValid(:,22) == 1);
    %size normalize and group
    meanFSC = mean(dataValid([ind0;ind1],13));
    cellsTF0{i} = dataValid(ind0,19)./dataValid(ind0,13)*meanFSC;
    cellsTF1{i} = dataValid(ind1,19)./dataValid(ind1,13)*meanFSC;
    cellsTF{i} = dataValid([ind0;ind1],19)./dataValid([ind0;ind1],13)*meanFSC;
    cellsRec0{i} = dataValid(ind0,20)./dataValid(ind0,13)*meanFSC;
    cellsRec1{i} = dataValid(ind1,20)./dataValid(ind1,13)*meanFSC;
    cellsRec{i} = dataValid([ind0;ind1],20)./dataValid([ind0;ind1],13)*meanFSC;
    numCells0(i) = length(ind0);
    numCells1(i) = length(ind1);
    
    cellsCerulean0{i} = dataValid(ind0,17)./dataValid(ind0,13)*meanFSC;
    cellsCerulean1{i} = dataValid(ind1,17)./dataValid(ind1,13)*meanFSC;
    
    cellsCeruleanNeon{i} = [dataValid(:,17),dataValid(:,18)]./dataValid(:,13)*meanFSC;
   
    meanTF0(i) = mean(dataValid(ind0,19));
    meanTF1(i) = mean(dataValid(ind1,19));
    meanRec0(i) = mean(dataValid(ind0,20));
    meanRec1(i) = mean(dataValid(ind1,20));
    
    medTF0(i) = median(dataValid(ind0,19));
    medTF1(i) = median(dataValid(ind1,19));
    medRec0(i) = median(dataValid(ind0,20));
    medRec1(i) = median(dataValid(ind1,20));
    
    varTF0(i) = var(dataValid(ind0,19));
    varTF1(i) = var(dataValid(ind1,19));
    varRec0(i) = var(dataValid(ind0,20));
    varRec1(i) = var(dataValid(ind1,20));
end

stepsTF = 10;
minTF = -100;
maxTF = 2000;
edgesTF = minTF:stepsTF:maxTF;
minRec = -60;
maxRec = 100;
edgesRec = minRec:20:maxRec;
minAF = -20;
maxAF = 100;
edgesAF = minAF:20:maxAF;


