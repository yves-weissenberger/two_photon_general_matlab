clear all


%% Load the data


%dayDir = '/Users/Yves/Desktop/20140131';  ROItype = 'boutons';                %for boutons
dayDir = '/Users/Yves/Desktop/20140304_Processed'; ROItype = 'cellBodies';    %for cortical cell bodies
%dayDir = '/Users/Yves/Desktop/2p_MGB_imaging/20131105';  ROItype = 'IC';       %for IC data
search_exp = '.*Area.*';
[Areafile_locs,~,~,~] = file_props(dayDir,search_exp);


numAreas = length(Areafile_locs);

areaIdx = 1;

for areaName=Areafile_locs
    Gnames=dir([fullfile(dayDir,areaName{:}) '/*.mat']);
    FileDir{areaIdx} =fullfile(dayDir,areaName{:},Gnames(3).name);
    areaIdx = areaIdx+1;
end

%for cortex 4,5,7 look good, 9 is interesting
%for IC, look at 2,4,5 for nice clusters
load(FileDir{4})


%% Clean up the data


if strcmp(ROItype,'cellBodies');
    popResponse = GRABinfo.NPcSnglTrlRsp;
    ToRemove = find(isnan(GRABinfo.NPcorrectedTraces(:,1)));
end


if strcmp(ROItype,'boutons');
    popResponse = GRABinfo.SnglTrlRspNew;
    ToRemove = find(isnan(GRABinfo.TracesNew(:,1)));
end


%if loading IC data
if strcmp(ROItype,'IC')
    popResponse = GRABinfo.SnglTrlRsp;
    ToRemove = find(isnan(GRABinfo.NPTraces(:,1)));
end

%Remove the traces that are suspect and count number of cells remaining
popResponse(ToRemove,:,:) = [];
numCells = size(popResponse,1);


%This is to make the number of ROIs similar to what is seen in the real
%data
if strcmp(ROItype,'boutons');
    popResponse = popResponse(randperm(numCells,30),:,:);
end

%%  Cluster the Population Responses

%set  parameters
reliableOnly = true;   %decision on whether to include only reliable responses
thresh = 0.15;       %The mean pariwise correlation threshold for considering a response reliable


%This is the function that does everything
[corrMat,reliableSounds,reliableCorrMat,ClustTree,order,clustStore,Clusters,HubertG,HubertG2,FClusterStore,freqOrder,freqLevel]...
    = Cluster_popResponse(popResponse,thresh,reliableOnly,'toPlot=1');



%% Generate Control data

%for varargin, one argument will be toPlot, which makes you show an overlay
%of the relation between responses of one neuron and that of other neurons
%. The second argument is controlType=x where x is either PCA, Uniform or FRA.
%This generates different types of control data
[ModelResponses] = Generate_Control_Data(popResponse,'controlType=FRA','toPlot=1');


%% Cluster the Model Population Responses

[CRTLcorrMat,CTRLreliableSounds,CTRLreliableCorrMat,CTRLClustTree,CTRLorder,CTRLclustStore,CTRLClusters,CTRLHubertG,CTRLHubertG2,CTRLFClusterStore,CRTLfreqOrder,CRTLfreqLevel] ...
    = Cluster_popResponse(ModelResponses,thresh,reliableOnly,'toPlot=1');

%% Plot the gap statistic

plot(HubertG - CTRLHubertG)