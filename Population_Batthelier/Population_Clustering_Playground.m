clear HubertG

%% Load ROIs

%Interesting, looks very different
%load('/Users/Yves/Desktop/20140131/Area02/(20140131_13_06_12)-_20140131_Area02_tones2_GRABinfo.mat')
%load('/Users/Yves/Desktop/20140131/Area01/(20140131_12_32_09)-_20140131_Area01_tones2_GRABinfo.mat')


ROItype = 'cellBodies';


%dayDir = '/Users/Yves/Desktop/20140131';   %for boutons
dayDir = '/Users/Yves/Desktop/20140304_Processed';   %for cortical cell bodies
%dayDir = '/Users/Yves/Desktop/2p_MGB_imaging/20131105';    %for IC data
search_exp = '.*Area.*';
[Areafile_locs,~,~,~] = file_props(dayDir,search_exp);


numAreas = length(Areafile_locs);
%% 

clustIdx = 1;

for areaName=Areafile_locs
    Gnames=dir([fullfile(dayDir,areaName{:}) '/*.mat']);
    FileDir{clustIdx} =fullfile(dayDir,areaName{:},Gnames(3).name);
    clustIdx = clustIdx+1;
end

%for cortex 4,5,7 look good, 9 is interesting
%for IC, look at 2,4,5 for nice clusters
load(FileDir{4})
%%

%If only want to use the reliable population responses, set r to 0
reliableOnly = true;

%% This just makes sure you process the data in the correct way


if strcmp(ROItype,'cellBodies');
    pop_response = GRABinfo.NPcSnglTrlRsp;
    ToRemove = find(isnan(GRABinfo.NPcorrectedTraces(:,1)));
end


if strcmp(ROItype,'boutons');
    pop_response = GRABinfo.SnglTrlRspNew;
    ToRemove = find(isnan(GRABinfo.TracesNew(:,1)));
end

%if using a model, select one of these two
%pop_response = ModelResponses;
%pop_response = MixAreaPop;

%if loading IC data
if strcmp(ROItype,'IC')
    pop_response = GRABinfo.SnglTrlRsp;
    ToRemove = find(isnan(GRABinfo.NPTraces(:,1)));
end

%Remove the traces that are suspect
pop_response(ToRemove,:,:) = [];

%% Calculate the correlation matrix between responses.

%Initialise some global parameters
numCells = size(pop_response,1);
numFreqs = size(pop_response,3);




for freq1Idx=1:numFreqs
    
    for freq2Idx=1:numFreqs
        
        %Each row of respFreqn is a 
        respFreq1 = squeeze(pop_response(:,:,freq1Idx))';
        respFreq2 = squeeze(pop_response(:,:,freq2Idx))';
        
        %This gives you the single trial correlation matrix; ie element
        %(m,n) of the matrix contains the correlation between population
        %responses to the mth presentation of respFreq1 and the nth
        %presentation of the respFreq2
        trlCorrMat = ((respFreq1*respFreq2') - numCells*mean(respFreq1,2)*mean(respFreq2,2)')./((numCells-1)*std(respFreq1,0,2)*std(respFreq2,0,2)');
        
        %This gives you the correlation matrix, ie the correlation between
        %responses to each of the stimuli
        corrMat(freq1Idx,freq2Idx) = median(trlCorrMat(:));
        
        
        
    end
    
end


%This code is checked against a more basic implementation, there are
%differences, but on the order of 10^-16, so probably rounding errors..



%% Now want to cluster based on the correlation matrix

%This for loop basically runs the program over different thresholds
threshIdx=3;

for threshIdx=1:3
    
maxColourMap = 0.8;


if reliableOnly==true
    
    
    %% This will select then only the reliable responses
    
    %The method for doing this is to look for population responses to one
    %stimulus that have a correlation above a certain value
    
    %First select the correlations between responses to the same stimulus
    D = diag(corrMat);
    
    %find the threshold level
    thresh = 0.05 + 0.05*threshIdx;
    
    
    %find the indices of the sounds that are reliable
    reliable_sounds = find(D>=thresh);
    %find the indices of the sounds that are unreliable
    unreliable_sounds = find(D<=thresh);
    
    
    %This contains the number of reliable sounds at each threshold
    numReliableSounds(threshIdx) = length(reliable_sounds);

    %reduce the correlation matrix to only those pairs of sounds that are
    %reliable
    reliableCorrMat = corrMat([reliable_sounds],[reliable_sounds]);
    
    %This returns returns a matrix Zreliable that encodes a tree of hierarchical 
    %clusters of the rows of the real matrix reliableCorrMat. Z is a 
    %(m ? 1)-by-3 matrix, where m is the number of reliable frequencies.
    %The first two columns contain the nodes that are linked, the third
    %contains the linkage distance. 
    Zreliable = linkage(reliableCorrMat,'ward');
    
    %Create the dendrogram by magic. orderR puts the order of stimuli in
    %the correlation matrix to the right place
    [Hreliable,Treliable,orderR] = dendrogram(Zreliable,0);
    
    freqOrderR = ceil(orderR./4);
    freqLevelR = rem(orderR,4)./10;
    clustRlblStore = reliableCorrMat(orderR,orderR);
    
    
    figure('Units','Normalized'),
    subplot(1,2,1)
    pcolor(flipud(reliableCorrMat'));
    
    rFreqOrder = ceil(reliable_sounds./4);
    rLevel = rem(reliable_sounds,4)./10;
    ax = gca;
    ax.XTickLabel = rFreqOrder+rLevel;
    ax.XTick = 1:length(rFreqOrder);
    ax.YTick = 1:length(rFreqOrder);
    ax.YTickLabel = rFreqOrder+rLevel;
    axis square
    subplot(1,2,2)
    pcolor(flipud(reliableCorrMat(orderR,orderR)'))
    ax = gca;
    ax.XTickLabel = freqOrderR+freqLevelR;
    ax.XTick = 1:length(freqOrderR);
    ax.YTick = 1:length(freqOrderR);
    ax.YTickLabel = freqOrderR+freqLevelR;
    colorbar()
    axis square
    
    
elseif reliableOnly==0
    
    %%
    
    
    %plot the images based on all responses
    Z = linkage(corrMat,'single');
    [H,T,OUTPERM] = dendrogram(Z,0);
    
    %order = OUTPERM;
    order = fliplr(OUTPERM);
    clustMeanStore = corrMat(order,order);
    
    figure,
    %imshow(corrMat,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
    pcolor(corrMat')
    colorbar()
    figure,
    %imshow(clustMeanStore,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
    pcolor(clustMeanStore')
    colorbar()
    
    
    
    
    
end


%% Continuation

%can use something like clusterIdx = cluster(Z,'maxclust',2), to look at
%the identities of the clusters in the data and use this for further
%processing. This basically says where to cut down the dendrogram, below
%which horizontal line of the dendrogram to cut.


for numClusters = 1:25
    
    %numClusters = 1 + jjjj;
    
    if reliableOnly==0
        
        ClustTree = Z;
        sequence = order;
        DataStore = clustMeanStore;
        
        
    elseif reliableOnly==1
        
        ClustTree = Zreliable;
        sequence = orderR;
        DataStore = clustRlblStore;
        
        
    end
    
    
    %Returns a list of cluster indices for each of the leaves
    clusterIdx = cluster(ClustTree,'maxclust',numClusters);
    
    
    ClusterMtx1 = repmat(clusterIdx(sequence),1,length(sequence));
    ClusterMtx2 = ClusterMtx1';
    
    
    Clusters = zeros(size(ClusterMtx1));
    FClusters = Clusters;
    
    
    for clustIdx=1:numClusters
        
        tempClustii = [];
        tempClustii = ClusterMtx1==clustIdx & ClusterMtx2==clustIdx;
        FClusters = FClusters + tempClustii;
        
        tempClustii = tempClustii.*clustIdx;
        Clusters = Clusters + tempClustii;
        
        
    end
    
    if numClusters <=3
     figure,
     pcolor(flipud(Clusters')); colorbar();
     ax = gca;
     ax.XTickLabel = freqOrderR+freqLevelR;
     ax.XTick = 1:length(freqOrderR);
     ax.YTick = 1:length(freqOrderR);
     ax.YTickLabel = freqOrderR+freqLevelR;
    end
    
    %% This calculates Huberts Gamma Statistic, a measure of how significant the clusters are
    
    
    % n = length(FClusters);
    %
    % FClusters=FClusters';
    % FClusters(1:n+1:n*n)=[];
    % FClusters=reshape(FClusters,n-1,n)';
    %
    % DataStore=DataStore';
    % DataStore(1:n+1:n*n)=[];
    % DataStore=reshape(DataStore,n-1,n)';
    
    
    %Calculates the mean_correlation
    meanReal = mean2(DataStore);
    
    %
    meanSim = mean2(FClusters);
    
    
    stdReal = std2(DataStore);
    
    
    stdSim = std2(FClusters);
    
    
    tempH = 0;
    
    for clustIdx = 1:(length(DataStore)-1)
        
        for jj = clustIdx+1:length(DataStore)
            
            tempH = tempH +  (DataStore(clustIdx,jj)-meanReal).*(FClusters(clustIdx,jj)-meanSim);
            
            
        end
        
    end
    
    HubertG(threshIdx,numClusters) = 1/((length(DataStore)*(length(DataStore)-1)*stdReal*stdSim))*tempH;
    
    
    
end

end



%% 

for clustIdx=1:3
    plot(HubertG(clustIdx,:))
    hold on
end



%% Other things



%% This is a check the part with the heading: Calculate the correlation matrix between responses.

freq1Idx = 4;
freq2Idx = 34;

respFreq1 = squeeze(pop_response(:,:,freq1Idx))';
respFreq2 = squeeze(pop_response(:,:,freq2Idx))';


trlCorrMat = ((respFreq1*respFreq2') - numCells*mean(respFreq1,2)*mean(respFreq2,2)')./((numCells-1)*std(respFreq1,0,2)*std(respFreq2,0,2)');


for clustIdx = 1:9
    
    for jj=1:9
        temp = corrcoef(respFreq1(clustIdx,:),respFreq2(jj,:));
        trlCheck(clustIdx,jj) = temp(2);
    end
end

%% Some other stuff

% %% Calculate the Firing Rate in response to different tones
%
% AcrossTrlRsp = squeeze(mean(GRABinfo.NPcSnglTrlRsp,2));
%
% AcrossCellRsp = squeeze(mean(AcrossTrlRsp,1));
%
% figure,
% imagesc(AcrossCellRsp(order));
%
%
%
%
% %% Look at the stimulus Lists
%
% numFreq = 25;
% numIntens = 4;
%
% stimS = repmat(1:numFreq,[numIntens 1]);
% stimS = stimS(:);
%
%
%
%
% subplot(2,1,2), imshow(stimS(order)',[1 25],'InitialMagnification',100,'Colormap',jet(255))