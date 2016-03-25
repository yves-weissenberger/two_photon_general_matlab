

r = 1;

%% I think this should calculate yes it does, creates the correlation matrix

%Ok the actual correlation thing works, it looks like it maaay be the
%indexing which is wrong

%magic happened, I think it works now; well, the output is kinda sensible.
%Maybe rounding errors are now causing the differences though I'm not sure.

numCells = size(GRABinfo.NPcSnglTrlRsp,1);

pop_response = GRABinfo.NPcSnglTrlRsp;


for kk=1:size(GRABinfo.NPcSnglTrlRsp,3)
    
    for k=1:size(GRABinfo.NPcSnglTrlRsp,3)
        
        xx = squeeze(GRABinfo.NPcSnglTrlRsp(:,:,kk))';
        y = squeeze(GRABinfo.NPcSnglTrlRsp(:,:,k))';
        
        
        corrmatrix = ((xx*y') - numCells*mean(xx,2)*mean(y,2)')./((numCells-1)*std(xx,0,2)*std(y,0,2)');
        
        mean_store(kk,k) = mean2(corrmatrix);
        
        
        
    end
    
end

% Ok, there is a problem that sometimes there are not all numbers
% represented uniquely. Some are not included in the cluster. When use
% cluster


%% plot the image based on only the reliable responses.
    maxColourMap = 0.6;


if r==1
    
    D = diag(mean_store);
    
    reliable_sounds = find(D>=0.2);
    
    reliable_store = mean_store([reliable_sounds],[reliable_sounds]);
    
    Zreliable = linkage(reliable_store,'single');
    [Hreliable,Treliable,OUTPERMreliable] = dendrogram(Zreliable,0);
    
    figure,
    imshow(reliable_store,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
    
    orderR = fliplr(OUTPERMreliable);
    clustRlblStore = reliable_store(orderR,orderR);
    
    
    figure,
    imshow(reliable_store(orderR,orderR),[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
    
    
    
elseif r==0
    
    %%
    
    
    
    %plot the images based on all responses
    Z = linkage(mean_store,'single');
    [H,T,OUTPERM] = dendrogram(Z,0);
    
    order = fliplr(OUTPERM);
    clustMeanStore = mean_store(order,order);
    
    figure,
    imshow(mean_store,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
    
    figure,
    imshow(clustMeanStore,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
    
    
    
    
    
    
end


%% Continuation

%can use something like clusterIdx = cluster(Z,'maxclust',2), to look at
%the identities of the clusters in the data and use this for further
%processing. This basically says where to cut down the dendrogram, below
%which horizontal line of the dendrogram to cut.


numClusters=2;

if r== 0
    
    ClustTree = Z;
    sequence = order;
    DataStore = clustMeanStore;
    
    
elseif r==1
    
    ClustTree = Zreliable;
    sequence = orderR;
    DataStore = clustRlblStore;
    
    
end


%Returns a list of cluster indices for each of the leaves
clusterIdx = cluster(ClustTree,'maxclust',numClusters);


%this variable is Huberts Gamma Statistic, a measure of how significant the
%clusters are





ClusterMtx1 = repmat(clusterIdx(sequence),1,length(sequence));
ClusterMtx2 = ClusterMtx1';

Clusters = ClusterMtx1.*ClusterMtx2;

    Clusters(find(Clusters==2)) = 0;
    Clusters(find(Clusters==3)) = 0;
    Clusters(find(Clusters==5)) = 0;
    Clusters(find(Clusters==6)) = 0;
    Clusters(find(Clusters==7)) = 0;
    Clusters(find(Clusters==8)) = 0;

    Clusters(find(Clusters==4)) = 2;
    Clusters(find(Clusters==9)) = 3;

figure,
imagesc(Clusters)



FClusters = Clusters;

    FClusters(find(FClusters==2)) = 1;
    FClusters(find(FClusters==3)) = 1;


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



meanReal = mean2(DataStore);
meanSim = mean2(FClusters);

stdReal = std2(DataStore);
stdSim = std2(FClusters);

tempH = 0;

for ii = 1:(length(DataStore)-1)
    
    for jj = ii+1:length(DataStore)

    tempH = tempH +  (DataStore(ii,jj)-meanReal).*(FClusters(ii,jj)-meanSim);
    
    
    end
    
end

HubertG = 1/((length(DataStore)*(length(DataStore)-1)*stdReal*stdSim))*tempH;


%sum(sum((DataStore - meanReal).*(FClusters - meanSim)))

%HubertG = ((1/stdReal*stdSim))*sum(sum((DataStore - meanReal).*(FClusters - meanSim)))



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