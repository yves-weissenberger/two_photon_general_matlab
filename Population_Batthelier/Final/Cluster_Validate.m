function [Clusters,HubertG,HubertG2,FClusterStore,freqOrder,freqLevel] = Cluster_Validate(ClustTree,clustStore,order,reliableSounds,toPlot)


for numClusters = 1:25
    
    
    %Returns a list of cluster indices for each of the leaves
    clusterIdx = cluster(ClustTree,'maxclust',numClusters);
    
    
    ClusterMtx1 = repmat(clusterIdx(order),1,length(order));
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
    
    %Store for each number of clusters the matrix
    FClusterStore(numClusters,:,:) = FClusters;
    freqOrder = ceil(reliableSounds./4);
    freqLevel = rem(reliableSounds,4)./10;
    
    if (numClusters <=4 && numClusters>1 && toPlot==true)
        
        figure('Units','normalized')
        pcolor(flipud(Clusters')); colorbar();
        ax = gca;
        ax.XTickLabel = flipud(freqOrder(order)+freqLevel(order));
        ax.XTick = 1:length(freqOrder);
        ax.YTick = 1:length(freqOrder);
        ax.YTickLabel = freqOrder(order)+freqLevel(order);
        ax.FontSize = 7.5;
    end
    
%% This calculates Huberts Gamma Statistic, a measure of how significant the clusters are
    
%From what I understand, the statistic works as follows    
 
    
    %Calculates the mean_correlation across the correlation matrix
    
    testLocs = find(FClusters~=0);
    IndicatorMtx = FClusters;
    IndicatorMtx(testLocs) = 1;
    %IndicatorMtx = ~IndicatorMtx*1;                                        %ignore this line, do NOT unhash it

    
    meanReal = mean2(clustStore);
    
    %
    meanSim = mean2(IndicatorMtx~=0);
    
    
    stdReal = std2(clustStore);
    
    
    stdSim = std2(IndicatorMtx~=0);
    
    tempH = 0;
    
    for ii = 1:(length(clustStore)-1)
        
        for jj = ii+1:length(clustStore)
            
            %take each little square, want to compare
            tempH = tempH +  (clustStore(ii,jj)-meanReal).*(IndicatorMtx(ii,jj)-meanSim);
            %In this equation, if you are looking outside of a cluster, the
            %second term will be positive. If you are looking inside of a
            %cluster it will be negative.
            
        end
        
    end
    
    HubertG(numClusters) = 1/((length(clustStore)*(length(clustStore)-1)*stdReal*stdSim))*tempH;
    temp2 = corrcoef(clustStore(:),IndicatorMtx(:));
    HubertG2(numClusters) = temp2(2);
    
    
end