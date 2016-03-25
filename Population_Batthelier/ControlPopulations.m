%% Tones

r=1;
kkk= 1;
ShuffleType = 2;

%Cortex
%load('/Users/Yves/Desktop/20140304_Processed/Area04/(20140304_16_47_27)-_20140304_Area04_tones1_GRABinfo.mat')
%load('/Users/Yves/Desktop/20131017_Processed/Area02/(20131017_17_31_47)-_Cortex_20131017_Area02_tones1_GRABinfo.mat')

%MGB
load('/Users/Yves/Desktop/Boutons_preprocessed/20140116/Area02/(20140116_15_16_17)-_MGB20140116_Area02_tones1_GRABinfo.mat')
%load('/Users/Yves/Desktop/Boutons_preprocessed/20140121/Area01/(20140121_13_20_45)-_MGB20140121_Area01_tones1_GRABinfo.mat')

clear HStore

argg = randperm(200);

%pop_response = ModelResponses(:,:,argg(1:15));
%pop_response = ModelResponses;
%pop_response = GRABinfo.NPcSnglTrlRsp;
%pop_response = MixAreaPop;
%pop_response = ModelTest +0.01;
pop_response = GRABinfo.SnglTrlRsp(argg(1:30),:,:);


for kkk=1:10
    %% Shuffle
    
    numCells = size(pop_response,1);
    
     pop_response = pop_response;
     shufflePop_response = pop_response;
    
    if ShuffleType ==1
        
        
        %Shuffle stuff around by using
        Numlist = randperm(prod(size(pop_response)));
        
        %Then reshaping this.
        NumList = reshape(Numlist,[size(pop_response)]);
        
        
        %Shuffle only along one dimension
        
        
        ShuffleDim = 1;
        
        [dunno sortList] = sort(NumList,ShuffleDim);
        
        
        %Now I think you should be able to use sortList to shuffle the elements
        %just along one dimension. Should be able to do the others by. NO THIS
        %DOESN'T WORK. NEED indices along all dimensions. Use repmat and create
        %arrays of ascending numbers.
        
        
        
        %for each stimulus
        for ii = 1:100
            
            %for each response to each stimulus
            for jj=1:9
                
                %randomise the order of stimulus presentations
                
                
                shufflePop_response(:,jj,ii) = shufflePop_response(sortList(:,jj,ii),jj,ii);
                
            end
        end
        
        
        
        
    elseif ShuffleType==2
        
       
        
        for iii = 1:(size(pop_response,1))
            
            Numlist = randperm(prod(size(pop_response(1,:,:))));
            NumList = reshape(Numlist,[size(pop_response(1,:,:))]);
            
            popRespi = squeeze(pop_response(iii,:,:));
            shufflePop_response(iii,:,:) = popRespi(NumList);
            
            
        end
        
        
        
        
        
        
        
    end
    
    
    %% This performs the clustering analysis on the data where the identity of neurons are shuffled stimulus by stimulus
    
    
    for kk=1:size(pop_response,3)
        
        for k=1:size(pop_response,3)
            
            %Shufflexx and Shuffley are simply
            Shufflexx = squeeze(shufflePop_response(:,:,kk))';
            Shuffley = squeeze(shufflePop_response(:,:,k))';
            
            
            ShuffleCorrMatrix = ((Shufflexx*Shuffley') - numCells*mean(Shufflexx,2)*mean(Shuffley,2)')./((numCells-1)*std(Shufflexx,0,2)*std(Shuffley,0,2)');
            
            ShuffleMeanStore(kk,k) = mean2(ShuffleCorrMatrix);
            
            
            
        end
        
    end
    
    
    %This is just a small aside to check whether there are any reliable
    %responses after shuffling
%     D = diag(ShuffleMeanStore);
%     reliable_sounds{kkk} = find(D>=0.2);
    
    
    if r==1
        
        D = diag(ShuffleMeanStore);
        
        [ValY Ind] = sort(D,'descend');
        
        reliable_sounds = Ind(1:34);
        
        reliable_store = ShuffleMeanStore([reliable_sounds],[reliable_sounds]);
        
        Zreliable = linkage(reliable_store,'single');
        [Hreliable,Treliable,OUTPERMreliable] = dendrogram(Zreliable,0);
        
        % figure,
        % imshow(reliable_store,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
        
        orderR = fliplr(OUTPERMreliable);
        clustRlblStore = reliable_store(orderR,orderR);
        
        
        % figure,
        % imshow(reliable_store(orderR,orderR),[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
        
        
        
    elseif r==0
        
        %%
        
        maxColourMap = 0.6;
        
        
        %plot the images based on all responses
        Z = linkage(ShuffleMeanStore,'single');
        [H,T,OUTPERM] = dendrogram(Z,0);
        
        order = fliplr(OUTPERM);
        clustMeanStore = ShuffleMeanStore(order,order);
        
%         %figure,
%         imshow(ShuffleMeanStore,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
%         
%         %figure,
%         imshow(clustMeanStore,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
        
        
        
        
        
        
    end
    
    
    
    
    %% Continuation
    
    %can use something like clusterIdx = cluster(Z,'maxclust',2), to look at
    %the identities of the clusters in the data and use this for further
    %processing. This basically says where to cut down the dendrogram, below
    %which horizontal line of the dendrogram to cut.
    
    for jjjj = 1:24
        
        numClusters = 1 + jjjj;
        
        if r==0
            
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
        
        
        
        ClusterMtx1 = repmat(clusterIdx(sequence),1,length(sequence));
        ClusterMtx2 = ClusterMtx1';
        
        Clusters = zeros(size(ClusterMtx1));
        FClusters = Clusters;
        
        
        for ii=1:numClusters
            
            tempClustii = [];
            tempClustii = ClusterMtx1==ii & ClusterMtx2==ii;
            FClusters = FClusters + tempClustii;
            
            tempClustii = tempClustii.*ii;
            Clusters = Clusters + tempClustii;
            
            
        end
        
        
%         figure,
%         imagesc(Clusters)
%         
        %% This calculates Huberts Gamma Statistic, a measure of how well the clusters are separated
        if r == 1
            DataStore = reliable_store;
        elseif r == 0
            DataStore = ShuffleMeanStore;
        end
        
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
        
        
        HStore(kkk,jjjj) = HubertG;
        
        
        
    end
    
end

