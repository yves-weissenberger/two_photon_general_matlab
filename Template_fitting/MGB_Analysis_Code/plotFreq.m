function [FRA_all,NeuronXYZ_store,BF_store,CF_store,CFband_Store,b21,b22,stats21,stats22,b3,bint3,r3,rint3,stats3,b,stats] = plotFreq(targetDir,Date,Type,oneD,MaxResp,cellBodies)


BF_store = [];s
NeuronXYZ_store = [];
NeuronXYZ_storeCF = [];
responsive_store = [];
kk = 0;
CF_store = [];
CFband_Store = [];
BiggestRespStore = [];

stimlist = logspace(log10(1.250),log10(80.000),25);

numAreas = length(dir([fullfile(targetDir,Date)])) - 3;

kkkk = 0;

for ii = [1:numAreas-3]
    %% Admin :/

    if ii>=10
        
        
        kk = [];
        
    end
    
    Area = sprintf('Area%d%d', kk,ii);
    
    %dir returns all files in a given directory. f = fullfile('myfolder','mysubfolder','myfile.m')
    %returns a string combining all the inputs
    Area_files = dir([fullfile(targetDir,Date,Area) '/*.mat']);
    
    %Load image to be processed
    GRABname=Area_files(2).name;
    load(fullfile(targetDir,Date,Area,GRABname));
    
    if cellBodies == 1
       pop_response = GRABinfo.NPcSnglTrlRsp;
    else
       pop_response = GRABinfo.SnglTrlRspNew;
    end
    
    numCells = size(pop_response,1);
    
    responsive = zeros(1,size(pop_response,1));
    
    %This resticts analysis to either areas low down or high up.
    if GRABinfo.xyzPosition(3) <= -200
    
    kkkk = kkkk +1;
    
    %% First find out what of your cells are actually responding
    
    %Remove the erroneously labelled ROIs
    if cellBodies == 1
    ToRemove = find(isnan(GRABinfo.NPcorrectedTraces(:,1)));
    pop_response(ToRemove,:,:) = [];
    else
    ToRemove = find(isnan(GRABinfo.TracesNew(:,1)));
    pop_response(ToRemove,:,:) = [];
    end
    
   
    
    
    responsive = zeros(1,size(pop_response,1));
    
    
    %This is the working script
    BF = [];
    BFint = [];
    CF = [];
    NeuronXYZ = [];
    CFband = [];
    
    for j=1:size(pop_response,1);
        p=anova1(squeeze(pop_response(j,:,:)),[],'off');
        
        if p<0.001;
            
            %FRA is a 4x25 matrix, of sound intensities against the
            %frequencies.
            FRA=flipud(reshape(mean(pop_response(j,:,:),2),4,25));
            FRAreshp(j,:)=FRA(:);
            FRA_all(:,:,j,kkkk)=FRA;
            
            %Returns the row (BFint) and column (BF) of the FRA in which the neuron has its best
            %response
            
            [BFint(j), BF(j)] = ind2sub(size(FRA), find(FRA == max(max(FRA))));
            
            responsive(j) = 1;
            
        end
    end
    
    %not entirely sure of this yet, check that the ANOVA thing works. Look at
    %plots, by removing ,[],'off'
    
    
    BF(find(BF==0)) = [];
    BFint(find(BFint==0)) = [];
    
    
    %% Find CF
    
    %This finds those cells that are responsive
    responsiveCells = find(responsive);
    
    
    jj = 1;
    
    
    for jj = 1:sum(responsive)
        
        Resp = reshape(squeeze(pop_response(responsiveCells(jj),:,:)),size(pop_response,2),4,25);
        
        %now have a 9x25x4 array and everything works out ok
        responses = squeeze(permute(Resp,[1,3,2]));
        
        
        % %% Script to show that this this permute and reshape still makes sense, plot FRA
        %
        %          FRA2 = flipud(squeeze(mean(responses,1))');
        %         imagesc(FRA2)
        
        
        pval=anova1(responses(:,:,1),[],'off');
        
        if pval <=0.05
            CF(jj) = find(FRA_all(1,:,responsiveCells(jj),kkkk)==max(FRA_all(1,:,responsiveCells(jj),kkkk)));
            CFband(jj) = 1;
        else
            pval=anova1(responses(:,:,2),[],'off');
            
            if pval <=0.05
                CF(jj) = find(FRA_all(2,:,responsiveCells(jj),kkkk)==max(FRA_all(2,:,responsiveCells(jj),kkkk)));
                CFband(jj) = 2;
                
            else
                
                pval=anova1(responses(:,:,3),[],'off');
                
                if pval<=0.05
                    CF(jj) = find(FRA_all(3,:,responsiveCells(jj),kkkk)==max(FRA_all(3,:,responsiveCells(jj),kkkk)));
                    CFband(jj) = 3;
                else
                    CF(jj) = find(FRA_all(4,:,responsiveCells(jj),kkkk)==max(FRA_all(4,:,responsiveCells(jj),kkkk)));
                    CFband(jj) = 4;
                end
                
                
            end
            
        end
        
        
        
    end
    
    %%
    BF_store = cat(2,BF_store,BF);
    BF_elements =  unique(BF_store);
    
    CF_store = cat(2,CF_store,CF);
    CF_elements = unique(CF_store);
    CFband_Store = cat(2,CFband_Store,CFband);
    
    responsive_store = cat(2,responsive_store,responsive);
    
%      for iii = 1:length(unique(CF_store))
%          CF_store2(find(CF_store == CF_elements(iii))) = stimlist(iii);
%      end
%     
%     for iii = 1:length(unique(BF_store))
%         BF_store(find(BF_store == BF_elements(iii))) = stimlist(iii);
%     end
    
    %% Find the positions of the neurons relative to the centre of the window
    
    if cellBodies == 1
    ROICentroids = [(GRABinfo.ROIprops(:).Centroid)];
    else
    ROICentroids = [(GRABinfo.ROIpropsNew(:).Centroid)];
    end
    
  
        
        %SHIT think there is a mistake in this script. I think it selects
        %the wrong responsive ROIs.
    NeuronXYZ = reshape(ROICentroids,2,numCells)';
    NeuronXYZ = NeuronXYZ + [ones(size(NeuronXYZ,1),1).*GRABinfo.xyzPosition(1),ones(size(NeuronXYZ,1),1).*GRABinfo.xyzPosition(2)];
    
    NeuronXYZ(:,3) = ones(size(NeuronXYZ,1),1).*GRABinfo.xyzPosition(3);
    
    
        NeuronXYZ(ToRemove,:) = [];
        
     
          %This gets rid of the neurons that are not significant with ANOVAs
          NeuronXYZ(find(responsive==0),:) = [];
          NeuronXYZ_store = cat(1,NeuronXYZ_store,NeuronXYZ);

          
          %Process CF stuff further
          NeuronXYZ_CF = NeuronXYZ;
          NeuronXYZ_CF(find(CFband==4),:) = [];
          NeuronXYZ_storeCF = cat(1,NeuronXYZ_storeCF, NeuronXYZ_CF);
    
    
    meanResp = squeeze(mean(pop_response,2));
    meanResp(find(responsive==0),:) = [];
    BiggestResp = max(meanResp,[],2);
    BiggestRespStore = cat(1,BiggestRespStore,BiggestResp);


    
    end
end

CF_store(find(CFband_Store==4))=[];
BF_storeCFmatch = BF_store;
BF_storeCFmatch(find(CFband_Store==4))=[];

if Type == 'cf'
    
    Store = CF_store;
    
elseif Type =='bf'
    
    Store = BF_store;
    
elseif Type == 'cb'
    
    figure,  scatter3(NeuronXYZ_storeCF(:,1),NeuronXYZ_storeCF(:,2),NeuronXYZ_storeCF(:,3),ones(length(CF_store),1)*70,CF_store,'fill')
    xlabel('Rostral to Caudal')
    ylabel('Medio-Lateral')
    zlabel('Depth')
    title('Characteristic Frequency')
    
    figure,  scatter3(NeuronXYZ_store(:,1),NeuronXYZ_store(:,2),NeuronXYZ_store(:,3),ones(length(BF_store),1)*70,BF_store,'fill')
    xlabel('Rostral to Caudal')
    ylabel('Medio-Lateral')
    zlabel('Depth')
    title('Best Frequency')

end


if Type ~= 'cb'
    figure,
    scatter3(NeuronXYZ_store(:,1),NeuronXYZ_store(:,2),NeuronXYZ_store(:,3),ones(length(BF_store),1)*70,Store,'fill')
    
end



X = [ones(size(NeuronXYZ_store,1),1),NeuronXYZ_store(:,1),NeuronXYZ_store(:,2),NeuronXYZ_store(:,3)];
[b,bint,r,rint,stats] = regress(BF_store',X);

if oneD == 1
    
    figure,
    X21 = [ones(size(NeuronXYZ_store,1),1),NeuronXYZ_store(:,1)];
    [b21,bint21,r21,rint21,stats21]  = regress(BF_store',X21);
    
    plot(NeuronXYZ_store(:,1),BF_store,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor', [0 0 0])
    hold on
    yFit =  b21(1) + b21(2)*NeuronXYZ_store(:,1);
    plot(NeuronXYZ_store(:,1),yFit,'r')
    title('Best Frequency')
        xlabel('Rostral to Caudal')

    
    hold off
    
    
    figure,
    X22 = [ones(size(NeuronXYZ_storeCF,1),1),NeuronXYZ_storeCF(:,1)];
    [b22,bint22,r22,rint22,stats22]  = regress(CF_store',X22);
    
    plot(NeuronXYZ_storeCF(:,1),CF_store,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor', [0 0 0])
    hold on
    yFit =  b22(1) + b22(2)*NeuronXYZ_storeCF(:,1);
    plot(NeuronXYZ_storeCF(:,1),yFit,'r')
    title('Characteristic Frequency')
        xlabel('Rostral to Caudal')

    
    hold off
    
end
    
    binW = 40;

    %This bins position into discrete locations
    RoundPosition = round(NeuronXYZ_store(:,1)/binW)*binW;
    
    BinLocations = unique(RoundPosition);
    numPositions = length(BinLocations);
    
    for iiii = 1:numPositions
        
        CellPos = find(RoundPosition==BinLocations(iiii));
        MeanBF(iiii) = round(mean(BF_store(find(RoundPosition==BinLocations(iiii)))));
        
    end
        
    figure,     
    plot(MeanBF,'o')
    


    
%% This part of the script works out the pairwise distances for Best Frequency


[Xpos Xorder] = sort(NeuronXYZ_store(:,1));

%PW stands for pairwise
BF_storePW = BF_store(Xorder);

PW_BFdifferenceStore = [];
PW_distanceStore = [];
NeuronXYZ_PW = NeuronXYZ_store(Xorder,1);

for iii = 1:length(Xorder)
    
    
    PW_BFdifference = abs(BF_storePW(iii) - BF_storePW(iii+1:length(Xorder)));
    PW_distance = NeuronXYZ_PW(iii,1) - NeuronXYZ_PW(iii+1:length(Xorder),1);

    PW_BFdifferenceStore = cat(2,PW_BFdifferenceStore,PW_BFdifference);
    PW_distanceStore = cat(1,PW_distanceStore,PW_distance);

    
end







X3 = [ones(size(PW_distanceStore,1),1),PW_distanceStore];
     [b3,bint3,r3,rint3,stats3]  = regress(PW_BFdifferenceStore',X3);

    
    figure,
    plot(PW_distanceStore,PW_BFdifferenceStore,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor', [0 0 0]);

    hold on
    yFit =  b3(1) + b3(2)*PW_distanceStore';
    plot(PW_distanceStore',yFit,'r')
    title('Pairwise Best Frequency difference')
        xlabel('Rostral to Caudal')

    
    hold off
    
    figure,
    [n,c]=hist3([PW_BFdifferenceStore',PW_distanceStore],[length(unique(PW_BFdifferenceStore)),300]);
    imagesc(n)
    [n,c]=hist3([PW_BFdifferenceStore',PW_distanceStore],[length(unique(PW_BFdifferenceStore)),300]);


%% Same as above, but for Characteristic Frequency
% [Xpos Xorder] = sort(NeuronXYZ_storeCF(:,1));
% 
% %PW stands for pairwise
% CF_storePW = CF_store(Xorder);
% 
% PW_CFdifferenceStore = [];
% PW_distanceStore = [];
% NeuronXYZ_PW = NeuronXYZ_storeCF(Xorder,1);
% 
% for iii = 1:length(Xorder)
%     
%     
%     PW_BFdifference = abs(CF_storePW(iii) - CF_storePW(iii+1:length(Xorder)));
%     PW_distance = NeuronXYZ_PW(iii,1) - NeuronXYZ_PW(iii+1:length(Xorder),1);
% 
%     PW_CFdifferenceStore = cat(2,PW_CFdifferenceStore,PW_BFdifference);
%     PW_distanceStore = cat(1,PW_distanceStore,PW_distance);
% 
%     
% end
% 
% X3 = [ones(size(PW_distanceStore,1),1),PW_distanceStore];
%      [b3,bint3,r3,rint3,stats3]  = regress(PW_CFdifferenceStore',X3);
% 
%     
%     figure,
%     plot(PW_distanceStore,PW_CFdifferenceStore,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor', [0 0 0]);
% 
%     hold on
%     yFit =  b3(1) + b3(2)*PW_distanceStore';
%     plot(PW_distanceStore',yFit,'LineWidth',1,'Color','r')
%     title('Pairwise Characteristic Frequency difference')
%         xlabel('Rostral to Caudal')
% 
%     
%     hold off
%     
%     figure,
%     [n,c]=hist3([PW_CFdifferenceStore',PW_distanceStore],[length(unique(PW_CFdifferenceStore)),300]);
%     imagesc(n)
%     [n,c]=hist3([PW_CFdifferenceStore',PW_distanceStore],[length(unique(PW_CFdifferenceStore)),300]);

end
%%X is caudorostral. Its -2000 is rostral; 2000 is caudal.

%% Peters Magic
% [n,c]=hist3([BFstore',CFstore'],[25,25]);
% imagesc(n)
% [n,c]=hist3([BFstore',CFstore'],[25,25]);

%look at bandwidths of an area, FRA shape. Compare distance of neurons.  