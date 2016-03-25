

%load('/Users/Yves/Desktop/20140304_Processed/Area04/(20140304_16_47_27)-_20140304_Area04_tones1_GRABinfo.mat')
%load('/Users/Yves/Desktop/20131017_Processed/Area02/(20131017_17_31_47)-_Cortex_20131017_Area02_tones1_GRABinfo.mat')


pop_response = GRABinfo.NPcSnglTrlRsp;

%I think this works, but am not sure. I don't think it does
mean_resp = squeeze(mean(pop_response,2));


meanMtx = permute(repmat(mean_resp,1,1,9),[1 3 2]);

meanList = reshape(meanMtx,prod(size(meanMtx)),1);

% %%
% 
% 
% %Think about what assumptions I want to and can make about responses. Q.
% %Does the variability of the response scale with the mean response? Simple,
% %look for a correlation coefficient between the two, for a single neuron.
% %Then do this for all neurons do see if you have the same relation
% 
% clear FinalStore
% clear normRespStore
% 
% FinalStore = zeros(size(pop_response,3),size(pop_response,1)*size(pop_response,2));
% 
% 
% %for all stimuli
% for ii = 1:size(pop_response,3)
%     
%     normRespStore = [];
%     
%     %for all neurons
%     for jj = 1:size(pop_response,1)
%         
%         
%         meanResp = mean(pop_response(jj,:,ii));
%         
%         maxResp = max(pop_response(jj,:,ii));
%         
%         normResp = (pop_response(jj,:,ii) - meanResp)./maxResp;
%         
%         normRespStore = cat(2,normRespStore,normResp);
%         
%     end
%     
%     FinalStore(:,ii) = normRespStore;
% end
% 
% FinalSlinear = reshape(FinalStore,1,prod(size(FinalStore)));
% 
% 
% %This will plot a histogram that plots the mean subtracted then normalised
% %responses of neurons.
% %If have large mean response would expect to have variance, would expect
% %variance proportional to the firing rate.
% figure, hist(FinalSlinear,500)
% 
% 
% 
% % This is kind of relevant
% %figure, hist(normrnd(0,1,100,1,1))

%% Look for whether response variability scales with the response amplitude for single neurons

pop_response = GRABinfo.NPcSnglTrlRsp;
%pop_response=ModelResponses;


%For test purposes, note that neuron 12 is very nicely responsive
neuronNum = 13;

%Get the order of responses, greatest to smallest for one neuron

SnglNeuronMean = squeeze(mean(pop_response,2));

%The std of the response of neuron m to stimulus n is given by
%SnglNeuronStd(m,n)
SnglNeuronStd = squeeze(std(pop_response,[],2));


%with regress, each row is an observation
[b,bint,r,rint,stats] = regress(SnglNeuronStd(neuronNum,:)',[ones(size(SnglNeuronMean,2),1),SnglNeuronMean(neuronNum,:)']);

plot(SnglNeuronMean(neuronNum,:),SnglNeuronStd(neuronNum,:),'o')
hold on

yFit =  b(1) + b(2)*SnglNeuronMean(neuronNum,:);
plot(SnglNeuronMean(neuronNum,:),yFit,'r')

hold off


%% Do the same for the whole population


NeuronsStdStore = reshape(SnglNeuronStd,prod(size(SnglNeuronStd)),1);
NeuronsMeanStore = reshape(SnglNeuronMean,prod(size(SnglNeuronMean)),1);


[b2,b2int,r2,r2int,stats2] = regress(NeuronsStdStore,[ones(size(NeuronsMeanStore,1),1),NeuronsMeanStore]);

%from plot, looks like the residuals are pretty linear
figure, 

hold on
plot(NeuronsMeanStore,yFit2,'r')
plot(NeuronsMeanStore,NeuronsStdStore,'o')

yFit2 =  b2(1) + b2(2)*NeuronsMeanStore;
plot(NeuronsMeanStore,yFit2,'r')
xlabel('Mean Response of Each Neuron to each stimulus')
ylabel('Std of the responses of each neuron')
legend('Line of Best Fit','Single Neuron, single Stimulus Mean & STD')
hold off



%Ok, so since there is a linear relationship between the Std of neural
%responses and their mean





%% Also want to find out the distribution of responses around the mean
%
%
%
%% First for one neuron

%load('/Users/Yves/Desktop/20140304_Processed/Area04/(20140304_16_47_27)-_20140304_Area04_tones1_GRABinfo.mat')

pop_response = GRABinfo.NPcSnglTrlRsp;
%pop_response = ModelResponses;



numCells = size(GRABinfo.NPcSnglTrlRsp,1);
numStim = size(GRABinfo.NPcSnglTrlRsp,3);


clear RespStLinear

ii = 1;


for jj = 1:numStim
    
    
    %response_store(1,:,jj) = (pop_response(ii,:,jj))/(mean(pop_response(ii,:,jj))*0.63);
    
    response_store(1,:,jj) = (pop_response(ii,:,jj)*100- mean(pop_response(ii,:,jj))*100)/(mean(pop_response(ii,:,jj))*0.63*100);
    
    
end

%This is now a 9x100 variable of repetitions in rows and stim in
%columns
Response_Store = squeeze(response_store);


%Ok, so this is the order of response size in the population
[Value OrderRespStore] = sort(mean(pop_response(ii,:,:),2),3,'ascend');


%This should order the responses, such that the mean response in
%column 1 is the smallest etc
Response_Store = Response_Store(:,OrderRespStore);



RespStLinear = reshape(Response_Store(:,66:100),1,prod(size(Response_Store(:,66:100))));


%This plots the range of responses of one neuron to all different stimuli,
%with the variation of the responses normalised across different stimuli
figure, hist(RespStLinear,10)








%% Now for all neurons

%Cell Bodies
%load('/Users/Yves/Desktop/20140304_Processed/Area04/(20140304_16_47_27)-_20140304_Area04_tones1_GRABinfo.mat')
pop_response = GRABinfo.NPcSnglTrlRsp;



%% Boutons
%load('/Users/Yves/Desktop/Boutons_preprocessed/20140121/Area01/(20140121_13_20_45)-_MGB20140121_Area01_tones1_GRABinfo.mat')
pop_response = GRABinfo.SnglTrlRsp;

%pop_response = ModelResponses;


numCells = size(pop_response,1);
numStim = size(pop_response,3);


response_store = zeros(1,9,100);

%First, get Std of responses independent of the mean response

clear RespStLinear2

for ii = 1:numCells
    
    
    
    for jj = 1:numStim
        
        
        %response_store(1,:,jj) = (pop_response(ii,:,jj))/(mean(pop_response(ii,:,jj))*0.63);
        
        response_store(1,:,jj) = (pop_response(ii,:,jj)*100- mean(pop_response(ii,:,jj))*100)/(mean(pop_response(ii,:,jj))*0.63*100);
        
        %This is now a 9x100 variable of repetitions in rows and stim in
        %columns
        Response_Store = squeeze(response_store);
        
        
        %Ok, so this is the order of response size in the population
        [Value OrderRespStore] = sort(mean(pop_response(ii,:,:),2),3,'ascend');
        
        
        %This should order the responses, such that the mean response in
        %column 1 is the smallest etc
        Response_Store = Response_Store(:,OrderRespStore);
        
        
        
    end
    
    
    
    
    
    RespStLinear2(ii,:) = reshape(Response_Store(:,70:100),1,prod(size(Response_Store(:,70:100))));
    
end


VarHistData = reshape(RespStLinear2,1,prod(size(RespStLinear2)));


clear pd xvalues pdf



[nelements,centers] = hist(VarHistData,60);

bar(centers,nelements/sum(nelements));


% hold on
% 
% pd = fitdist(VarHistData','Normal');
% xvalues = [-15:1:15];
% pdf = pdf(pd,xvalues);
% plot(xvalues',pdf,'LineWidth',2)
% 
% hold off
% 
% 


%% Fit the Variability to key distribution Parameters

moments = {mean(VarHistData),std(VarHistData),skewness(VarHistData),kurtosis(VarHistData)};
%rng('default');  % For reproducibility
[r,type] = pearsrnd(moments{:},1000,1);
figure, hist(r,30)



%% Fit Smart Model

ModelResponses = zeros(numCells,9,100);

%NOT SURE WHY, BUT THESE PERMUTATIONS DON'T SEEM TO WORK
sFRAall = permute(reshape(smFRA_store,100,size(smFRA_store,3)),[2 1]);
%FRAall = permute(reshape(FRAall,100,size(FRAall,3)),[2 1]);


    %pop_response = GRABinfo.NPcSnglTrlRsp;

    FRAStore= squeeze(mean(pop_response,2));
    
    
    

for ii = 1:size(sFRAall,1)
    
    %Selects the FRA of a given cell
    CellFRA = smFRA_store(:,:,ii);
    
    %reshapes this FRA to be linear
    CellLinFRA = FRAStore(ii,:);
    %CellLinFRA = reshape(permute(CellFRA,[2 1]),1,100);
    
    


    
    %Run through all the stimuli
    for jj=1:size(CellLinFRA,2)
        
        
        %This is supposed to multiply the mean responses of a given cell to a given tone by a series of random numbers 
        ModelResponses(ii,:,jj) = pearsrnd(CellLinFRA(jj),CellLinFRA(jj)+0.3,moments{3},moments{4},1,9,1);

        
    end
    
end


ModelResponses(find(isnan(ModelResponses))) = 0;


%% Simple Dumb Model


ModelResponses = zeros(size(FRAall,1),9,size(FRAall,2));

%This creates the sFRAall vector, where
%sFRAall = permute(reshape(smFRA_store,100,size(smFRA_store,3)),[2 1]);

%Run through all the cells
for ii = 1:size(sFRAall,1)
    
    %Selects the FRA of a given cell
    CellFRA = smFRA_store(:,:,ii);
    
    %reshapes this FRA to be linear
    CellLinFRA = reshape(permute(CellFRA,[1 2]),1,100);


    
    %Run through all the stimuli
    for jj=1:size(sFRAall,2)
        
        
        %This is supposed to multiply the mean responses of a given cell to a given tone by a series of random numbers 
        ModelResponses(ii,:,jj) = CellLinFRA(jj).*normrnd(0.5,2,1,9,1);
        
        %NewFRAStore = 

        
    end
    
end



%% This creates a population of neurons that just have a line down their FRA


%Want to create a simple model population

clear simFRA_all

simFRA_all = zeros(4,25,24);

for ii = 1:24
    
    FRATemplate = zeros(4,25);
    
    %This gives overlapping R-FIELDS
    %FRATemplate(:,ii:ii+1) = 1;
    
    
     FRATemplate(:,ii) = 1;

    simFRA_all(:,:,ii) = FRATemplate;
    
end



%Run through all the cells
for ii = 1:size(simFRA_all,3)
    
    %Selects the FRA of a given cell
    CellFRA = simFRA_all(:,:,ii);
    
    %reshapes this FRA to be linear
    CellLinFRA = reshape(permute(CellFRA,[1 2]),1,100);


    
    %Run through all the stimuli
    for jj=1:size(CellLinFRA,2)
        
        
         %This is supposed to multiply the mean responses of a given cell to a given tone by a series of random numbers 
        %ModelResponses(ii,:,jj) = pearsrnd(CellLinFRA(jj),CellLinFRA(jj),moments{3},moments{4},1,9,1);
        
        %Now generate a non-stochastic response
        ModelResponses(ii,:,jj) = CellLinFRA(jj);
        

        

        
    end
    
end


%Generates noise necessary so do not have 0 responses
%ModelResponses = ModelResponses + rand(size(ModelResponses))/10;
ModelResponses = ModelResponses + 0.1;




%% Ok this should create the same model population as in Batthelier et al

%ie one with firing rates enclosed in the hypercube of the real data, and
%also aligned to its principal components to fit the overall structure of
%the data.


%% Cell Bodies
%load('/Users/Yves/Desktop/20140304_Processed/Area04/(20140304_16_47_27)-_20140304_Area04_tones1_GRABinfo.mat')
pop_response = GRABinfo.NPcSnglTrlRsp;


% %% Boutons
% 
% load('/Users/Yves/Desktop/Boutons_preprocessed/20140116/Area02/(20140116_15_16_17)-_MGB20140116_Area02_tones1_GRABinfo.mat')
% pop_response = GRABinfo.SnglTrlRsp(1:200,:,:);


%Note that is may or may not already have the correct shape. Try both ways.
%If you want to capture the shape of the population responses, then its the
%wrong way around. If you want to capture signle neuron responses, it is
%the right way.
popLin = reshape(pop_response,size(pop_response,1),prod([size(pop_response,2), size(pop_response,3)]))';


%Poplin is arranged so that an observation is single complete population response,
%and a variable is the calcium response of one neuron
[coeff,score,latent,tsquared,explained] = pca(popLin,'Centered','on');

%This calculates orthogonal versions of the PCA coefficients
%coefforth = inv(diag(std(popLin)))*coeff;

%This works out the explained variance, not used any further
VarExplained = cumsum(latent)./sum(latent);

%This should exactly recreate my data, except for the variance and stuff of
%the data
ModelTest = score/coeff;
ModelTest = reshape(ModelTest,size(pop_response,1),size(pop_response,2), size(pop_response,3));


%This gives the correct answer, to within rounding error or something like
%that. SanityCheck should basically be zero
cscores = popLin*coeff - repmat((mean(popLin,1)*coeff),900,1);
SanityCheck = cscores - score;

%Now do the other sanity check, try to recreate popLin from the scores. Ok,
%now this works as well
pLimSim = (score + repmat((mean(popLin,1)*coeff),900,1))/coeff;

%Additionally want to have the variance of single neuron firing. So
%actually match the distribution in real space.

%Ok, now want to find the boarders of the data in PCA space, individually 
%for each cell. Specifically, in doing this, approximate the structure in
%the population responses 


%Rather than below, where the absolute limits are taken, tries to remove
%outliers, so that the resultant population response is still more similar.
sortedScore = sort(score,1,'ascend');
MaxbordersPCA = sortedScore(round(length(sortedScore)*0.98),:);
MinbordersPCA = sortedScore(round(length(sortedScore)*0.02),:);
 
% MaxbordersPCA = max(score,[],1);
% MinbordersPCA = min(score,[],1);

%Actually, I think this is wrong, this is not how to generate the
%pop_responses. No, because it is the coeffs that make sure it mixes
%across. I hope/think matrix division is the same
for ii = 1:size(MaxbordersPCA,2)

ModelResponsesPCA(:,ii) = (MaxbordersPCA(ii)-MinbordersPCA(ii)).*rand(900,1) + MinbordersPCA(ii);

end

%This is the key step!
ModelResp = ((ModelResponsesPCA)/coeff)';

ModelResponses = reshape(ModelResp,size(ModelResp,1),9,100);


%% This generates the simple uniform distirbution model

popLin = reshape(pop_response,size(pop_response,1),prod([size(pop_response,2), size(pop_response,3)]))';

%This is simple absolute max
% MaxbordersUni = max(popLin);
% MinbordersUni = min(popLin);

%This inlcudes an attempt to remove outliers
sortedpopLin = sort(popLin,1,'ascend');
MaxbordersUni = sortedpopLin(round(length(sortedpopLin)*0.95),:);
MinbordersUni = sortedpopLin(round(length(sortedpopLin)*0.02),:);



for ii = 1:size(MaxbordersUni,2)

ModelResponses(ii,:) = (MaxbordersUni(ii)-MinbordersUni(ii)).*rand(900,1) + MinbordersUni(ii);

end


%% Helps visualise this final step
figure,
scatter3(ModelResponses(1,:),ModelResponses(12,:),ModelResponses(14,:),'r')
hold on
scatter3(pop_response(1,:),pop_response(12,:),pop_response(14,:),'b')
hold off

numwindows = round(size(ModelResponses,1) + 0.5);

numrows = round(sqrt(numwindows)+0.5);

figure,

NUMBER = 1;

for ii = 1:numwindows-1
    
    subplot(numrows,numrows,ii), 
    plot(ModelResponses(NUMBER,:),ModelResponses(ii,:),'o','MarkerEdgeColor','b','MarkerSize',6)
    hold on
    plot(pop_response(NUMBER,:),pop_response(ii,:),'.','MarkerEdgeColor','r','MarkerSize',6)
    hold off
    
end

legend('ModelResponses','Population Responses')


