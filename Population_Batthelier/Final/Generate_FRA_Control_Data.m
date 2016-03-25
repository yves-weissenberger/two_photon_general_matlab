function [ModelResponse] = Generate_FRA_Control_Data(popResponse,toPlot)


%% First extract distribution Parameters from the Raw Data

SnglNeuronMean = squeeze(mean(popResponse,2));
SnglNeuronStd = squeeze(std(popResponse,[],2));


NeuronsStdStore = reshape(SnglNeuronStd,numel(SnglNeuronStd),1);
NeuronsMeanStore = reshape(SnglNeuronMean,numel(SnglNeuronMean),1);

[param,paramConfInt,residuals,residualsConfInt,stats] = regress(NeuronsStdStore,[ones(size(NeuronsMeanStore,1),1),NeuronsMeanStore]);

stdMeanFit =  param(1) + param(2)*NeuronsMeanStore;

if toPlot==true
    %from plot, looks like the residuals are pretty linear
    figure, 

    hold on
    plot(NeuronsMeanStore,stdMeanFit,'r')
    plot(NeuronsMeanStore,NeuronsStdStore,'o')

    plot(NeuronsMeanStore,stdMeanFit,'r')
    xlabel('Mean Response of Each Neuron to each stimulus')
    ylabel('Std of the responses of each neuron')
    legend('Line of Best Fit','Single Neuron, single Stimulus Mean & STD')
    title('Relation of Mean to STD of Responses in the population')
    hold off
end


%Ok, so since there is a linear relationship between the Std of neural
%responses and their mean


%%  Get the FRAs and do stuff

[fraMeanStore,fraStdStore,pVals,fraMeanAll,fraStdAll] = getFRAs(popResponse);

numCells = size(popResponse,1);

numFreqs = size(fraMeanStore,3);
numLevels =  size(fraStdStore,2);

numStim = size(fraMeanAll,2);


ModelResponse = zeros(size(popResponse));


for cellIdx = 1:numCells
    
    for stimIdx = 1:numStim
        
        ModelResponse(cellIdx,:,stimIdx) = normrnd(fraMeanAll(cellIdx,stimIdx),fraStdAll(cellIdx,stimIdx),[1,9]);

    end
    
end








end






















%% This is the old code, something simpler is implemented above

% %intialise the 
% response_store = zeros(1,9,100);
% 
% %First, get Std of responses independent of the mean response
% for cellIdx = 1:numCells
%     
%     
%     
%     for stimIdx = 1:numStim
%         
%         
%         %response_store(1,:,jj) = (pop_response(ii,:,jj))/(mean(pop_response(ii,:,jj))*0.63);
%         
% 
%         response_store(1,:,stimIdx) = (popResponse(cellIdx,:,stimIdx)*100- mean(popResponse(cellIdx,:,stimIdx))*100)/(mean(popResponse(cellIdx,:,stimIdx))*0.63*100);
%         test1(:,stimIdx) = (popResponse(cellIdx,:,stimIdx)*100- mean(popResponse(cellIdx,:,stimIdx))*100)/(mean(popResponse(cellIdx,:,stimIdx))*0.63*100);
%         %This is now a 9x100 variable of repetitions in rows and stim in
%         %columns
%         Response_Store = squeeze(response_store);
%         
%         
%         %Ok, so this is the order of response size in the population
%         [Value OrderRespStore] = sort(mean(popResponse(cellIdx,:,:),2),3,'ascend');
%         
%         
%         %This should order the responses, such that the mean response in
%         %column 1 is the smallest etc
%         Response_Store = Response_Store(:,OrderRespStore);
%         
%         
%         
%     end
%     
%     
%     
%     
%     
%     RespStLinear2(cellIdx,:) = reshape(Response_Store(:,70:100),1,prod(size(Response_Store(:,70:100))));
%     
% end
% 
% 
% VarHistData = reshape(RespStLinear2,1,prod(size(RespStLinear2)));
% 
% 
% clear pd xvalues pdf
% 
% 
% 
% [nelements,centers] = hist(VarHistData,60);
% 
% bar(centers,nelements/sum(nelements));
% 
% 
% % hold on
% % 
% % pd = fitdist(VarHistData','Normal');
% % xvalues = [-15:1:15];
% % pdf = pdf(pd,xvalues);
% % plot(xvalues',pdf,'LineWidth',2)
% % 
% % hold off
% % 
% % 
% 
% 
% %% Fit the Variability to key distribution Parameters
% 
% moments = {mean(VarHistData),std(VarHistData),skewness(VarHistData),kurtosis(VarHistData)};
% %rng('default');  % For reproducibility
% [r,type] = pearsrnd(moments{:},1000,1);
% figure, hist(r,30)
% 
% 
% 
% %% Fit Smart Model
% 
% ModelResponses = zeros(numCells,9,100);
% 
% %NOT SURE WHY, BUT THESE PERMUTATIONS DON'T SEEM TO WORK
% sFRAall = permute(reshape(smFRA_store,100,size(smFRA_store,3)),[2 1]);
% %FRAall = permute(reshape(FRAall,100,size(FRAall,3)),[2 1]);
% 
% 
%     %pop_response = GRABinfo.NPcSnglTrlRsp;
% 
%     FRAStore= squeeze(mean(popResponse,2));
%     
%     
%     
% 
% for cellIdx = 1:size(sFRAall,1)
%     
%     %Selects the FRA of a given cell
%     CellFRA = smFRA_store(:,:,cellIdx);
%     
%     %reshapes this FRA to be linear
%     CellLinFRA = FRAStore(cellIdx,:);
%     %CellLinFRA = reshape(permute(CellFRA,[2 1]),1,100);
%     
%     
% 
% 
%     
%     %Run through all the stimuli
%     for stimIdx=1:size(CellLinFRA,2)
%         
%         
%         %This is supposed to multiply the mean responses of a given cell to a given tone by a series of random numbers 
%         ModelResponses(cellIdx,:,stimIdx) = pearsrnd(CellLinFRA(stimIdx),CellLinFRA(stimIdx)+0.3,moments{3},moments{4},1,9,1);
% 
%         
%     end
%     
% end
% 
% 
% ModelResponses(find(isnan(ModelResponses))) = 0;
% 
