function [ModelResponses] = Generate_Control_Data(popResponse,varargin)


if nargin>1
    for field_num=1:(nargin-1)
        
        currentFN = char(varargin(field_num));
        %if it finds a field with the name celltype
        if length(regexp(currentFN,'controlType='))>0
            %find out what sample rate
            sr = regexp(currentFN,'controlType=','split');
            %and assign it to a variable
            controlType = sr(2);
        end
        
        %same as above
        if ~isempty(regexp(currentFN,'toPlot='))
            %find out what sample rate
            sr = regexp(currentFN,'toPlot=','split');
            %and assign it to a variable
            toPlot = str2double(sr(2));
        end
       
            
    end
     
else
    controlType = 'PCA';
    toPlot = true;
    
end


if strcmp(controlType,'PCA')
    [ModelResponses] = Generate_PCA_Control_Data(popResponse);
elseif strcmp(controlType,'Uniform')
    [ModelResponses] = Generate_Uniform_Control_Data(popResponse);
elseif strcmp(controlType,'FRA')
    [ModelResponses] = Generate_FRA_Control_Data(popResponse,toPlot);
end




if toPlot==true
    
    
%% Visualise this final step
  

    numWindows = round(size(ModelResponses,1) + 0.5);

    numRows = round(sqrt(numWindows)+0.5);

    figure('Units','normalized')

    baseNeuronIdx = 1;
    
    %This plots, trial by trial the relation between the response
    %amplitudes of the different neurons
    for subPlotNr = 1:numWindows-1

        subplot(numRows,numRows,subPlotNr), 
        plot(ModelResponses(baseNeuronIdx,:),ModelResponses(subPlotNr,:),'o','MarkerEdgeColor','b','MarkerSize',6)
        hold on
        plot(popResponse(baseNeuronIdx,:),popResponse(subPlotNr,:),'.','MarkerEdgeColor','r','MarkerSize',6)
        hold off
        if subPlotNr==floor(numRows/2)+1
            title('Each Boxplot is the Trial by Trial Relation between the Firing Rate of one Neuron and of Another Neuron from the sample')
        end

    end

    legend('ModelResponses','Population Responses')

end