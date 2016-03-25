function [corrMat,reliableSounds,reliableCorrMat,ClustTree,order,clustStore,Clusters,HubertG,HubertG2,FClusterStore,freqOrder,freqLevel] = Cluster_popResponse(popResponse,thresh,reliableOnly,varargin)





if nargin>3
        for field_num=1:(nargin-3)
        
        currentFN = char(varargin(field_num));
            if ~isempty(regexp(currentFN,'toPlot=','ONCE'))
                %find out what sample rate
                sr = regexp(currentFN,'toPlot=','split');
                %and assign it to a variable
                toPlot = str2double(sr(2));
            end
        end
else
    
    toPlot=false;
            
end

     

%% Create the correlation Matrix


[corrMat] = Pop_corr_mat(popResponse);
%% Extract Clusters from the data

[reliableSounds,reliableCorrMat,ClustTree,order,clustStore] = Clusters_from_CorrMat(corrMat,reliableOnly,thresh,toPlot);
 
 

%% Create Cluster Matrix and calculate a Goodness of Fit

[Clusters,HubertG,HubertG2,FClusterStore,freqOrder,freqLevel] = Cluster_Validate(ClustTree,clustStore,order,reliableSounds,toPlot);



end