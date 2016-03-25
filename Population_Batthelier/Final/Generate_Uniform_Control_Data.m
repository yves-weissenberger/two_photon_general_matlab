%% This generates the simple uniform distirbution model
function [ModelResponses] = Generate_Uniform_Control_Data(pop_response)


popLin = reshape(pop_response,size(pop_response,1),prod([size(pop_response,2), size(pop_response,3)]))';

%This is simple absolute max
% MaxbordersUni = max(popLin);
% MinbordersUni = min(popLin);

%This inlcudes an attempt to remove outliers
sortedpopLin = sort(popLin,1,'ascend');
MaxbordersUni = sortedpopLin(round(length(sortedpopLin)*0.95),:);
MinbordersUni = sortedpopLin(round(length(sortedpopLin)*0.02),:);


%initialise the matrix
ModelResponses = zeros(size(pop_response,1),size(pop_response,2),size(pop_response,3));

for ii = 1:size(MaxbordersUni,2)

    ModelResponses(ii,:,:) = (MaxbordersUni(ii)-MinbordersUni(ii)).*rand(9,100) + MinbordersUni(ii);

end


end