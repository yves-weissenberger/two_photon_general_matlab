function [ModelResponses] = Generate_PCA_Control_Data(pop_response)



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

%ModelResponsesPCA(:,ii) = (MaxbordersPCA(ii)-MinbordersPCA(ii)).*normrnd(0,0.2,900,1) + MinbordersPCA(ii);
ModelResponsesPCA(:,ii) = (MaxbordersPCA(ii)-MinbordersPCA(ii)).*rand(900,1) + MinbordersPCA(ii);


end

%This is the key step!
ModelResp = ((ModelResponsesPCA)/coeff)';

ModelResponses = reshape(ModelResp,size(ModelResp,1),9,100);

end