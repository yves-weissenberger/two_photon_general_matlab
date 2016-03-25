load('/Users/Yves/Desktop/20140304_Processed/Area04/(20140304_16_47_27)-_20140304_Area04_tones1_GRABinfo.mat')


%% Preamble to find reliable Sounds

pop_response = GRABinfo.NPcSnglTrlRsp;


%Remove the erroneously labelled ROIs
ToRemove = find(isnan(GRABinfo.NPcorrectedTraces(:,1)));
pop_response(ToRemove,:,:) = [];

numCells = size(pop_response,1);



for kk=1:size(pop_response,3)
    
    for k=1:size(pop_response,3)
        
        xx = squeeze(pop_response(:,:,kk))';
        y = squeeze(pop_response(:,:,k))';
        
        
        corrmatrix = ((xx*y') - numCells*mean(xx,2)*mean(y,2)')./((numCells-1)*std(xx,0,2)*std(y,0,2)');
        
        mean_store(kk,k) = mean2(corrmatrix);
        
        
        
    end
    
end

    
    D = diag(mean_store);
    
    thresh = 0.2;
    
    reliable_sounds = find(D>=thresh);
    unreliable_sounds = find(D<thresh);
    



%% First, want to restrict FRAs only to the reliable sounds

%This is just for plotting
pop2 = pop_response;
%pop2(:,:,unreliable_sounds) = -0.7;


pop3 = pop2;
pop3(:,:,unreliable_sounds) = [];

%pop3 is a 2-d matrix, with the first dimension being the cell identity and
%the second one the stimulus signifier. Mean responses are contianed in the
%matrix.
pop3 = squeeze(mean(pop3,2));

pop4 = reshape(pop2,size(pop2,1),900);
pop5 = reshape(pop3,size(pop3,1),prod([size(pop3,2), size(pop3,3)]));

%Basically the idea has to do with how many unique representations of
%stimuli there are. 
%Ok, so from simulating random responses on the FRAs show that its not a
%real population pattern. Unless the FRA arises because the population
%pattern in present. Not the case, see the contribution of cortical vs
%thalamus auditory paper.
%Ok, so it doesn't arise because of that. Then arises because of the way
%the neurons are tuned. What property of the distribution of tuning
%properties gives rise to clusters.
%If all neurons respond to all stimuli, then 1 cluster fits the data best.
%If all neurons respond to separate stimuli, then 25 clusters fit the data
%best. If half the neurons respond to one set of stimuli, and the other
%neurons respond to the other set,

%Ok, so in a sense, want to find the n sounds that elicit responses from
%different neurons according to some threshold. More simply, want a way of
%quantifying how different the neuronal populations are that respond to
%each stimulus. Well, simple, plot the populations response to each
%stimulus in the appropriate dimensional space, with the mean response
%being the place where the response is, and look at the distance between
%them.
%ha!, I think this is literally what the clustering algorithm does. With
%the pdist function thing.
%Where do you want to attack, this idea of discrete representation.
%How is this relevant for processing in the brain?

%Maybe try to do PCA though the 100-D positions and look at whether can
%depict in 2/3 dimensions

%No, this is not the same as you did in Model responses, or not quite
%because there, look at the 


%First, just plot out pop3 in PCA space. Each row in pop3 is an observation and
%each column is a variable. So in this case, want to have have the neuronal
%responses as variables, so want to have a 
[coeff,score,latent,tsquared,explained] = pca(pop5');


%each column in coeff corresponds to the coefficients mapping the real data
%onto one principal component.

%Want to plot the mean population response to each stimulus in numCells
%dimensional space. No, just plot the mean response in numCell dimensional
%space.
%Want to see whether can find any 2 or three stimuli that are represented
%by unique sets of neurons. Should be able to plot a difference FRA for
%each sound.
%Take the population response to one stimulus and subtract it from the
%population response to another stimulus.
%Using this metric, if two sounds are represented by entirely different
%populations.
%Actually, first want to normalise the entire response matrix (i.e. pop3) to between -1
%and 1.

pop3norm = (pop3 - (max(max(pop3)) - abs(min(min(pop3))))/2);
pop3norm = pop3norm./max(max(abs(pop3norm)));


%Now work out the difference in the population response to each of the
%sounds.

for ii = 1:(size(pop3norm,2)-1)        
    corrPops(ii) = mean(pop3norm(:,1) - pop3norm(:,ii+1));
end

%Its literally exactly what they did. Its all about finding responses that
%are similar to one one another but different to other ones. That is
%literally the simplest way of finding clusters.

