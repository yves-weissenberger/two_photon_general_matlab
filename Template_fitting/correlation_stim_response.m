%% This Script outputs the distribution of correlations between time post stimulus onset and ?F/F
  

function [corrz,percentile,std_deviation] = correlation_stim_response(traces,ISI)


num_ROIs = size(traces,1);

%This is the number of stimuli presented
num_stim = 13500/ISI;

%This generates an array in which the 1st index is bouton_ID, the second
%index is a given post-stimulus trace and the third is an index of which
%post-stimulus trace is in question
ordered_traces = reshape(traces,num_ROIs,ISI,num_stim);

x_vals = 1:45;


corrz = zeros(num_ROIs,num_stim);


for jj=1:num_ROIs

    if jj==1
        tic;
    elseif jj==2
        one_loop = toc;
        expected_time = round(one_loop*num_ROIs);
        sprintf('Expected Duration is %d seconds',uint8(expected_time))

    end
    
   
    
    for ii=1:num_stim


        temp = corrcoef(ordered_traces(jj,:,ii),x_vals);
        corrz(jj,ii) = temp(2,1);

    end


    
end

idx_offset = 580;
num_plots = 49;
std_deviation = std(corrz,1,2);
%NB low kurtosis means more data is in the tails
kurt_data = kurtosis(corrz,1,2);

scatter(std_deviation,kurt_data,8,'MarkerEdgeColor','g')
xlabel('STD')
ylabel('Kurtosis')
title('STD vs Kurtosis of the Histogram of Correlations between Trial by Trial post-Stimulus Time and\DeltaF/F')


%The reason this isn't working is because array_order tells you the
%the which element should go into the current position, rather than given
%an index in the old one, which index will it go to in the sorted one.
[~,temp_order] = sort(std_deviation,'descend');
[~,array_order] = sort(temp_order,'ascend');

percentile = 1-((array_order)./length(array_order));


for ii=1:num_plots
    

    subplot(sqrt(num_plots),sqrt(num_plots),ii)
    over_title= sprintf('B%d SD = %.2f; %.1f percentile ',uint16(ii + idx_offset),std_deviation(ii+idx_offset),100*percentile(ii+idx_offset));
    hist(corrz(idx_offset+ii,:),5)
    title(over_title);
    xlim([-0.8,0.8])
    ylim([0,100])
end

%This just plot os SDs
figure, plot(std(corrz,1,2),'o')






%% This part of the scipt plots a histogram of the correlation values.
min_SD = min(std_deviation);
max_SD = max(std_deviation);

num_bars = 30;

Hy_vals = hist(std_deviation,num_bars);

Hx_vals= linspace(min_SD,max_SD,num_bars);

figure(21), bar(Hx_vals,Hy_vals./sum(Hy_vals));


%% Maybe try to do PCA of the histograms

%First need to create all the histograms

clear hist_store
%in coeff, the coefficients for 1 principal component are along a column
for ii=1:size(corrz,1)

    hist_store(ii,:) = hist(corrz(ii,:),10);

end


[coeff,score,latent] =  pca(hist_store);

%each column of score corresponds to the projections of all observations
%onto one pca axis
figure(20),
scatter3(score(:,1),score(:,2),score(:,3),'o','MarkerFaceColor','k')
xlabel('PC1')
ylabel('PC2')

frac_var_explained =cumsum(latent/sum(latent));

num_princomp = length(latent);
figure(24), 
bar(0:num_princomp,[0,frac_var_explained'])
ylabel('Cumulative Variance Explained')
xlabel('Principal Components Included')
xlim([0,10])
ylim([0,1])

end


%% Notes %%

%Mean of correlation coefficients of a time series is not the same as 
%correlation of the mean time series



%Just by observarion and double checking, it looks like traces with low
%kurtosis are real responsive traces, but that might not be true

% ab = rand(1,5)
%
% ab =
% 
%     0.1731    0.5174    0.9953    0.7076    0.0806
% 
% [p,q]=sort(ab,'descend')
% 
% p =
% 
%     0.9953    0.7076    0.5174    0.1731    0.0806
% 
% 
% q =
% 
%      3     4     2     1     5
% 
% desired = [4,3,1,2,5]