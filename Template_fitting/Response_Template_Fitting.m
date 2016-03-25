%% This script fits a template to the calcium data

function [template_fit,list_idx,AAA] = Response_Template_Fitting(traces)

%% Notes
%shorten the template because of short ISIs
%%


%load('/Users/Yves/Desktop/TestGRABinfo')

temp = load('/Users/Yves/Documents/MATLAB/G6_exp_fit');
%temp = load('/Users/Yves/Documents/MATLAB/G6_bouton_fit');


global parameter
parameter = temp.a1;


%time in frames
detection_criterion_p = caTransient_detect(traces,0);

%Threshold for detection of a response
thresh = 3;


above_thresh  = detection_criterion_p>thresh;
cross_count = sum(above_thresh,2);
%plot(cross_count,'o');

%55 causes problems, for some reason
test_NR = 1;

[~, list_idx] = sort(cross_count,'descend');
best_bouton = list_idx(test_NR);

%most_responsive

%figure out what this is
half_win = floor(121/2);

%time indices
t = 0:69;


%% This section is checked and works to find the largest contiguous value

best_bouton
best_trace = traces(best_bouton,:);

%This is a list of indices where the criterion is exceeded
event_idxs = find(detection_criterion_p(best_bouton,:)>3);

%This is the separation between the indices
separation = diff(event_idxs);

%This returns the indices from the separation matrix that are not
%contiguous
boundaries = [1,find(separation>half_win)+1,length(separation)+2];


for ii=1:length(boundaries)-1
    
    %These are the indices of a contiguous set of values
    contig_crossings = event_idxs(boundaries(ii):boundaries(ii+1)-1);
    
    % now want to take the indices for each event 
    [~,d] = max(best_trace(contig_crossings));
    
    %This is the index of the best event for a series of contiguous events
    max_detected(ii) = d+boundaries(ii)-1;
end

event_idxs = event_idxs(max_detected);


separate_crossings = find(sum(reshape(detection_criterion_p(best_bouton,21:13430),45,13410/45)>3,1));
stim_responses = reshape(detection_criterion_p(best_bouton,21:13430),45,13410/45);
response_traces = stim_responses(:,separate_crossings);

for ii=1:length(event_idxs)
AAA(ii,:) = traces(best_bouton,event_idxs(ii):event_idxs(ii)+69);
end


new_avg_resp = mean(AAA,1);



%% This just plots the traces
figure()
plot(1:13430,detection_criterion_p(best_bouton,:),'r');
hold on
plot(1:13500,best_trace,'k');
plot(0:45:13410,3,'Marker','o','MarkerSize',10, ...
     'MarkerFaceColor','r','MarkerEdgeColor','r')

TEXT = sprintf('Number of Crossings = %d',length(separate_crossings)); %cross_count(list_idx(test_NR)));
title(TEXT)
%% Fit the function



full_temp_function = @(norm,rise,decay,offset,x) norm.*((1-exp(-x/rise)).*exp(-x/decay))+offset;

%This is where everything stops working
[template_fit,b1,c1] = fit([0:length(t)-1]',new_avg_resp',full_temp_function, ...
                            'StartPoint', [parameter.norm,parameter.rise,parameter.decay,0]);
          


end


%%
%for ii=1:899
%separate_crossings(ii) = length(find(sum(reshape(detection_criterion_p(ii,21:13430),45,13410/45)>3,1)));
%end
