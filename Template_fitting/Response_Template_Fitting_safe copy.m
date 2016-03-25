%% This script fits a template to the calcium data

function [template_fit,list_idx,response_trace_store,response_trace_criterion_store] = Response_Template_Fitting_safe(traces)

%% Notes
%shorten the template because of short ISIs
%%


load('/Users/Yves/Desktop/TestGRABinfo')

temp = load('/Users/Yves/Documents/MATLAB/G6_exp_fit');
%temp = load('/Users/Yves/Documents/MATLAB/G6_bouton_fit');


global parameter
parameter = temp.a1;


scanrate = 29.6912;

%time in frames
[detection_criterion_p,scale_p,~] = caTransient_detect(traces,0,69);

%Threshold for detection of a response
thresh = 3.5;

%detect the
above_high  = detection_criterion_p>thresh;
cross_high = sum(above_high,2);

above_low = detection_criterion_p>(thresh-0.5);
cross_low = sum(above_low,2);
cross_low(find(cross_high>1),:) = 0;

above_thresh = above_high;
cross_count = cross_high;
%%
%This is the real script

%plot(cross_count,'o');

%There is a real problem with 11
test_NR = 2;%length(find(cross_count>10))-10;


[~, list_idx] = sort(cross_count,'descend');
bouton_idx = list_idx(test_NR);

%most_responsive



%time indices
t = 0:69;
len_template = length(t);
%Define how big the trace segments you cut out will be
pre_event = 50;
post_event = 100;


%% This section is checked and works to find the largest contiguous value



%This selects the best trace
best_trace = traces(bouton_idx,:);

%This finds the length of the trace
len_trace = length(best_trace);

%This is a list of indices, in the indices, of the original trace where the
%criterion is exceeded
above_crit = find(above_thresh(bouton_idx,:));

%This is the separation between the indices
separation = diff(above_crit);
separation


%This is an arbitrary choice of window over which to separate events. It is
%used to to say that if events exceeding the criterion are separated more
%than half_win, then lets call it two separate events.
half_win = 5;

%This returns the indices from the matrix 'separation' that are not
%contiguous
boundaries = [1,find(separation>half_win)+1,length(separation)+1];

%Control Plot
if length(boundaries)<=3
        error('Nothing Crosses The Threshold. Your Trace Sucks :/')
end

%For all windows defined by two boundaries (note in the first line of code
%here use boundary_idx and (boundary_idx+1)
for boundary_idx=1:length(boundaries)-1
    
    %These are the indices of a contiguous set boundary crossings in the
    %indices of the original trace
    contig_crossings = above_crit(boundaries(boundary_idx):boundaries(boundary_idx+1)-1);
    
    %% This is what was in the code originally, don't know why was looking at the gradient
    % now want to take the indices for each event and find the 
    %gradients = conv([1,1,0,-1,-1],best_trace(contig_crossings));
    %[~,d] = max(gradients);
    %%
    if length(contig_crossings)>1
        [~,d] = max(best_trace(contig_crossings));
    else
        d = 1;
    end
    %This is the index of the best crossing
    max_detected(boundary_idx) = d+boundaries(boundary_idx)-1;
end


%This selects the best fitting events in each window.
event_idxs = above_crit(max_detected);

%removes trace that will be out of bound
event_idxs(find(event_idxs<=pre_event+1)) = [];
event_idxs(find(event_idxs>=len_trace-post_event-1-len_template)) = [];

%% Unsure what this bit of code is for, not necessary for the main bit, I think
%[~,separate_crossings] = find_responsive(detection_criterion_p,'single',404);
% separate_crossings = find(sum(reshape(detection_criterion_p(bouton_idx,21:13430),45,13410/45)>3,1));
% stim_responses = reshape(detection_criterion_p(bouton_idx,21:13430),45,13410/45);
% response_traces = stim_responses(:,separate_crossings);
%% 

for ii=1:length(event_idxs)
response_trace_store(ii,:) = traces(bouton_idx,event_idxs(ii)-pre_event:event_idxs(ii)+post_event);

response_trace_criterion_store(ii,:) = max(detection_criterion_p(bouton_idx,event_idxs(ii)-pre_event:event_idxs(ii)+post_event));
end

new_avg_resp = mean(response_trace_store,1);



%% This just plots the traces

stim_times = (1:45:13500)/scanrate;

load('/Users/Yves/Documents/MATLAB/tones_outDat.mat')

color = outDat.trialOrder;

figure()
plot((1:13431)/scanrate,detection_criterion_p(bouton_idx,:),'r');
hold on
plot((1:13500)/scanrate,best_trace,'k');
scatter(stim_times,3*ones(1,length(stim_times)),40,color,'filled')
legend('Detection Criterion','Trace','Tone')
h = colorbar();
%set(h, 'XTick', [0, 50, 100])
%h.Label.String = 'Frequency (kHz)';


xlabel('Time (s)')

title_TEXT = sprintf('Number of Crossings = %d',length(boundaries)-2);
title(title_TEXT)
hold off
%% Fit the function

%need to interpolate
interp_new = interp(new_avg_resp,3);
full_temp_function = @(norm,rise,decay,offset_y,x) norm.*((1-exp(-x/rise)).*exp(-x/decay))+offset_y;



%This is where everything stops working
[template_fit,~,~] = fit([-pre_event:post_event]',new_avg_resp',full_temp_function, ...
                            'StartPoint', [parameter.norm,parameter.rise,parameter.decay,0], ...
                            'Lower',[-Inf,2,2,-Inf],'Upper',[+Inf,15,40,+Inf]);
          


end


%%
%for ii=1:899
%separate_crossings(ii) = length(find(sum(reshape(detection_criterion_p(ii,21:13430),45,13410/45)>3,1)));
%end
