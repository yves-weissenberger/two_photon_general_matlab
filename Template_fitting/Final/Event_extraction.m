%% This script fits a template to the calcium data

function [template_fit,list_idx,response_trace_store,response_trace_criterion_store,response_trace_scale_store,event_idxs] = Event_extraction(traces,detection_criterion_p,scale_p,standard_error_p,extraction_threshold,len_template,test_NR,toPlot)

%% Notes
%This script basically makes sure that one only takes one instance of each
%event to the further analysis. This is necessary because each event will
%exceed the detection criterion at multiple time points

%shorten the template because of short ISIs. This doesn't work. 
%%


load('/Users/Yves/Documents/MATLAB/two_photon_general_matlab/Template_fitting/TestGRABinfo')

temp = load('/Users/Yves/Documents/MATLAB/two_photon_general_matlab/Template_fitting/G6_exp_fit.mat');
global parameter
parameter = temp.a1;


scanrate = 29.6912;



%detect the
above_high  = detection_criterion_p>extraction_threshold;
cross_high = sum(above_high,2);

above_low = detection_criterion_p>(extraction_threshold-0.5);
cross_low = sum(above_low,2);
cross_low(find(cross_high>1),:) = 0;

above_thresh = above_high;
cross_count = cross_high;
%%
%This is the real script

%plot(cross_count,'o');

%There is a real problem with 11
%test_NR = 11;%length(find(cross_count>10))-10;


[~, list_idx] = sort(cross_count,'descend');
bouton_idx = list_idx(test_NR);

%most_responsive



%time indices
t = 0:len_template;
%Define how big the trace segments you cut out will be
pre_event = 35;
post_event = 35;


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


%This is an arbitrary choice of window over which to separate events. It is
%used to to say that if events exceeding the criterion are separated more
%than half_win, then lets call it two separate events.
half_win = 10;

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
    
    %%
    if length(contig_crossings)>1
        %THIS STEP NEEDS IMPROVAL, WHAT IS GOING ON HERE                                              %THIS IS GOING WRONG!!! 
        %[~,d] = max(best_trace(contig_crossings));
        [~,d] =min(detection_criterion_p(bouton_idx,contig_crossings));
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


%% 

for ii=1:length(event_idxs)
response_trace_store(ii,:) = traces(bouton_idx,event_idxs(ii)-pre_event:event_idxs(ii)+post_event);

response_trace_criterion_store(ii,:) = max(detection_criterion_p(bouton_idx,event_idxs(ii)-pre_event:event_idxs(ii)+post_event));
response_trace_scale_store(ii,:) = max(scale_p(bouton_idx,event_idxs(ii)-pre_event:event_idxs(ii)+post_event));
end

new_avg_resp = mean(response_trace_store,1);



%% This just plots the traces

stim_times = (1:45:13500)/scanrate;

load('/Users/Yves/Documents/MATLAB/tones_outDat.mat')
color = outDat.trialOrder;

if strcmp(toPlot,'yesPlot')
    figure()
    plot((1:length(detection_criterion_p))/scanrate,detection_criterion_p(bouton_idx,:),'r');
    hold on
    plot((1:length(best_trace))/scanrate,best_trace,'k');
    scatter(stim_times,3*ones(1,length(stim_times)),40,color,'filled')
    legend('Detection Criterion','Trace','Tone')
    h = colorbar();
    %set(h, 'XTick', [0, 50, 100])
    %h.Label.String = 'Frequency (kHz)';


    xlabel('Time (s)')

    title_TEXT = sprintf('Number of Crossings = %d',length(boundaries)-2);
    title(title_TEXT)
    hold off
end
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
