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
detection_criterion_p = caTransient_detect(traces,0,69);

%Threshold for detection of a response
thresh = 2;

%% Interlude to visualise crossings above/below certain things


detection_criterion_p(:,1:28) = 0;
detection_criterion_p(:,end-28:end) = 0;

for ii=1:size(detection_criterion_p,1)
    
    high_idxs = find(detection_criterion_p(ii,:)>(thresh+0.5));
    
    for jj=-25:25
        detection_criterion_p(ii,high_idxs+jj) = 0;
    end
        
end
%%
%This is the real script
above_thresh  = detection_criterion_p>thresh;
cross_count = sum(above_thresh,2);
max(cross_count)
%plot(cross_count,'o');

%There is a real problem with 19
test_NR = 1; %max(find(cross_count>5))

[~, list_idx] = sort(cross_count,'descend');
best_bouton = list_idx(test_NR);

%most_responsive

%figure out what this is
half_win = floor(121/2);

%time indices
t = 0:69;

%Define how big the trace segments you cut out will be
pre_event = 50;
post_event = 100;


%% This section is checked and works to find the largest contiguous value

best_trace = traces(best_bouton,:);
len_trace = length(best_trace);

%This is a list of indices where the criterion is exceeded
event_idxs = find(detection_criterion_p(best_bouton,:)>thresh);

%This is the separation between the indices
separation = diff(event_idxs);

%This returns the indices from the separation matrix that are not
%contiguous
boundaries = [1,find(separation>half_win)+1,length(separation)+2];

if length(boundaries)<=3
        error('Nothing Crosses The Threshold. Your Trace Sucks :/')
end

for ii=1:length(boundaries)-1
    
    %These are the indices of a contiguous set of values
    contig_crossings = event_idxs(boundaries(ii):boundaries(ii+1)-1);
    
    % now want to take the indices for each event 
    gradients = conv([1,1,0,-1,-1],best_trace(contig_crossings));
    [~,d] = max(gradients);
    %This is the index of the best crossing
    max_detected(ii) = d+boundaries(ii)-1;
end

event_idxs = event_idxs(max_detected);

%removes trace that will be out of bound
event_idxs(find(event_idxs<=pre_event+1)) = [];
event_idxs(find(event_idxs>=len_trace-post_event-1)) = [];


%[~,separate_crossings] = find_responsive(detection_criterion_p,'single',404);

separate_crossings = find(sum(reshape(detection_criterion_p(best_bouton,21:13430),45,13410/45)>3,1));
stim_responses = reshape(detection_criterion_p(best_bouton,21:13430),45,13410/45);
response_traces = stim_responses(:,separate_crossings);


for ii=1:length(event_idxs)
response_trace_store(ii,:) = traces(best_bouton,event_idxs(ii)-pre_event:event_idxs(ii)+post_event);
response_trace_criterion_store(ii) = max(detection_criterion_p(best_bouton,event_idxs(ii)-pre_event:event_idxs(ii)+post_event));
end

new_avg_resp = mean(response_trace_store,1);



%% This just plots the traces

stim_times = (1:45:13500)/scanrate;

load('/Users/Yves/Documents/MATLAB/tones_outDat.mat')

color = outDat.trialOrder;

figure()
plot((1:13431)/scanrate,detection_criterion_p(best_bouton,:),'r');
hold on
plot((1:13500)/scanrate,best_trace,'k');
scatter(stim_times,3*ones(1,length(stim_times)),40,color,'filled')
legend('Detection Criterion','Trace','Tone')
h= colorbar();
%set(h, 'XTick', [0, 50, 100])
%h.Label.String = 'Frequency (kHz)';


xlabel('Time (s)')

TEXT = sprintf('Number of Crossings = %d',length(separate_crossings));
title(TEXT)
%% Fit the function

%need to interpolate
interp_new = interp(new_avg_resp,3);
full_temp_function = @(norm,rise,decay,offset_y,x) norm.*((1-exp(-x/rise)).*exp(-x/decay))+offset_y;



%This is where everything stops working
[template_fit,b1,c1] = fit([0:150]',new_avg_resp',full_temp_function, ...
                            'StartPoint', [parameter.norm,parameter.rise,parameter.decay,0], ...
                            'Lower',[-Inf,2,2,-Inf],'Upper',[+Inf,15,40,+Inf]);
          


end


%%
%for ii=1:899
%separate_crossings(ii) = length(find(sum(reshape(detection_criterion_p(ii,21:13430),45,13410/45)>3,1)));
%end
