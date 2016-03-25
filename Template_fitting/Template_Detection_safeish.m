%% This script detects events in calcium traces using the method described by


%       "Detection of Spontaneous Synaptic Events with an Optimally Scaled Template"
%       by Clements and Bekkers (1997)

%Following additional procedures described in 

%       "Learning-related fine-scale specificity imaged in motor cortex circuits of behaving mice"
        %Komiyama et al (2010)

% "For event detection (Supplementary Fig. 9), we detected segments of the 
% trace that exceeded the baseline by twice the standard deviation of the 
% baseline (?) for at least two successive frames. ? was calculated from 
% selecting values below baseline and appending values of the flipped sign,
% to construct an estimate of the baseline fluctuations without contaminations
% from true calcium signals"

%% Initialise Parameters



% This is the data file to load
%load('/Users/Yves/Desktop/2p_MGB_imaging/Boutons_preprocessed/20140116/Area01/(20140116_14_54_50)-_MGB20140116_Area01_tones2_GRABinfo.mat')
load('/Users/Yves/Desktop/2p_MGB_imaging/Boutons_preprocessed/20140121/Area05/(20140121_15_53_29)-_MGB20140121_Area05_tones1_GRABinfo.mat')


scanrate= GRABinfo.scanFrameRate;

%Bouton 334 has to be one of the best
responsiveBoutons = [582,334,191,91,166];

%Length of the window
len_window = 121;
half_win = floor(len_window/2);


%% 


%now detect all the events that exceed the threshold
responsive_trace = GRABinfo.Traces(responsiveBoutons(2),:);
len_Trace = size(GRABinfo.Traces,2);

%regions with large activity
high_dF = find(responsive_trace>2);

separation = diff(high_dF);

close_idx = find(separation<half_win) + 1;

high_dF(close_idx) = [];


high_dF(find(high_dF<34)) = [];
high_dF(find(high_dF>len_Trace-34)) = [];



%Number of responses crudely detected
numResp = length(high_dF);




responses = zeros(numResp,len_window);



%now find the maxima in each of these regions
for kk=1:numResp  
    responses(kk,:)= responsive_trace(high_dF(kk)-half_win:high_dF(kk)+half_win);    
end

average_response = mean(responses,1);



%This was found my inspection, would want to automate it via the gradient
%thing
average_response = average_response(52:end);


%% This then fits the function

%This is from the GCaMP6 paper
time_to_peak = 0.3*scanrate;
decay_time = 1*scanrate;
init_amp = 2;

%reverse
rise=time_to_peak;
decay = decay_time;

%these are the parameters
init_params = [init_amp,rise,decay];



%equation describing the template 
%template = norm*(1-exp(-t/rise)*exp(-t/decay));

%anonymouse function handle to the template
%making decay smaller makes the decay time faster same with rise
%temp_fnctn = @(norm,rise,decay,x) norm.*(1-exp(-x/rise)).*exp(-x/decay);

temp_fnctn = @(norm,x) norm.*(1-exp(-x/rise)).*exp(-x/decay);

%(norm,rise,decay,x)
initial_values = [init_amp,time_to_peak,decay_time];

%This is where everything stops working
[a,b,c] = fit([0:length(average_response)-1]',average_response',temp_fnctn, ...
              'Algorithm','Levenberg-Marquardt','StartPoint', init_amp);

plot(a,[0:length(average_response)-1],average_response)


%The template time course parameters are estimated from  few large synaptic
%events that are detected manually.


time_around = 2;


%% Create an anonymous function of the template

for iii=1:825

% Loop across all boutons
tic
analyse_trace = GRABinfo.Traces(iii,:);


t = 0:69;
if time_around==1
    template =  a.norm.*(1-exp(-t/rise)).*exp(-t/decay);
else
    template =  a1.norm.*(1-exp(-t/a1.rise)).*exp(-t/a1.decay);
end


len_template = length(template);

%% Go through the trace and find the scaling factor etc for each event


for ii=1:(length(analyse_trace)-len_template)
    
    data = analyse_trace(ii:ii+69);
    N = length(template);

    scale(ii) = (sum(template.*data) - (sum(template)*sum(data/N)))/ ...
                (sum(template.^2) - (sum(template)*sum(template/N)));
            
    offset(ii) = (sum(data) - scale(ii)*sum(template))/N;
    
    %Standard Error of Fit
%      SSE(ii) = sum(data.^2) + scale(ii)^2 * sum(template.^2) + N + offset(ii)^2 ...
%            -2*(scale(ii)*sum(template.*data) + offset(ii)*sum(data) -  ...
%             scale(ii)*offset(ii)*sum(template));
    
       
    SSE(ii) = sum((data - (template*scale(ii) + offset(ii))).^2);
       
    standard_error(ii) = (SSE(ii)/(N-1))^(1/2);
    
end


detection_criterion = scale./standard_error;

cross3_store(iii) = length(find(detection_criterion>=2.5));

end
%% A Succesful attempt to make the above parallel

tic
analyse_trace = GRABinfo.Traces(iii,:);


t = 0:69;
if time_around==1
    template =  a.norm.*(1-exp(-t/rise)).*exp(-t/decay);
else
    template =  a1.norm.*(1-exp(-t/a1.rise)).*exp(-t/a1.decay);
end


N = length(template);

num_boutons = size(GRABinfo.Traces,1);


template_parallel = repmat(template,num_boutons,1);


%initialise the storage vectors
scale_parallel = zeros(num_boutons,length(analyse_trace)-N);
offset_parallel = zeros(num_boutons,length(analyse_trace)-N);
SSE_parallel = zeros(num_boutons,length(analyse_trace)-N);
standard_error_p = zeros(num_boutons,length(analyse_trace)-N);




    scale(ii) = (sum(template.*data) - (sum(template)*sum(data/N)))/ ...
                (sum(template.^2) - (sum(template)*sum(template/N)));
            

for ii=1:(length(analyse_trace)-N)
    
    data_parallel = GRABinfo.Traces(:,ii:ii+69);

    scale_parallel(:,ii) = (sum(template_parallel.*data_parallel,2) - (sum(template_parallel,2).*sum(data_parallel/N,2)))./ ...
                (sum(template_parallel.^2,2) - (sum(template,2).*sum(template_parallel/N,2)));
            
    offset_parallel(:,ii) = (sum(data_parallel,2) - scale_parallel(:,ii).*sum(template_parallel,2))./N;
    

    SSE_parallel(:,ii) = sum((data_parallel - (template_parallel.*repmat(scale_parallel(:,ii),1,N) + repmat(offset_parallel(:,ii),1,N))).^2,2);
       
    
end


    standard_error_p = (SSE_parallel./(N-1)).^(1/2);
    detection_criterion_p = scale_parallel./standard_error_p;
    

%% Ok, now want a plot of the number of 

%thresh = 2.5;


for ii = 0:100
    
    thresh = ii/20;

    above_thresh  = detection_criterion_p>thresh;
    cross_count = sum(above_thresh,2);

    count_responsive(ii+1) = length(find(cross_count));
    

end
%% Plot the number of boutons crossing each threshold.
gradient_cur = conv(count_responsive,[-1,0,1],'same');
gradient_cur(find(gradient_cur<0)) = 0;
gradient_cur(end) = 0;

figure(77),
h1 = axes;


curve1 = plot(5:-1/20:0,fliplr(count_responsive));
hold on


curve2 = plot(5:-1/20:0,fliplr(gradient_cur),'g');

legend('# Responsive Cells','gradient')


xlabel('Threshold')
%ylabel('# Responsive Cells')
title('Number of Responsive Cells at Varying Thresholds')
set(h1,'Xdir','reverse')

%So don't forget, want to find the mean gradient. Maybe just smooth it over
%or something like that.

%Ok, look at the distribution of the number of threshold crossings. Or that
%might not be informative, because whe have a responsive cell, in a lot of
%cases, not all the responses cross the threshold.

%%
onwards=1;

if onwards==1
%%

% %% Double Check!!
% detection_criterion(find(detection_criterion==inf)) = 0;
% detection_criterion(find(detection_criterion==-inf)) = 0;
% 
% detection_criterion = 5*detection_criterion/max(detection_criterion);


%% Diagnostic Plots

time = [1:70]./scanrate;

%This is the index of the position
start_pos = 4700;

plot_show=1;

if plot_show==1
    figure,
    for ii=1:25

        pos_idx = start_pos+ii*5;


        title(int2str(pos_idx))
        subplot(5,5,ii)

    %     plot(time,analyse_trace(pos_idx:pos_idx+69),'o','MarkerSize',4, ...
    %          'MarkerFaceColor','k','MarkerEdgeColor','k','Marker','*')


        plot(time,analyse_trace(pos_idx:pos_idx+69),'b')
        ylim([-1,1.5])
        set(gca,'YTick',[-0.75:0.75:1.5])



        hold on
        plot(time,template,'k')
        plot(time,template*scale(pos_idx),'r','LineWidth' ...
             ,1.5,'LineStyle','-');


        %This just plots a horizontal line at the rise time
        plot(rise/GRABinfo.scanFrameRate,[-1:1.5],'LineStyle','-','Color','g')


        if ii==25
            legend('Raw Data','Unscaled Template','Fitted Template')
        end
        xlabel('Time (s)','FontSize',14)
        ylabel('\DeltaF/F','FontSize',14)
        grid()

    end



    %% Test Plots to see what an event is
    figure(45),
    plot(GRABinfo.Traces(334,:),'b');
    hold on;
    plot(1:13430,detection_criterion_p(334,:),'k','LineWidth',1.5);
    plot(0:45:13410,3,'Marker','o','MarkerSize',10, ...
        'MarkerFaceColor','r','MarkerEdgeColor','r')
    plot(0:13430,2,'g','LineWidth',2)
    plot(0:13430,1.5,'r','LineWidth',2)
end


%% This section is checked and works to find the largest contiguous value

%This is a list of indices where the criterion is exceeded
event_idxs = find(detection_criterion>3);

%This is the separation between the indices
separation = diff(event_idxs);

%This returns the indices from the separation matrix that are not
%contiguous
boundaries = [1,find(separation>half_win)+1,length(separation)+2];

for ii=1:length(boundaries)-1
    
    %These are the indices of a contiguous set of values
    contig_crossings = event_idxs(boundaries(ii):boundaries(ii+1)-1);
    
    % now want to take the indices for each event 
    [c,d] = max(detection_criterion(contig_crossings));
    
    %This is the index of the best event for a series of contiguous events
    max_detected(ii) = d+boundaries(ii)-1;
end

event_idxs = event_idxs(max_detected);
toc
%%

num_plots = min(16,2*round(length(event_idxs)-0.5));

if num_plots == 16
    num_draw = num_plots;
    row_cols = [4,4];
    
else
    
    if num_plots<=7
       num_draw= num_plots;
       row_cols = [1,num_draw];
    elseif num_plots==8
       num_draw= num_plots;
       row_cols = [4,2];


    elseif num_plots ==12||13
       num_draw=10;
       row_cols=factor(num_draw);

    elseif isprime(num_plots)
        num_draw=num_plots-1;
        row_cols=factor(num_draw);

    else
        num_draw=num_plots;
        row_cols = factor(num_plots);
    
    end
end

num_draw=16;

for ii=1:num_draw
    
        subplot(row_cols(1),row_cols(2),ii)
            
        plot(-10:59,analyse_trace(event_idxs(ii)-10:event_idxs(ii)+59),'b')
        hold on
        plot(0:69,template*scale(event_idxs(ii))+offset(event_idxs(ii)),'k');
        
        plot(-10,-2:0.1:2,'r-')
        plot(0,-2:0.1:2,'g-')

        hold off
        
        
        AAA(ii,:) = analyse_trace(event_idxs(ii):event_idxs(ii)+69);
end

 
%%
% 
new_avg_resp = mean(AAA,1);


full_temp_function = @(norm,rise,decay,offset,x) norm.*((1-exp(-x/rise)).*exp(-x/decay))+offset;

%This is where everything stops working
[a1,b1,c1] = fit([0:length(average_response)-1]',new_avg_resp',full_temp_function, ...
              'StartPoint', [init_params,0]);

figure,
plot(a1,[0:length(average_response)-1],new_avg_resp)

time_around =2;

end
%% Notes %%


%The problem with just using the correlation between post-stimulus time and
%activity is that neuropil correction sucks and there are lots of motion
%artifacts



% How to organise the data? Want the ultimate output to be.


% Think of the curve fitting for bouton vs axons as well. Look at decay vs
% rise constants or the mean traces from each.