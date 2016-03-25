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


responsiveBoutons = [12,91]; %[39,159];

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



%% Create an anonymous function of the template

analyse_trace = GRABinfo.Traces(91,:);


t = 0:69;

template =  a.norm.*(1-exp(-t/rise)).*exp(-t/decay);
len_template = length(template);

%% Go through the trace and find the scaling factor etc for each event


for ii=1:(length(analyse_trace)-len_template)
    
    data = analyse_trace(ii:ii+69);
    N = length(template);

    scale(ii) = (sum(template.*data) - (sum(template)*sum(data/N)))/ ...
                (sum(template)^2 - (sum(template)*sum(template/N)));
            
    offset(ii) = (sum(data) - scale(ii)*sum(template))/N;
    
    %Standard Error of Fit
%      SSE(ii) = sum(data.^2) + scale(ii)^2 * sum(template.^2) + N + offset(ii)^2 ...
%            -2*(scale(ii)*sum(template.*data) + offset(ii)*sum(data) -  ...
%             scale(ii)*offset(ii)*sum(template));
    
       
    SSE(ii) = sum((data - (template*scale(ii) + offset(ii))).^2);
       
    standard_error(ii) = (SSE(ii)/(N-1))^(1/2);
    
end


detection_criterion = scale./standard_error;

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




%% Plots the entire event trace together with the 

plot(analyse_trace(1:13430),'b');
hold on;
plot(1:13430,detection_criterion,'k','LineWidth',1.5);
plot(0:13430,2,'g','LineWidth',2)
plot(0:13430,1.5,'r','LineWidth',2)



%% This section is checked and works to find the largest contiguous value

%This is a list of indices where the criterion is exceeded
event_idxs = find(detection_criterion>2);

%This is the separation between the indices
separation = diff(event_idxs);

%This returns the indices from the separation matrix that are not
%contiguous
boundaries = [1,find(separation>half_win)+1,length(separation)+2];

for ii=1:length(boundaries)-1
    
    %These are the indices of a contiguous set of values
    contig_crossings = event_idxs(boundaries(ii):boundaries(ii+1)-1);
    
    % now want to take the indices for each event 
    [a,b] = max(detection_criterion(contig_crossings));
    
    %This is the index of the best event for a series of contiguous events
    max_detected(ii) = b+boundaries(ii)-1;
end

event_idxs = event_idxs(max_detected);

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

new_avg_resp = mean(AAA,1);

%This is where everything stops working
[a1,b1,c1] = fit([0:length(average_response)-1]',new_avg_resp',temp_fnctn, ...
              'Algorithm','Levenberg-Marquardt','StartPoint', [a.norm,a.rise,a.decay]);

plot(a1,[0:length(average_response)-1],average_response)


