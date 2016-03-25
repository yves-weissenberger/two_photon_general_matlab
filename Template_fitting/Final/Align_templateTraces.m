function [aligned_traces,aligned_gradS] = Align_templateTraces(interp_c2,toPlot,varargin)
% Algin_traces returns aligned versions of input traces and optionally
% gradients.
%   'interp_c2' is the vector of unaligned traces
%   'varargin' !must! be used to specify what slope is being aligned to. It
%   should contain an argument 'SampleRate=xx' where xx is the samplerate
%   of the signal in Hz and 'SlopeDuration=yy', where yy is the duration of
%   the slope that is to be aligned to. It is often better to overestimate
%   the duration of the slope as underestimating it may cause noisy
%   estimates of slope position. 
%   For GCaMP6 imaging in boutons 'SlopeDuration=0.4' works well.
%
%
%   Example Align_Traces(interp_c2,'SampleRate=150','SlopeDuration=0.4');

nargin

if nargin>3

    for field_num=1:(nargin-2)
        
        currentFN = char(varargin(field_num));
        %if it finds a field with the name celltype
        if length(regexp(currentFN,'SampleRate='))>0
            
            %find out what sample rate
            sr = regexp(currentFN,'SampleRate=','split');
            %and assign it to a variable
            sample_rate = str2double(sr(2));
        end
        %Same syntax as above   
         if length(regexp(currentFN,'SlopeDuration='))>0
            slp_time = regexp(currentFN,'SlopeDuration=','split');
            slope_duration_s = str2double(slp_time(2));
            
         end

            
     end
    
    
    
else
    
    error('You have not put in enough input arguments')

end



%% This Calculates the gradients of your traces

t = [0:size(interp_c2,2)-1];

len_filter = ceil(slope_duration_s*sample_rate);


%initialise/clear the variables
interp_c2_grad = [];
interp_c2_gradS = [];

for ii=1:size(interp_c2,1)
    interp_c2_grad = conv(interp_c2(ii,:),[ones(1,len_filter),0,(-1)*ones(1,len_filter)],'same');
    interp_c2_gradS(ii,:) =3*smooth(interp_c2_grad,11)./max(interp_c2_grad);
end  


%% First Want to Calculate the Position of Maximal Gradient in the Traces, which is Alignment Point
offset = 20;
[~,max_grad] = max(interp_c2_gradS(:,offset:offset+200),[],2);
%plot(max_grad,'o')

%% Plot the good traces

if strcmp(toPlot,'yesPlot')
    num_traces = size(interp_c2,1);
    num_plots = floor(sqrt(num_traces))^2;

    figure(09)
    for ii=1:num_plots

        subplot(sqrt(num_plots),sqrt(num_plots),ii)
        plot(t(offset:end-100),interp_c2_gradS(ii,offset:end-100),'g')
        hold on
        plot(t(offset:end-100),interp_c2(ii,offset:end-100),'k')
        scatter(offset+max_grad(ii),1,55,'fill','r')
        ylim([-2,4])

        hold off
    end
end
%% Machinery that actually aligns traces accoridng to the point of their Steepest Gradient

%peaks are the point used for alignment

%Find, across traces the latest peak
late_peak = max(max_grad);
%Find, across traces the earliest peak
early_peak = min(max_grad);

%this intermediate result is used to calculate the number of zeros after
%the trace that is required
total_disp = late_peak - early_peak;

%initialise (or clear) variables
aligned_traces =[]; 
noAlign_interp_c2 = [];

pre_zeros = [];
post_zeros = [];

%For all traces
for ii=1:size(interp_c2,1)
    %calculate the number of zeros required to align trace ii with the
    %trace with the latest peak gradient
    pre_zeros = late_peak - max_grad(ii);
    
    %calculate the number of zeros required at the end so that all traces
    %end up the same length.
    post_zeros = total_disp-pre_zeros;
    
    
    %create aligned trace ii, sandwiched between zeros
    aligned_traces(ii,:) = cat(2,zeros(1,pre_zeros),interp_c2(ii,offset:end),zeros(1,post_zeros));
    aligned_gradS(ii,:) = cat(2,zeros(1,pre_zeros),interp_c2_gradS(ii,offset:end),zeros(1,post_zeros));
    %for control purposes/plots, create equal sized matrix of unaligned
    %traces
    noAlign_interp_c2(ii,:) = cat(2,zeros(1,55),interp_c2(ii,offset:end),zeros(1,total_disp-55));
end
   
%% Create diagnostic plots

if strcmp(toPlot,'yesPlot')
    t2 = [0:size(aligned_traces,2)-1];   %
    size(t2)
    size(noAlign_interp_c2)
    figure(07)
    title('Aligned Single Trial Responses used for Averaging to Make Example Trace')
    for ii=1:4

        subplot(sqrt(num_plots),sqrt(num_plots),ii)
        plot(t2,aligned_traces(ii,:),'k')
        hold on

        ylabel('\DeltaF/F','FontSize',10)
        xlabel('Time (s)','FontSize',10)

        ylim([-2,4])

        hold off
    end
end


%Create an Overlay of the single traces in grey as well as plotting the aligned and
%unaligned mean trace.

if strcmp(toPlot,'yesPlot')
    figure(08)

    for ii=1:size(aligned_traces,1) 

        curve10 = plot(aligned_traces(ii,:),'Color',[0.4,0.4,0.4]);
        hold on

    end
        plot(mean(aligned_traces,1),'k','LineWidth',2)
        plot(mean(noAlign_interp_c2,1),'r','LineWidth',2)
        hold off
end


if strcmp(toPlot,'yesPlot')
% Plot aligned and unaligned mean traces as well as a single trace   
    figure(14)
    
    trace_to_plot = round(rand(1)*num_traces);
    plot(t2,mean(aligned_traces,1),'k','LineWidth',2)
    hold on
    plot(t2,mean(noAlign_interp_c2,1),'r','LineWidth',2)
    plot(t2,aligned_traces(trace_to_plot,:),'g','LineWidth',1.5)
    legend('Aligned Average','Unaligned Average','Single Trace')
    title('Effect of Alignment of Single Trial Traces on Mean Response','FontSize',14)
    ylabel('\DeltaF/F','FontSize',13)
    xlabel('Time (s)','FontSize',13)
    hold off
end

end


%% Notes

%% To do List

%Make the beginning and ends not zeros but the same value as the last
%available ones so extrapolate