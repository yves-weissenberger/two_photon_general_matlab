function [params_fit] = fit_interpTraces(response_trace_store,response_trace_criterion_store,threshold,template_fit,toPlot)

% This function interpolates and then fits to multiple multiple multiple of
% the interpolated values

upSR_rate = 5;

%This is the threshold for inclusion
%threshold = 5;
scanrate = 29.6912;


len_orig_trace = size(response_trace_store,2);

new_avg_resp = mean(response_trace_store,1);

interp_resp = interp(new_avg_resp,upSR_rate);
t = [0:length(interp_resp)-1]; %/(GRABinfo.scanFrameRate*upSR_rate);
num_responses = size(response_trace_store,1);

%This was just found by trial and error. Maybe should justify...
len_filter = 30;

good_resp = [];

%threshold for considering a good response
good_resp = find(response_trace_criterion_store>=threshold)';

loop_var1 = [];
loop_var1(1,:) = good_resp;
loop_var1(2,:) = 1:size(good_resp,2);

%initialise/clear the variables
interp_c2 = [];


for ii=loop_var1
    interp_c2(ii(2),:) = interp(response_trace_store(ii(1),:),upSR_rate);
end  

%% Align the interpolated traces


[aligned_traces,aligned_gradS] = Align_templateTraces(interp_c2,toPlot,'SampleRate=150','SlopeDuration=0.4');

%% This Detects the double peaked responses

grdnt_thresh = 1;

grad_thresh = aligned_gradS>grdnt_thresh;
cross_grad_thresh = diff(grad_thresh,1,2);

num_good_resp = size(good_resp,2);

for resp_idx=1:num_good_resp
num_crosses(resp_idx) = max(length(find(cross_grad_thresh(resp_idx,:)==1)),length(find((cross_grad_thresh(resp_idx,:).*-1)==1)));
end

%This variable is one if it only crosses once and 0 else
double_cross=num_crosses~=1;

check_late = sum(aligned_gradS(:,end/2:end)>grdnt_thresh,2);
late_cross = check_late~=0;

num_traces = size(aligned_traces,1);
num_rows = ceil(sqrt(num_traces));

if strcmp(toPlot,'yesPlot')
    for trace_idx=1:num_traces
        subplot(num_rows,num_rows,trace_idx)
        plot(aligned_traces(trace_idx,:),'k')
        hold on
        plot(aligned_gradS(trace_idx,:),'g')
        plot([0,800],[1,1],'--b')
        ylim([-2,5])

        if double_cross(trace_idx) || late_cross(trace_idx)==1
            plot([0,800],[5,-2],'r','LineWidth',1)
        end

    end
end

to_Remove = (double_cross + late_cross')>=1;
%% This removes double peaks
aligned_traces(to_Remove,:) = [];
aligned_gradS(to_Remove,:) = [];


%% This finds the range of values around which to search for a start point


aligned_mean = mean(aligned_traces,1);
len_filter_2 = 20;

aligned_mean_grad = conv(aligned_mean,[ones(1,len_filter_2),0,(-1)*ones(1,len_filter_2)],'same');


search_start_range = [50,5];
[~,start_point] = max(aligned_mean_grad);
search_start_point = start_point-search_start_range(1):start_point+search_start_range(2);
first_search_idx = start_point-search_start_range(1);

% Old Position of section: "plot to show that find the correct start point
% in the data"
%% Find some Parameters
tmax = length(aligned_mean)/(upSR_rate*scanrate);
length_t = length(aligned_mean);

t = linspace(0,tmax,length_t);

full_temp_function = @(norm,rise,decay,offset_y,x) norm.*((1-exp(-x/rise)).*exp(-x/decay))+offset_y;

best_finish = length_t;

fit_error = [];


%% Find the best start point

for ii=search_start_point
    
    try
        [a3,b3,c3] = fit((t(ii:best_finish)-t(ii))',aligned_mean(ii:best_finish)',full_temp_function, ...
                     'StartPoint', [4,0.3,1,aligned_mean(ii)-0.2], 'Lower',[-1000,-1000,-1000,aligned_mean(ii)-0.2], ...
                     'Upper',[+1000,+1000,+1000,aligned_mean(ii)+0.2]);
    catch
        a3store{ii-first_search_idx+1} = template_fit;
        start_fit_error(ii-first_search_idx+1) = +inf;
    end
        
    a3store{ii-first_search_idx+1} = a3;
    start_fit_error(ii-first_search_idx+1) = b3.rmse/length(t(ii:best_finish));%b3.sse/length(t(ii:end));
end

%Find the index of the lowest error relative to the indices searched
[~,low_fin_Error] = min(start_fit_error);
%plot(search_start_point, start_fit_error)

%Find the index of the lowest error relative to the indices of the entire
%trace
best_start = low_fin_Error + first_search_idx;

%% Does Control plots for best start point

if strcmp(toPlot,'yesPlot')
pre_start = best_start-20;

    for ii=2:2:40
        subplot(5,4,ii/2)
        plot((t(pre_start:best_finish)-t(pre_start))',aligned_mean(pre_start:best_finish)','g*')
        hold on
        plot((t(pre_start+ii:best_finish)-t(pre_start))',aligned_mean(pre_start+ii:best_finish)','b+')
        %create plot of the fit curve
        param_bstart = a3store{ii};
        template_bstart = full_temp_function(param_bstart.norm,param_bstart.rise,param_bstart.decay,param_bstart.offset_y,(t(pre_start+ii:best_finish)-t(pre_start+ii)));
        plot((t(pre_start+ii:best_finish)-t(pre_start))',template_bstart','r');
        hold off
        xlabel('Time (s)')
        ylabel('\DeltaF/F')
        legend('Full Data Set','Fit Data Set','Fit Curve')

        title(sprintf('Start Point = %d',pre_start+ii))
        if (pre_start+ii)== low_fin_Error + first_search_idx
            text(0.5,0,'Best Fit','FontSize',16,'Color','r')
        end
    end
end


%% Provide Final Fit
best_fit = best_start;
best_finish = length_t;


[a4,b4,c4] = fit((t(best_fit:best_finish)-t(best_fit))', aligned_mean(best_fit:best_finish)',full_temp_function, ...
                'StartPoint', [4,0.3,1,0], 'Lower',[-Inf,-Inf,-Inf,0],'Upper',[+Inf,+Inf,+Inf,0.5]);

%best start point in non_downsampled coordinates
orig_best = floor(best_fit/upSR_rate);
%plot(a4)


%%

template = full_temp_function(a4.norm,a4.rise,a4.decay,a4.offset_y,(t(best_fit:best_finish)-t(best_fit)));
[~,peak_idx] =max(template);

params = [a4.norm,t(peak_idx),a4.decay,a4.offset_y,b4.rsquare];
params_print = sprintf(' Norm = %.3g \n Rise Time, t_{peak} (s) = %.3g \n Decay Time, t_{1/2} (s) = %.3g \n Offset = %.3g \n r^{2} = %.3g',params);

figure(22),


plot((t(best_fit-80:best_finish)-t(best_fit)),aligned_mean(:,best_fit-80:best_finish),'*g')
hold on
plot((t(best_fit:best_finish)-t(best_fit)),aligned_mean(:,best_fit:best_finish),'*r')
plot((t(best_fit:best_finish)-t(best_fit)),template)
set(gca,'FontSize',12)
text(2,1,params_print,'FontSize',12)
xlabel('Time (s)','FontSize',14)
ylabel('\DeltaF/F','FontSize',14)
h=legend('Full Data','Fit Data','Fit Curve');
hold off

%On this plot want the REAL non interpolated data as well, to see how well
%it fits that.


%% Now want to shift the parameters back into the space used for fitting
upSR_rate = 5;
scanrate = 29.6912;

params_fit = {};
params_fit.norm = a4.norm;
params_fit.rise = a4.rise*scanrate;
params_fit.decay = a4.decay*scanrate;
params_fit.offset = a4.offset_y;


%% Notes


%% Additional Code

%% bar plot
% interp_c2_gradS_reshape = reshape(interp_c2_gradS,[1,prod(size(interp_c2_gradS))]);
% interp_c2_gradS_reshape(find(interp_c2_gradS_reshape>0)) = [];
% interp_c2_gradS_reshape = abs(interp_c2_gradS_reshape);
% 
% num_bars = 50;
% 
% hist_gradS = hist(interp_c2_gradS_reshape,num_bars);
% total_gradS = sum(hist_gradS);
% norm_hist_gradS = hist_gradS/total_gradS;
% x_vals = linspace(min(interp_c2_gradS_reshape),max(interp_c2_gradS_reshape),num_bars);
% 
% outliers = find(cumsum(norm_hist_gradS)>0.95);
% outlier_pos = outliers(1);
% 
% figure(17)
% plot([outlier_pos,outlier_pos],[0,max(norm_hist_gradS)],'r','LineWidth',2)
% hold on
% bar(x_vals,norm_hist_gradS,'histc')
% xlim([0,outlier_pos+2])
% hold off

%% Fit the best end point
% 
% end_fit_error = [];
% for ii=600:630
%     [a3,b3,c3] = fit((t(best_start:ii)-t(ii))',fit_trace(best_start:ii)',full_temp_function, ...
%                     'StartPoint', [4,0.3,1,fit_trace(ii)-0.2], 'Lower',[-1000,-1000,-1000,-1000], ...
%                      'Upper',[+1000,+1000,+1000,fit_trace(ii)+0.2]);
%                  
%     a3store{ii-159} = a3;
%     end_fit_error(ii-599) = b3.rmse/length(t(best_start:ii));%b3.sse/length(t(ii:end));
%     ii
% end
% 
% 
% plot(600:630,end_fit_error)
% 
% [~,low_fin_Error] = min(end_fit_error);
% best_finish = low_fin_Error + 600;


%% plot to show that find the correct start point in the data
% figure(13)
% plot(aligned_mean)
% hold on
% plot(aligned_mean_grad/max(aligned_mean_grad));
% plot(search_start_point,zeros(1,length(search_start_point)),'g','LineWidth',4)
% legend('Trace','Gradient of Trace','Search for Startpoint')
% hold off

%% Bugs

%Problems with 
