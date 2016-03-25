% This function interpolates and then fits to multiple multiple multiple of
% the interpolated values

upSR_rate = 5;
scanrate = 29.6912;

len_orig_trace = size(c2,2);

new_avg_resp = mean(c2,1);

interp_resp = interp(new_avg_resp,upSR_rate);
t = [0:length(interp_resp)-1];%/(GRABinfo.scanFrameRate*upSR_rate);
num_responses = size(c2,1);

%This was just found by trial and error. Maybe should justify...
len_filter = 30;

good_resp = [];

%threshold for considering a good response
good_resp = find(d2>=8);

loop_var1 = [];
loop_var1(1,:) = good_resp;
loop_var1(2,:) = 1:size(good_resp,2);

%initialise/clear the variables
interp_c2 = [];
interp_c2_grad = [];
interp_c2_gradS = [];

for ii=loop_var1
    interp_c2(ii(2),:) = interp(c2(ii(1),:),upSR_rate);
    interp_c2_grad(ii(2),:) = conv(interp_c2(ii(2),:),[ones(1,len_filter),0,(-1)*ones(1,len_filter)],'same');
    interp_c2_gradS(ii(2),:) =3*smooth(interp_c2_grad(ii(2),:),11)./max(interp_c2_grad(ii(2),:));


end  

%% This Detects the double peaked responses

%
grad_thresh = interp_c2_gradS>1;
cross_grad_thresh = diff(grad_thresh,1,2);

num_good_resp = size(good_resp,2);

for ii=1:num_good_resp
num_crosses(ii) = max(length(find(cross_grad_thresh(ii,:)==1)),length(find((cross_grad_thresh(ii,:).*-1)==1)));
end

%This variable is one if it only crosses oncev and 0 else
single_cross=num_crosses~=1;




for ii=1:22
    subplot(6,6,ii)
    plot(interp_c2(ii,:),'k')
    hold on
    plot(interp_c2_gradS(ii,:),'g')
    plot([0,800],[1,1],'--b')
    ylim([-2,5])

    if single_cross(ii)==1
        plot([0,800],[5,-2],'r','LineWidth',5)
    end

end


%% This removes double peaks
interp_c2(single_cross,:) = [];
interp_c2_grad(single_cross,:) = [];
interp_c2_gradS(single_cross,:) = [];

num_traces = size(interp_c2,1);



%Looks like it succesfully detects double peaked response
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

%%
offset = 100;
[~,max_grad] = max(interp_c2_gradS(:,offset:offset+200),[],2);
%plot(max_grad,'o')

%% Plot the good traces
num_plots = floor(sqrt(num_traces))^2;

figure(09)
for ii=1:num_plots
    
    subplot(sqrt(num_plots),sqrt(num_plots),ii)
    plot(t(offset:end-100),interp_c2_gradS(ii,offset:end-100),'g')
    hold on
    plot(t(offset:end-100),interp_c2(ii,offset:end-100),'k')
    scatter(offset+max_grad(ii),1,55,'fill','r')
    title(d2(good_resp(ii)))
    ylim([-2,4])

    hold off
end

%% Align the good single trace responses

%peaks are the point used for alignment
late_peak = max(max_grad);
early_peak = min(max_grad);
total_disp = late_peak - early_peak;


align_interp_c2 =[]; %= zeros(size(interp_c2,1),335);%(late_index-early_index+size(interp_c2(:,offset:end-100),2)));

pre_zeros = [];
post_zeros = [];

for ii=1:size(interp_c2,1)
    pre_zeros = late_peak - max_grad(ii);
    post_zeros = total_disp-pre_zeros;
    align_interp_c2(ii,:) = cat(2,zeros(1,pre_zeros),interp_c2(ii,offset:end),zeros(1,post_zeros));
    noAlign_interp_c2 = cat(2,zeros(1,55),interp_c2(ii,offset:end),zeros(1,total_disp-55));
    %length(cat(2,zeros(1,pre_zeros),zeros(1,post_zeros)))
end
    
%%

t2 = [0:size(align_interp_c2,2)-1]/(upSR_rate*scanrate);

figure(07)
title('Aligned Single Trial Responses used for Averaging to Make Example Trace')
for ii=1:4
    
    subplot(sqrt(num_plots),sqrt(num_plots),ii)
    plot(t2,align_interp_c2(ii,:),'k')
    hold on

    ylabel('\DeltaF/F','FontSize',10)
    xlabel('Time (s)','FontSize',10)

    title(d2(good_resp(ii)) )
    ylim([-2,4])

    hold off
end



figure(08)

for ii=1:num_traces 
    
    curve10 = plot(align_interp_c2(ii,:),'Color',[0.4,0.4,0.4]);
    hold on
    
end
    plot(mean(align_interp_c2,1),'k','LineWidth',2)
    plot(mean(noAlign_interp_c2,1),'r','LineWidth',2)
    hold off

    figure(14)
    trace_to_plot = round(rand(1)*num_traces);
    plot(t2,mean(align_interp_c2,1),'k','LineWidth',2)
    hold on
    plot(t2,mean(noAlign_interp_c2,1),'r','LineWidth',2)
    plot(t2,align_interp_c2(trace_to_plot,:),'g','LineWidth',1.5)
    legend('Aligned Average','Unaligned Average','Single Trace')
    title('Effect of Alignment of Single Trial Traces on Mean Response','FontSize',14)
    ylabel('\DeltaF/F','FontSize',13)
    xlabel('Time (s)','FontSize',13)
    hold off
    

%% This finds the range of values around which to search for a start point


fit_trace = mean(align_interp_c2,1);
len_filter_2 = 20;

aligned_mean_grad = conv(fit_trace,[ones(1,len_filter_2),0,(-1)*ones(1,len_filter_2)],'same');


search_start_range = [50,5];
[~,start_point] = max(aligned_mean_grad);
search_start_point = start_point-search_start_range(1):start_point+search_start_range(2);
first_search_idx = start_point-search_start_range(1);
%% plot to show that find the correct start point in the data
figure(13)
plot(fit_trace)
hold on
plot(aligned_mean_grad/max(aligned_mean_grad));
plot(search_start_point,zeros(1,length(search_start_point)),'g','LineWidth',4)
legend('Trace','Gradient of Trace','Search for Startpoint')
hold off


%% Find the best start point
tmax = length(align_interp_c2)/(upSR_rate*scanrate);
length_t = length(align_interp_c2);

t = linspace(0,tmax,length_t);

full_temp_function = @(norm,rise,decay,offset_y,x) norm.*((1-exp(-x/rise)).*exp(-x/decay))+offset_y;

best_finish = length_t;

fit_error = [];


for ii=search_start_point
    
    try
        [a3,b3,c3] = fit((t(ii:best_finish)-t(ii))',fit_trace(ii:best_finish)',full_temp_function, ...
                     'StartPoint', [4,0.3,1,fit_trace(ii)-0.2], 'Lower',[-1000,-1000,-1000,fit_trace(ii)-0.2], ...
                     'Upper',[+1000,+1000,+1000,fit_trace(ii)+0.2]);
    catch
        a3store{ii-first_search_idx+1} = a3;
        start_fit_error(ii-first_search_idx+1) = +inf;
        ii
    end
        
    a3store{ii-first_search_idx+1} = a3;
    start_fit_error(ii-first_search_idx+1) = b3.rmse/length(t(ii:best_finish));%b3.sse/length(t(ii:end));
end

%%
[~,low_fin_Error] = min(start_fit_error);
%plot(search_start_point, start_fit_error)

best_start = low_fin_Error + first_search_idx;

%% Does Control plots for best start point
want_subplots=0;

%best_start=121;
if want_subplots==1
pre_start = best_start-20;

aga = 1;
if aga ==1
for ii=2:2:40
    subplot(5,4,ii/2)
    plot((t(pre_start:best_finish)-t(pre_start))',fit_trace(pre_start:best_finish)','g*')
    hold on
    plot((t(pre_start+ii:best_finish)-t(pre_start))',fit_trace(pre_start+ii:best_finish)','b+')
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
end
% %% Fit the best end point
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

%%
best_fit = best_start;
best_finish = length_t;


[a4,b4,c4] = fit((t(best_fit:best_finish)-t(best_fit))',mean(align_interp_c2(:,best_fit:best_finish),1)',full_temp_function, ...
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


plot((t(best_fit-80:best_finish)-t(best_fit)),mean(align_interp_c2(:,best_fit-80:best_finish),1),'*g')
hold on
plot((t(best_fit:best_finish)-t(best_fit)),mean(align_interp_c2(:,best_fit:best_finish),1),'*r')
plot((t(best_fit:best_finish)-t(best_fit)),template)
set(gca,'FontSize',12)
text(2,1,params_print,'FontSize',12)
xlabel('Time (s)','FontSize',14)
ylabel('\DeltaF/F','FontSize',14)
h=legend('Full Data','Fit Data','Fit Curve');
hold off

%On this plot want the REAL non interpolated data as well, to see how well
%it fits that.


%% Bugs

%Problems with 
