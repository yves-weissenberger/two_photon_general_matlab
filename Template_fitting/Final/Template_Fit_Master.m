clear all

%Load the set of data code is tested on
load('/Users/Yves/Documents/MATLAB/two_photon_general_matlab/Template_fitting/TestGRABinfo')

%This loads the order of the tones
load('/Users/Yves/Documents/MATLAB/two_photon_general_matlab/Template_fitting/tones_outDat.mat')


%This just loads an initial estimate of the time course of the fit
temp = load('/Users/Yves/Documents/MATLAB/two_photon_general_matlab/Template_fitting/G6_exp_fit.mat');

global parameter
parameter = temp.a1;
clear temp

temp2 = load('/Users/Yves/Documents/MATLAB/two_photon_general_matlab/Template_fitting/plotColors.mat');
global plotColors
plotColors = temp2.plotColors;
clear temp2





%% Lets try low pass filtering everything first

%These are the traces from the code
traces = GRABinfo.Traces;

scanrate = GRABinfo.scanFrameRate;

Nqst = scanrate/2;


% bands is a vector of pairs of normalized frequency points, specified in the
% range between 0 and 1, where 1 corresponds to the Nyquist frequency. The 
%frequencies must be in increasing order.
%cutoff_freq = 8*scan

%just picked something that works.
bands = [0,0.55,0.7];  %/scanRate;
bands(end+1) = 1;
ampl = [1,1,0,0];


%low_pass = firpm(122,bands,ampl);
[b,a]=butter(2,.1); 

for bouton_num = 1:size(traces,1)
    temp0 = filter(b,a,traces(bouton_num,:));
    %temp0 = upfirdn(traces(bouton_num,:),low_pass,1,1);
    traces(bouton_num,:) = temp0(61:end-62);
end


%% Set some parameters


%This is not the actual idx of the bouton, but the bouton with the n-most
%threshold crossings. To obtain the actual idx of the bouton, in order to
%look at it from traces use bouton_idx = b2(bouton_NR);
bouton_NR = 5;

 len_template = 65;
 extraction_threshold = 3.5;

%extraction_threshold = 3.5;
%len_template = 12;


%% Fit the template to the trace and obtain the detection and scale criteria
tic
[detection_criterion_p,scale_p,standard_error_p] = caTransient_detect(traces,0,len_template);
toc
%% This makes sure that each detected event is only considered once.

tic
%To plot things change the final variable to 'yesPlot'. To not plot, any other string
[template_fit,b2,response_trace_store,response_trace_criterion_store,response_trace_scale_store,event_times] = Event_extraction(traces,detection_criterion_p,scale_p,standard_error_p,extraction_threshold,len_template,bouton_NR,'noPlot');
toc
%% This fits the interpolated traces

%if this step returns an error message, it is most likely that nothing fits
%the template

tic 
inclusion_thresh = 5;

[param_fit] = fit_interpTraces(response_trace_store,response_trace_criterion_store,inclusion_thresh,template_fit,'yesPlot');
toc

%% Recreate the trace from the event times and the anonymous function and look at residuals


full_temp_function = @(norm,rise,decay,offset_y,x) norm.*((1-exp(-x/rise)).*exp(-x/decay))+offset_y;

len_trans = 200;
t = [0:len_trans];


Transient_shape = full_temp_function(parameter.norm,parameter.rise,parameter.decay,0,t);

%This should be the noise reduced trace
fit_trace = zeros(1,size(traces,2));

for event_nr = 1:size(response_trace_store,1)
    event_idx = event_times(event_nr);
    if event_idx<length(fit_trace)-len_trans
        fit_trace(event_idx:event_idx+len_trans) = fit_trace(event_idx:event_idx+len_trans) + Transient_shape.*response_trace_scale_store(event_nr);
        %fit_trace(event_idx:event_idx+len_trans) = Transient_shape.*response_trace_scale_store(event_nr);

    end
end



%%


stim_times = (1:45:13500)/scanrate;

color = log(outDat.trialOrder);

fig = figure('Units','normalized');
tt = [1:13500]/GRABinfo.scanFrameRate;

bouton_idx = b2(bouton_NR);
hold on

plot(tt,traces(bouton_idx,:),'ro-')
plot(tt,fit_trace,'b','LineWidth',2)
xlabel('Time (s)')
ylabel('\Delta F/F')

ll = legend('Data','Fit');
ll.Position = [0.6,0.83,0.2,0.05];
scatter(stim_times,2*ones(1,length(stim_times)),40,color,'filled')

r = xcorr(Transient_shape(1:40),traces(3,:));
plot(tt,r(1:13500)./30,'g','LineWidth',2)
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
grid();
h = colorbar();
h.Ticks = linspace(0,max(color),5);
h.TickLabels = round(logspace(0,log10(80),5));

hold off


%% Best segmentation

normV1 = parameter.norm*5;
normV2 = parameter.norm;
normV3 = parameter.norm;


Tshape1 = full_temp_function(normV1,parameter.rise,parameter.decay,0,t);
Tshape2 = full_temp_function(normV2,parameter.rise,parameter.decay,0,t);
Tshape3 = full_temp_function(normV3,parameter.rise,parameter.decay,0,t);



t1 = cat(2,Tshape1,zeros(1,90));
t2 = cat(2,zeros(1,45),Tshape2,zeros(1,45));
t3 = cat(2,zeros(1,90),Tshape3);



for ii = 1:100

    for trl = 1:9
        %From inspection noise looks like it goes up and down by 0.2DF/F
        %noise = (rand(3,length(t1))-0.5)./2.5;
        noise = normrnd(0,0.2,3,length(t1));




        transDiff = (t2+noise(1,:)) - ((t1+noise(2,:))+(t3+noise(3,:)));

        %transDiff(transDiff<=0) = 0;


        store(trl,:) = transDiff;

    end


    meanTrace = mean(store,1);

    firstcrossTEMP1 = find(meanTrace(20:end)>0.2);
    firstCROSS = firstcrossTEMP1(1)+20; 

    secondcrossTEMP1 = find(meanTrace(60:end)<0.2);
    secondCROSS = round(mean(secondcrossTEMP1(1)+secondcrossTEMP1(2)))+60;
    
    cumsumAUC(ii,:) = cumsum(meanTrace);
    
    store1(ii) = firstCROSS;
    store2(ii) = secondCROSS;

end


%% Create a histogram of first and second crossings

max1 = max(store1);
min1 = min(store1);
bins1 = min1:1:max1;

max2 = max(store2);
min2 = min(store2);
bins2 = min2:1:max2;

binnedCross1 = histcounts(store1,bins1);
binnedCross2 = histcounts(store2,bins2);


clrs = plotColors{2}./255;
subplot(2,1,1)
bar1 = bar(bins1(2:end),binnedCross1./sum(binnedCross1));
bar1.FaceColor = clrs(1,:);
bar1.EdgeColor = 'w';

subplot(2,1,2)
bar2 = bar(bins2(2:end),binnedCross2./sum(binnedCross2));
bar2.FaceColor = clrs(2,:);
bar2.EdgeColor = 'w';



%% Work out the AUC for various combinations of things

AUC = [];
for strt=45:90
    
    for windowLen = 1:60
        
        AUC(strt-44,windowLen,:) = sum(cumsumAUC(:,strt+windowLen-1) - cumsumAUC(:,strt));
        
    end
    
    
end


[value, location] = max(AUC(:));

[bestStartFrame,bestWindowSize] = ind2sub(size(AUC),location);


im1 = imagesc(AUC);
xlabel('Number of Stimulus frames');
ylabel('Starting Frame (after Stimulus Onset)');
set(gca,'FontWeight','bold');
set(gca,'FontSize',12);
%% This plots the raw data, just a visualisation

load('/Users/Yves/Documents/MATLAB/plotColors')
clrs = plotColors{3}./255;

subplot(2,1,1)
hold on
plot(t1,'LineWidth',2,'Color',clrs(1,:),'LineStyle','--')
plot(t1+noise(1,:),'LineWidth',1.5,'Color',clrs(1,:))
plot(t2,'LineWidth',2,'Color',clrs(2,:),'LineStyle','--')
plot(t2+noise(2,:),'LineWidth',2,'Color',clrs(2,:))

plot(t3,'LineWidth',2,'Color',clrs(3,:),'LineStyle','--')
plot(t3+noise(3,:),'LineWidth',2,'Color',clrs(3,:))

plot(t1+t2+t3+noise(1,:),'LineWidth',1.5,'Color','k')
legend('mean transient 1','transient 1','mean transient 2','transient 2', ...
       'mean transient 3','transient 3','Simulated trace')

grid()
hold off


subplot(2,1,2)
plot(meanTrace,'LineWidth',2,'Color','k')
ylim([-1,1])
legend('Mean Difference Between Noisy Second and other Noisy Traces')
grid()




%From looking at this plot, want to go from 55 to 90

%% This just loads the colors used for plotting

figure()




colors = plotColors{2}./255;

residuals = GRABinfo.Traces(bouton_idx,:) - fit_trace;

meanRs = mean(residuals);
sigmaRs = std(residuals);
bins = [-4:0.05:5];
binned_residuals = histcounts(residuals,bins);

hBars = bar(bins(1:end-1),binned_residuals/sum(binned_residuals),1,'histc');
hBars.FaceAlpha = 1;
hBars.FaceColor = colors(1,:);
hBars.EdgeColor = 'w';



hold on
gauss_eqn = @(meanRs,sigmaRs,x_vals) (1./(sigmaRs.*sqrt(2*pi))).*exp(-((x_vals - meanRs).^2)/(2*(sigmaRs.^2)));

gauss_resids = gauss_eqn(meanRs,sigmaRs,bins);
ylabel('Frequency %')
xlabel('Residual Size')

plot(bins,gauss_resids./sum(gauss_resids),'Color',colors(2,:),'LineWidth',2)



grid()


%% 

figure()
trace1 = traces(b2(2),:);
meanT1 = mean(trace1);
stdT1 = std(trace1);

trace2 = traces(b2(20),:);
meanT2 = mean(trace2);
sdtT2 = std(trace2);

bins2 = [-5:0.05:5];
binned_vals = histcounts(trace2,bins2);

binned_vals((length(bins)+3)/2:length(bins)) = binned_vals((length(bins)+3)/2:length(bins)) - fliplr(binned_vals(1:(length(bins)-1)/2));

hBars = bar(bins2(1:end-1),binned_vals/sum(binned_vals),1,'histc');
hBars.FaceAlpha = 1;
hBars.FaceColor = colors(1,:);
hBars.EdgeColor = 'w';


