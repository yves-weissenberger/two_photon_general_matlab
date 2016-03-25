%% This is a function that returns the event detection criteria for the time course of arbitrary number of boutons


%NB: in traces, each row is the time course of 1 bouton
function [detection_criterion_p,scale_p,standard_error_p] = caTransient_detect(traces,sensitivity_plot,varargin)


%This is the fit for boutons with GCaMP6m
temp = load('/Users/Yves/Documents/MATLAB/G6_exp_fit');
%temp = load('/Users/Yves/Documents/MATLAB/G6_bouton_fit');

global parameter
parameter = temp.a1;

clear temp
if exist('varargin','var')
    len_template = varargin{1};
else
    len_template = 69;
end

global t
t = 0:len_template;

num_boutons = size(traces,1);


%% Extect




template =  parameter.norm.*(1-exp(-t/parameter.rise)).*exp(-t/parameter.decay) + parameter.offset;


N = length(template);

trace_length = size(traces,2);

template_parallel = repmat(template,num_boutons,1);


%Initialise the vectors to store data
scale_p = zeros(num_boutons,trace_length-N);
offset_p = zeros(num_boutons,trace_length-N);
SSE_p = zeros(num_boutons,trace_length-N);
standard_error_p = zeros(num_boutons,trace_length-N);
            
%Fill the vector with data
for ii=1:(trace_length-len_template)
    
    data_p = traces(:,ii:ii+len_template);

    scale_p(:,ii) = (sum(template_parallel.*data_p,2) - (sum(template_parallel,2).*sum(data_p/N,2)))./ ...
                (sum(template_parallel.^2,2) - (sum(template,2).*sum(template_parallel/N,2)));
            
    offset_p(:,ii) = (sum(data_p,2) - scale_p(:,ii).*sum(template_parallel,2))./N;
    

    SSE_p(:,ii) = sum((data_p - (template_parallel.*repmat(scale_p(:,ii),1,N) + repmat(offset_p(:,ii),1,N))).^2,2);
       
    
end

scale_p(find(scale_p<0.25)) = 0;

standard_error_p = (SSE_p./(N-1)).^(1/2);

detection_criterion_p = scale_p./standard_error_p;

%% Plot the number of boutons crossing each threshold.

if sensitivity_plot==1
    
    max_crosses = 300;
    
    for jj=1:max_crosses
        
        for ii = 0:100

        thresh = ii/20;

        above_thresh  = detection_criterion_p>thresh;
        cross_count = sum(above_thresh,2);

        count_responsive(jj,ii+1) = length(find(cross_count>=jj));


        end
        
    end
    
    gradient_cur = conv(count_responsive(1,:),[-1,0,1],'same');
    gradient_cur(find(gradient_cur<0)) = 0;
    gradient_cur(end) = 0;

%% First plot    
    figure(77),
    h1 = axes;
    
    curve1 = plot(5:-1/20:0,fliplr(count_responsive(1,:)));
    hold on

    curve2 = plot(5:-1/20:0,fliplr(gradient_cur),'g');
    
    curve3 = plot(5:-1/20:0,num_boutons*ones(1,length(5:-1/20:0)),'--r');

    legend('# Responsive Cells','gradient','Number of ROIs')


    xlabel('Threshold')
    %ylabel('# Responsive Cells')
    title('Number of Responsive Cells at Varying Thresholds')
    set(h1,'Xdir','reverse')
    hold off
%% Second (more informative plot)
    color = zeros(max_crosses,3);
    color(:,3) = linspace(0,1,max_crosses);
    color(:,1) = linspace(1,0,max_crosses);
    figure(4)
    h2 = axes;
    for jj=1:max_crosses
        plot(5:-1/20:0,fliplr(count_responsive(jj,:)),'Color',color(jj,:));
        hold on
    end
    h_clr = colorbar();
    clr_map_Sz = size(colormap());
    clr_map = zeros(clr_map_Sz);
    clr_map(:,3) = linspace(0,1,64);
    clr_map(:,1) = linspace(1,0,64);
    colormap(clr_map);
    caxis([1,max_crosses]);
    h_clr.Label.String = 'Number of Threshold Crossings to be Counted Responsive';
    
    %ylabel('# Responsive Cells')
    title('Number of Responsive Cells at Varying Thresholds')
    set(h2,'Xdir','reverse')
    xlabel('Threshold')

    hold off

    
    %So don't forget, want to find the mean gradient. Maybe just smooth it over
    %or something like that.

    %Ok, look at the distribution of the number of threshold crossings. Or that
    %might not be informative, because whe have a responsive cell, in a lot of
    %cases, not all the responses cross the threshold.
    
end

end


