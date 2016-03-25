  
function [ROI_XYZ_aligned] = extract_pos(short_GRABinfo,varargin)
% extract_pos returns the xyz positions of boutons in the short_GRABinfo
% file
%   'short_GRABinfo' should be a structure that contains at least the
%   fields .xyzPosition and .ROIprops.
%   'varargin' can be used to specify what type of structure is imaged by
%   specifying either 'ROItype=bodies' or 'ROItype=boutons'. Plot
%   determines whether a 3-D plot will be created and keywords are yes + no
%   NeuronXYZ returns a three column vector of X,Y and Z position, relative 
%   to the centre of the entire imaging field on a given day (not a single
%   FOV). Each row then contains the complete coordinates of one ROI
%   each row contains the coordinates
%   
%
%   Example extract_pos(GRABinfo,'celltype=bodies','plot=yes')


if nargin>1

    for field_num=1:(nargin-1)
        
        currentFN = char(varargin(field_num));
        %if it finds a field with the name celltype
        if length(regexp(currentFN,'ROItype='))>0
            
            %find out what celltype it is
            ROI_type_cell = regexp(currentFN,'ROItype=','split');
            %and assign it to a variable
            ROI_type = char(ROI_type_cell(2));
        end
        %Same syntax as above   
         if length(regexp(currentFN,'plot='))>0
            to_plot = regexp(currentFN,'plot=','split');
            plot = char(to_plot(2));
            
         end

            
        end
    
    
    
else
    
    ROI_type = 'bodies';
    plot='no';

end
    %%

    
%% Find the positions of the neurons relative to the centre of the window
    
    if strcmp(ROI_type,'bodies')
    ROICentroids = [(short_GRABinfo.ROIprops(:).Centroid)];
    elseif strcmp(ROI_type,'boutons')
    ROICentroids = [(short_GRABinfo.ROIpropsNew(:).Centroid)];
    else
        error('Entered Invalid celltype')
    end
     
    numCells = length(ROICentroids)/2;
        
    %SHIT think there is a mistake in this script. I think it selects
    %the wrong responsive ROIs.
    
    %This reshapes so that one the first column contains the x positions of
    %the neurons cented on the current ROI and the second column contains
    %the y positions
    ROIXYZ = reshape(ROICentroids,2,numCells)';
    
    x_offset = short_GRABinfo.xyzPosition(1);
    y_offset = short_GRABinfo.xyzPosition(1);
    
    %This is now aligned to the centre of the cranial window for a given
    %experiment
    ROI_XYZ_aligned(:,1) = ROIXYZ(:,1)+x_offset;
    ROI_XYZ_aligned(:,2) = ROIXYZ(:,2)+y_offset;

    
   % NeuronXYZ = NeuronXYZ + FOV_offset;
    
    ROI_XYZ_aligned(:,3) = short_GRABinfo.xyzPosition(3);
    
    %if plot
    if strcmp(plot,'yes')
        scatter3(ROI_XYZ_aligned(:,1),ROI_XYZ_aligned(:,2),ROI_XYZ_aligned(:,3))
    end
    
    
end


%% Notes

%Create a function that goes through the function arguments automatically,
%checking for each of them.

%% Below is an alternative way to process arguments in pairs.
%# define defaults at the beginning of the code so that you do not need to
%# scroll way down in case you want to change something or if the help is
%# incomplete
% options = struct('single','no','secondparameter',magic(3));
% 
% %# read the acceptable names
% optionNames = fieldnames(options);
% 
% %# count arguments
% nArgs = length(varargin);
% if round(nArgs/2)~=nArgs/2
%    error('EXAMPLE needs propertyName/propertyValue pairs')
% end
% 
% for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
%    inpName = lower(pair{1}); %# make case insensitive
% 
%    if any(strcmp(inpName,optionNames))
%       %# overwrite options. If you want you can test for the right class here
%       %# Also, if you find out that there is an option you keep getting wrong,
%       %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
%       options.(inpName) = pair{2};
%    else
%       error('%s is not a recognized parameter name',inpName)
%    end
% end