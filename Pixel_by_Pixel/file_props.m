
function [file_locs,num_divisions,division_Pixels,numParts] = file_props(targetDir,search_exp,varargin)


%% This is function argument handling
if nargin>2

    for field_num=1:(nargin-2)
        
        currentFN = char(varargin(field_num));
        %if it finds a field with the name celltype
        if length(regexp(currentFN,'available_memory='))>0
            
            %find out what sample rate
            sr = regexp(currentFN,'available_memory=','split');
            %and assign it to a variable
            available_memory = str2double(sr(2));
        end
        %Same syntax as above   
         if length(regexp(currentFN,'frameSize='))>0
            slp_time = regexp(currentFN,'frameSize=','split');
            slp_time(2)
            frameSize = str2num(slp_time{2});
            
         end

            
     end
    
end

%% This determines file directories

dirFiles = dir(targetDir);

for file_idx=1:length(dirFiles)
    
    isRelevant(file_idx) = ~isempty(regexp(dirFiles(file_idx).name,search_exp));
end

relevant_nums = find(isRelevant);
numParts = length(relevant_nums);

ii = 1;
totalSize = 0;
for idx = relevant_nums
    
    file_locs{ii} = dirFiles(idx).name;
    totalSize = totalSize +dirFiles(idx).bytes;

    ii = ii+1;
end


%% Here determine the number fo divisions and their size in case it is a raw file

if (exist('available_memory','var') && exist('frameSize','var') )
    
    maxGB_singleFile = available_memory/2;
    max_mem_bytes = maxGB_singleFile*10^9;
    
    %This works out the nearest power of two that fits, so can divide the
    %process into this number of parts
    num_divisions = 2^nextpow2(totalSize./max_mem_bytes);

    division_Pixels = frameSize(2)/num_divisions;
else
    num_divisions = NaN;
    division_Pixels = NaN;
end

