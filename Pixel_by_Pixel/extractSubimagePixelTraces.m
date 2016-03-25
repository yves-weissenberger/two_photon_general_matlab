
function [traces] = extractSubimagePixelTraces(file_locs,frameSize,numFrames,num_divisions,division_nr,targetDir)


num_files = length(file_locs);

traces = zeros(frameSize(1),frameSize(2)/num_divisions,num_files*numFrames,'uint16');
for file_idx = 1:num_files
    
    %get the file location
    fname = fullfile(targetDir,file_locs{file_idx});


    %open the file as read only, read in little endian order

    %initialise the data
    segment_trace = zeros(frameSize(1),frameSize(2)/num_divisions,numFrames,'uint16');

    pixel_offset = (frameSize(1)*frameSize(2)/num_divisions)*(division_nr-1);
    offset = pixel_offset*2;
    
    
    
    if exist('fid','var');
        fid=fclose(fid);
    end
    
    fid=fopen(fname,'r','l');
    
    for frame_num=1:numFrames


        skip_uint16s = 512^2;
        skip_bytes = skip_uint16s*2;
        fseek(fid,offset+skip_bytes*(frame_num-1),'bof');
%         if frame_num<4
%             ftell(fid)
%         end
        segment_trace(:,:,frame_num) = fread(fid,[frameSize(1),frameSize(2)/num_divisions],'uint16=>uint16');
    end
    
    traces(:,:,1+numFrames*(file_idx-1):file_idx*numFrames) = segment_trace;


end