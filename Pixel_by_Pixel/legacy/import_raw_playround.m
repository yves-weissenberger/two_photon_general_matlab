

fname = '/Users/Yves/Desktop/Area01/(20140121_13_16_08)-_MGB20140121_Area01_noise_reg.tif';
InfoImage=imfinfo(fname);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(fname, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   FinalImage(:,:,i)=TifLink.read();
end
TifLink.close();


%% Method1 for the whole frame extraction
tic 
%fname2 = '/Users/Yves/Desktop/20140121_regRaw/Area01/(20140121_13_16_08)-_MGB20140121_Area01_noise_reg.raw';
fname3 = '/Users/Yves/Desktop/20140121_regRaw/Area01/(20140121_13_20_45)-_MGB20140121_Area01_tones1_reg.raw';


fname=fname2;

if ~isempty(strfind(fname,'noise'))
    numFrames=1800;
else
    numFrames=13500;
end

fid=fopen(fname,'r','l');
%img_line = zeros(1,numFrames*512.^2,'uint16');
img_line = fread(fid, numFrames*512.^2, 'uint16=>uint16').';

%time_series = reshape(img_line,512,512,numFrames);
%clear img_line
toc



%% method2 for the whole frame extraction


%fname2 = '/Users/Yves/Desktop/20140121_regRaw/Area01/(20140121_13_16_08)-_MGB20140121_Area01_noise_reg.raw';
fname3 = '/Users/Yves/Desktop/20140121_regRaw/Area02/(20140121_13_52_14)-_MGB20140121_Area02_noise_reg.raw';


fid=fopen(fname3,'r','l');
tic
if ~isempty(strfind(fname3,'noise'))
    numFrames=1800;
else
    numFrames=13500;
end

time_series2 = zeros(512,512,numFrames,'uint16');
for frame_num=1:numFrames
    
    if frame_num<10
        readPos = ftell(fid)/2;
        sprintf('now reading frame %d',readPos/512^2)
    end
    time_series2(:,:,frame_num) = reshape(fread(fid,512^2,'uint16'),512,512);
    %time_series2(:,frame_num) = fread(fid,512^2,'uint16');

end

toc

%% Understanding


%fname2 = '/Users/Yves/Desktop/20140121_regRaw/Area01/(20140121_13_16_08)-_MGB20140121_Area01_noise_reg.raw';
fclose(fid);
fid=fopen(fname2,'r','l');




%Ok, so unit16 data is stored as 2 bytes per uint16 number

%To obtain the size of a file, go
fileRef = dir(fname2);
file_size = fileRef.bytes;

%Because of the statement above, numel(time_series)*2==file_size

%% Final answer this extracts the correct pixels


fclose(fid);
fid=fopen(fname,'r','l');

one_pixel_trace = zeros(1,numFrames,'uint16');

pixel_offset = 0;
offset = pixel_offset*2;
tic
for frame_num=1:numFrames
    
    
    skip_uint16s = 512^2;
    skip_bytes = skip_uint16s*2;
    fseek(fid,offset+skip_bytes*(frame_num-1),'bof');
    if frame_num<=5
        ftell(fid)
    end
    one_pixel_trace(frame_num) = fread(fid,1,'uint16=>uint16');
end
toc





    %fseek(fid,skip_bytes*(frame_num-1),'bof');

%% Or in one line extract the single trial_trace. This is faster for single line

fclose(fid);
fid=fopen(fname2,'r','l');

skip_uint16s = 512^2;
skip_bytes = skip_uint16s*2;

tic
fseek(fid,0,'bof');
%I have no idea why you have to subtract two...
one_pixel_trace2 = fread(fid,1800,'uint16=>uint16',skip_bytes-2,'l');

toc

%% Now go for parallel



%fclose(fid);
fid=fopen(fname3,'r','l');

segment_trace = zeros(512,64,numFrames,'uint16');

pixel_offset = 512;
offset = pixel_offset*2;
tic
for frame_num=1:numFrames
    
    
    skip_uint16s = 512^2;
    skip_bytes = skip_uint16s*2;
    fseek(fid,offset+skip_bytes*(frame_num-1),'bof');
    if frame_num<=5
        ftell(fid)
    end
    segment_trace(:,:,frame_num) = fread(fid,[512,64],'uint16=>uint16');
end
toc
