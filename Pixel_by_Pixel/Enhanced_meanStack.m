function [Stack] = Enhanced_meanStack(fname)


%tic
if ~isempty(strfind(fname,'noise'))
    numFrames=1800;
else
    numFrames=13500;
end


time_series2 = zeros(512,512,numFrames,'uint16');

fid=fopen(fname,'r','l');
for frame_num=1:numFrames
    
%     if frame_num<10
%         readPos = ftell(fid)/2;
%         sprintf('now reading frame %d',readPos/512^2)
%     end
    time_series2(:,:,frame_num) = reshape(fread(fid,512^2,'uint16'),512,512);
    %time_series2(:,frame_num) = fread(fid,512^2,'uint16');

end

fclose('all');
%toc


Stack = Visible_Boutons(mean(time_series2,3),0);