function [img_line,time_series] = import_raw_as_mat(file_loc)


%% Direct Method
tic 
fid=fopen(file_loc);
img_line = zeros(1,13500*512.^2);
img_line = fread(fid, 13500*512.^2, 'uint16=>uint16').';

time_series = reshape(img_line,512,512,13500);
toc

end