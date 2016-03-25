%% PixelWise Two-Photon Data Analysis

%Yves Weissenberger 2015

%This is the master file that calls all subroutines which should be
%included in the Pixel_by_Pixel folder. Careful when further processing
%output of the subroutines because number formats are specified to control
%RAM usage and MATLAB often defaults matrices to double...

%% Specify some Parameters

%The directory that the area folders, which in turn contain the .raw files
%are in
targetDir = '/Users/Yves/Desktop/20140121_regRaw/Area04';


%Specify the size of a single frame that is to be loaded from the data.
frameSize = '512,512';

%specify the number of frames per file
numFramesperFile = 13500;

%Specify what the largest amount of data (in GB) that you are willing to 
%commit to RAM is. Write it as a string
available_memory = '4';

%% Load file data and divide loading 

%Here type the regular expression that allows you to select only the
%relevant files in the target directory
search_exp = '.*tones.*.raw';

%This works for videos whose pixel numbers can specified as 2^n where n is
%an integer. NB numParts is the number of divisions the data is split over


[file_locs,num_divisions,division_Pixels,numParts] = file_props(targetDir,search_exp,strcat('available_memory=',available_memory),strcat('frameSize=',frameSize));



%% Extract the data and calculate the preferred frequency

%initialise the data matrix
PixelWise_selectivity_Matrix = zeros(frameSize(1),frameSize(2),3,'double');


% Open a waitbar
wb = waitbar(0,'Please Wait'); 


try
    
    for division_nr=1:num_divisions

        %Extract the data
        [Pixeltraces] = extractSubimagePixelTraces(file_locs,frameSize,numFramesperFile,num_divisions,division_nr,targetDir);
        
        
       

        t_f01 = 25;
        t_f02= 40500;


        tic
        [d1,d2,d3] = size(Pixeltraces);


        %This works
        interm1 = permute(squeeze(reshape(Pixeltraces,d1,d2,t_f02,d3./t_f02)),[1,2,4,3]);
        clear testarray

        interm2 = reshape(interm1,d1,d2,d3./t_f02,t_f02/t_f01,t_f01);
        clear interm1
        interm3 = min(mean(interm2,5),[],4);
        clear interm2

        %this is the old one
        %temp1 = repmat(Pixels_f0,[1,t_f0,1]);

        temp1 = repmat(interm3,[1,t_f02,1]);
        clear interm3
        f0_mtx = reshape(temp1,d1,d2,[]);
        clear temp1
        
        Pixeltraces2 = (Pixeltraces - uint16(f0_mtx))./uint16(f0_mtx);

        clear temp1
        
        
        
        
        
        

        %missing bit is calculating deltaF/F for all the Pixeltraces

        %Calculate the preferred frequency and frequency selectivity
        %indices of the traces and return them combined as final_Img. This
        %is a custom routine for the MGB imaging study. So should be
        %replaced/updated if used for something else.
        [ImgSegment,~,top_val] = PixelWise_Selectivity(Pixeltraces2);

        %This is an RGB array: size(final_Mat) = [frameSize(1),frameSize(2),3]
        PixelWise_selectivity_Matrix(:,1+division_Pixels*(division_nr-1):division_Pixels*division_nr,:) = ImgSegment;

        %This extract the value from each from each frame to give a normalising
        %constant later
        top_val_Store(division_nr) = top_val;

        %just to show progress
        waitbar(division_nr/ num_divisions,wb,sprintf('now completed %.d%% of the processing',round(100*division_nr/num_divisions)))
    end

catch
    %close the waitbar
    close(wb)
    %Print an error message
    sprintf('Something went wrong :(')
end

%close the waitbar
close(wb)



%% Create the MeanStack Image

%specify the file_name
fname = fullfile(targetDir,file_locs{1});

%extract the meanStack of one file from the ...reg.raw data_files
[MeanStack] = Enhanced_meanStack(fname);

%% Shoooww me seee daaataaa

RGB_image = PixelWise_selectivity_Matrix./repmat(max(top_val_Store),512,512,3);

figure('Units','normalized'),

%img1 = imshow(MeanStack); hold on
img2 = image(RGB_image); hold off

%If this is on and below is off, can see overlay
%img2.AlphaData = 1;

%If is this on and above is off, have colorbar
colormap(jet(25));
cbar = colorbar();
cbar.Ticks = 1:6;
cbar.TickLabels = round(logspace(0,log10(80),6));
cbar.Label.String = 'Preferred Frequency (kHz)';
cbar.FontWeight = 'bold';
cbar.FontSize = 11;



%% This script can be used to threshold the RGB image
%The three planes of hsv_image are hue, saturation, and value components

RGBnew = RGB_image;
% for ii = 1:3
% RGBnew(:,:,ii) = medfilt2(RGB_image(:,:,ii),[1,1]);
% end

hsv_image = rgb2hsv(RGBnew);


%To view image with no depence of intensity on the selectivity index, set
%thresh >=1 and interm1(interm1<thresh) = 1;

thresh = 0.1; 

interm1 = imadjust(hsv_image(:,:,3));
interm1(interm1<thresh) = 0;
hsv_image(:,:,3) = MeanStack.*medfilt2(interm1,[1,1]);%interm1;

RGB_image2 = hsv2rgb(hsv_image);

figure()
%img3 = imshow(MeanStack); hold on
img4 = image(RGB_image2); hold off

%img4.AlphaData = 0.35;



