% Semi Automated ROI detection Script
%Yves Weissenberger Jan 2014; yves.weissenberger@magd.ox.ac.uk

%% Open relevent files
function  Semi_automated_ROI(targetDir, Date, Area, subplots)

 clear Xx
 clear Yy
 clear B


%dir returns all files in a given directory. f = fullfile('myfolder','mysubfolder','myfile.m')
%returns a string combining all the inputs
Area_files = dir([fullfile(targetDir,Date,Area) '/*.mat']);

%Load image to be processed
GRABname=Area_files(1).name;
load(fullfile(targetDir,Date,Area,GRABname));

%Deafults accepting output of the function to no
acceptstr='n';

A = GRABinfo.AllMeanStack;


%% Preprocessing: normalises A by its greatest value to have a distribution of pixel intensities between 0 and 1.
normalising_factor = max(max(A));
A = A./normalising_factor;


%Just makes the image brighter, shouldn't make a computational difference
A_I = A + 30;

% In case A_I has more than 2 dimension reduces the dimensionality to 2 to
% avoid errors in output from other functions
A_I = A_I(:,:,1); 


% morphological opening (imopen) is a matlab function that performs and erosion
%of the image followed by a dilation using the structuring element (strel) specified as
%the second input to the imopen function. This erodes away most of the
%foreground objects whose radius is considerably less than 15.
background_full = imopen(A_I,strel('disk',15));

%Subtract background from original image to enhance the foreground edges
%and map the intensity values between 0 and 1.
A_I2 = A_I - background_full;
A_I3 = imadjust(A_I2);



%%


image_dimensions = size(GRABinfo.AllMeanStack);

%Note that using imshow, rows and columns of pixels define columns and rows
%of the image, respecitvely.
xmax = image_dimensions(2);
ymax = image_dimensions(1);

% Specifies the radius of the area around the points you click that will be
% used to create the subimages.
radius = 17;

%Returns a handle for the figure that is used to select the regions of
%interest
roifig=figure('Position',get(0,'ScreenSize')); imshow(A_I3);

ii = 0;

%while this figure is open
while ishandle(roifig)


% Select your circle
[x, y] = ginput(1);

ii = ii + 1;


if x ~= 0
    
Xx(ii) = round(x);
Yy(ii) = round(y);
viscircles([x, y], radius,'EdgeColor','b','LineWidth',0.5);
end

% If you press enter, ginput will output [], so Xx will not be extended
% meaning ii will be greater than length(Xx) and the figure will be closed
% ending the while loop
if ii - length(Xx) > 0    
close gcf
end

x = [];
y = [];

end

% Check that the selected regions of interest are entirely within the
% bounds of the image, otherwise they are moved.
toosmall_x = find(Xx<=radius);
Xx(toosmall_x) = radius+1;

toosmall_y = find(Yy<=radius);
Yy(toosmall_y) = radius+1;

toobig_x = find(Xx>=(xmax - radius));
Xx(toobig_x) = (xmax - radius) - 1;

toobig_y = find(Yy>=(xmax - radius));
Yy(toobig_y) = (ymax - radius) - 1;


%% Selects your preliminary ROIs around the points that you clicked
ii = length(Xx);
    
   

 for j = 1:ii
 B(:,:,j) = A_I(Yy(j)-radius:Yy(j)+radius,Xx(j)-radius:Xx(j)+radius);
 end
% % Corresponding matrix element
% %figure; imshow(B); 
I3 = B;

%% The same preprocessing that was done on the original of the entire image above is done on the preliminary ROIs that were clicked
I = B + 30;

% morphological opening is basically an erosion followed by a dilation
for j = 1:length(Xx)
background(:,:,j) = imopen(I(:,:,j),strel('disk',15));
end

% Display the Background Approximation as a Surface
%set(gca,'ydir','reverse');
I2 = I - background;

for j = 1:ii
I3(:,:,j) = imadjust(I2(:,:,j));
end



%% Identify Central Dark Circles in the image
%subplot(1,2,1), imshow(I_BW)


for j = 1:ii
[CentersDark, RadiiDark, Metric] = imfindcircles(I3(:,:,j),[5 7], 'Sensitivity',1,'Objectpolarity', 'dark','Edgethreshold',0.001,'Method','Twostage');

if RadiiDark ~= 0
centersDark(j,1) = CentersDark(1,1);
centersDark(j,2) = CentersDark(1,2);
radiiDark(j) = RadiiDark(1);
Circle_strengths(j) = Metric;
end


% Put some code to detect the bright radii here as well, if no dark ones
% are detected; as an elseif
end

%% Create the corresponding mask
[m, n, o] = size(I3);

% Select the annulus of supposed bright pixels.
[rr, cc] = meshgrid(1:m);

for j = 1:ii
mask1(:,:,j) = sqrt((rr-centersDark(j,1)).^2+(cc-centersDark(j,2)).^2)<=radiiDark(j) + 4;
mask2(:,:,j) = sqrt((rr-centersDark(j,1)).^2+(cc-centersDark(j,2)).^2)>=radiiDark(j);
end
mask = mask1.*mask2;


%% Layover mask

after_mask = I3.*mask;


%%
for j = 1:ii;

z = after_mask(:,:,j);
a = find(z);
b = mean2(z(a));
bb = max(max(z(a)));
c = im2bw(I3(:,:,j), (b + bb)/2);
d = z | c;


e = imopen(d, strel('disk',1));



[L, num] = bwlabel(e,4);

for i = 1:num
    
    obj(i) = length(find(L==i));
end

[size_obj, index] =  max(obj(:));

cell = (L==index);

cells(:,:,j) = cell;

end

if subplots == 1
    
for j = 1:ii
figure, subplot(1,2,2), imshow(cells(:,:,j)); subplot(1,2,1), imshow(I3(:,:,j));
end

end

%% Create the complete mask

Size_Mask = size(A_I3);
Final_Mask = zeros(Size_Mask);

for j = 1:ii
    
  Final_Mask(Yy(j)-radius:Yy(j)+radius,Xx(j)-radius:Xx(j)+radius)  =  Final_Mask(Yy(j)-radius:Yy(j)+radius,Xx(j)-radius:Xx(j)+radius) + cells(:,:,j);
    
end


%% Display the results the of the process
    %% Overlay the Mask onto the real image
    
% Two colour image


 % E is the thing you want to overlay on the original.
 E = A_I3; 
 imshow(E, 'InitialMag', 'fit')
 
 I = Final_Mask.*0.1; 
 imshow(I, 'InitialMag', 'fit')
 
 imshow(E, 'InitialMag', 'fit')
 % Make a truecolor all-green image.
 green = cat(3, zeros(size(E)), ones(size(E)), zeros(size(E)));
 hold on 
 h = imshow(green); 
 hold off
 
  % Use our influence map as the 
 % AlphaData for the solid green image.
 set(h, 'AlphaData', I)




%% Save Final Results    

acceptstr=input('Accept ROIs? y/n  :','s');


    if strcmpi(acceptstr,'y');
        GRABinfo.bwcells = Final_Mask;
    
        GRABinfo.ROImaps = bwlabel(Final_Mask);

        for j=1:size(Area_files,1);

            GRABname=Area_files(j).name;
            load(fullfile(targetDir,Date,Area,GRABname));
            GRABinfo.bwcells=Final_Mask;
            GRABinfo.ROImaps = bwlabel(Final_Mask);
            save(fullfile(targetDir,Date,Area,GRABname),'GRABinfo');

        end;
        
    end
    

end
 
 %% Breakcells have a look. bwcells is used by the ROI selection GUI that we already have. ROImap is the file labelled with numbered ROIs
    