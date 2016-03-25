% Semi Automated ROI detection Script
%Yves Weissenberger Jan 2014; yves.weissenberger@magd.ox.ac.uk

%% Open relevent files
function  [ROImaps,Xpil,Ypil,num, ii] = Semi_automated_ROIfnctn(targetDir, Date, Area, subplots)

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


%% Preprocessing.
normalising_factor = max(max(A));
A = A./normalising_factor;



%A_I = A + 30;
A_I = A(:,:,1); 


% morphological opening is basically an erosion followed by a dilation
background_full = imopen(A_I,strel('disk',15));

% Display the Background Approximation as a Surface
%set(gca,'ydir','reverse');
A_I2 = A_I - background_full;
A_I3 = imadjust(A_I2);



% A_I3 = Visible_Boutons(A);
% A_I = A_I3;


% % %% Remove boutons??

% 
% R = rangefilt(A_I3,ones(3));
% imshow(R)
% 
% [centersBright, radiiBright] = imfindcircles(R,[1 3], 'Sensitivity',0.80,'Objectpolarity', 'bright','Edgethreshold',0.2,'Method','Twostage');
% viscircles(centersBright, (radiiBright),'EdgeColor','b','LineWidth',0.5);
% 
% 

%%
xmax = 512;
ymax = 512; %size(GRABinfo.AllMeanStack);
radius = 13;
diameter = 2*radius + 1;


roifig=figure('Position',get(0,'ScreenSize')); imshow(A_I3);

lengthEveryThing = 0;

while ishandle(roifig)


% Select your circle
[x, y] = ginput(1);

lengthEveryThing = lengthEveryThing + 1;


if x ~= 0
    
Xx(lengthEveryThing) = round(x);
Yy(lengthEveryThing) = round(y);
viscircles([x, y], radius,'EdgeColor','b','LineWidth',0.5);
end

if lengthEveryThing - length(Xx) > 0    
close gcf
end

x = [];
y = [];

end

% Check that everything is within bounds
toosmall_x = find(Xx<=radius);
Xx(toosmall_x) = radius+1;

toosmall_y = find(Yy<=radius);
Yy(toosmall_y) = radius+1;

toobig_x = find(Xx>=(xmax - radius));
Xx(toobig_x) = (xmax - radius) - 1;

toobig_y = find(Yy>=(xmax - radius));
Yy(toobig_y) = (ymax - radius) - 1;


%%
% y = point_coordinates(1);
% x = point_coordinates(2);

    lengthEveryThing = length(Xx);
    
   

 for jj = 1:lengthEveryThing
 B(:,:,jj) = A_I3(Yy(jj)-radius:Yy(jj)+radius,Xx(jj)-radius:Xx(jj)+radius);
 
 end
% % Corresponding matrix element
% %figure; imshow(B); 
I3 = B;

%%
%I = B + 30;
% I = B;
% 
% % morphological opening is basically an erosion followed by a dilation
% for j = 1:length(Xx)
% background(:,:,j) = imopen(I(:,:,j),strel('disk',15));
% end
% 
% % Display the Background Approximation as a Surface
% %set(gca,'ydir','reverse');
% I2 = I - background;
% 
% for j = 1:ii
% I3(:,:,j) = imadjust(I2(:,:,j));
% end



%% Identify Circles
%subplot(1,2,1), imshow(I_BW)

warning('off','all');
warning;

for jj = 1:lengthEveryThing
[CentersDark, RadiiDark] = imfindcircles(I3(:,:,jj),[2 8], 'Sensitivity',1,'Objectpolarity', 'dark','Edgethreshold',0.001,'Method','Twostage');

if RadiiDark ~= 0
centersDark(jj,1) = CentersDark(1,1);
centersDark(jj,2) = CentersDark(1,2);
radiiDark(jj) = RadiiDark(1);
end



% Put some code to detect the bright radii here as well, if no dark ones
% are detected; as an elseif
end

warning('on','all');
warning;


%% Create the corresponding mask
[m, n, o] = size(I3);

% Select the annulus of supposed bright pixels.
[rr, cc] = meshgrid(1:m);

CircCentre_small = radius;

CentreCircle_small = sqrt((rr-CircCentre_small).^2+(cc-CircCentre_small).^2)<=1;


for jj = 1:lengthEveryThing
    
    ROIcentre = (I3(:,:,jj).*CentreCircle_small);
    ROInonzero = find(ROIcentre);
    
    if mean(ROIcentre(ROInonzero)) >=0.8
        
        mask1(:,:,jj) = regiongrowing(I3(:,:,jj),17,17,0.09);
        mask2 = ones(size(mask1));
        
        
    else
        mask1(:,:,jj) = sqrt((rr-centersDark(jj,1)).^2+(cc-centersDark(jj,2)).^2)<=radiiDark(jj) + 4;
        mask2(:,:,jj) = sqrt((rr-centersDark(jj,1)).^2+(cc-centersDark(jj,2)).^2)>=radiiDark(jj);
        
    end
    
end
mask = mask1.*mask2;


%% Layover mask

after_mask = I3.*mask;


%%
for jj = 1:lengthEveryThing;

%Find the mean pixel intensity in the mask
z = after_mask(:,:,jj);
a = find(z);
b = mean2(z(a));

%find the maximum pixel intensity in the mask area
bb = max(max(z(a)));

%threshold
c = im2bw(I3(:,:,jj), (1.7*b + 0.3*bb)/2);

if mean2(I3(:,:,jj).*CentreCircle_small) <0.8

    d = z | c;

else
    
    d = c;

end


e = imopen(d, strel('disk',1));

%%

%POTENTIAL PROBLEM, might be a reason why cells are sometimes not
%selected!!
[L, num] = bwlabel(e,4);


%%


% This is the script that removes small side objects.
for ii = 1:lengthEveryThing
    
    obj(ii) = length(find(L==ii));
end

[size_obj, index] =  max(obj(:));

cell = (L==index);

cells(:,:,jj) = cell;

end

if subplots == 1
    
for jj = 1:lengthEveryThing
figure, subplot(1,2,2), imshow(cells(:,:,jj)); subplot(1,2,1), imshow(I3(:,:,jj));
end

end

%% Create the complete mask

Size_Mask = size(A_I3);
Final_Mask = zeros(Size_Mask);

for jj = 1:lengthEveryThing
    
  Final_Mask(Yy(jj)-radius:Yy(jj)+radius,Xx(jj)-radius:Xx(jj)+radius)  =  Final_Mask(Yy(jj)-radius:Yy(jj)+radius,Xx(jj)-radius:Xx(jj)+radius) + jj.*cells(:,:,jj) + 0.0005;
  
  Final_Mask(find(Final_Mask==0.0005)) = 0;

    
end


%% Create the Neuropil Mask

%as output basically want an extra GRABinfo thing that has the [y x]
%coordinates of the ROImasks for each cell

Ypil = zeros((2*diameter)^2,lengthEveryThing);
Xpil = zeros((2*diameter)^2,lengthEveryThing);

Xxpil = Xx;
Yypil = Yy;

% Check that everything is within bounds
toosmall_xpil = find(Xx<=2*radius);
Xxpil(toosmall_xpil) = 2*radius+1;

toosmall_ypil = find(Yy<=2*radius);
Yypil(toosmall_ypil) = 2*radius+1;

toobig_xpil = find(Xx>=(xmax - 2*radius));
Xxpil(toobig_xpil) = (xmax - 2*radius) - 1;

toobig_ypil = find(Yy>=(xmax - 2*radius));
Yypil(toobig_ypil) = (ymax - 2*radius) - 1;


%This could go wrong, if x and y are zero for a given lookup.
for jj = 1:lengthEveryThing
    
    
  empty_Mask = zeros(size(Final_Mask));
  
  cellpill = imdilate(Final_Mask(Yypil(jj)-2*radius:Yypil(jj)+2*radius,Xxpil(jj)-2*radius:Xxpil(jj)+2*radius),strel('disk',4));
  cellpill = imfill(cellpill,'holes');
  neuropil(:,:,jj) = ones(size(cellpill)) - cellpill;
  
  %This allows easy lookup of the x and y positions of the neuropil mask
  empty_Mask(Yypil(jj)-2*radius:Yypil(jj)+2*radius,Xxpil(jj)-2*radius:Xxpil(jj)+2*radius) = neuropil(:,:,jj);
  
  [y,x] = find(empty_Mask==1);

  Ypil(1:length(y),jj) = y;
  Xpil(1:length(x),jj) = x;
  
end


%% Code to work backwards to check that the Neuropil mask thing works
% yp = Ypil(:,2);
% xp = Xpil(:,2);
% b = find(xp);
% 
% b = find(xp);
% c = find(yp);
% g = sub2ind(size(empty_Mask),yp(c),xp(b));
% aa = zeros(512);
% aa(g) = 1;
% imshow(aa)

%% Crate the real final mask, removing overlapping regions

%This finds pixels with overlapping labels
overlapPix = Final_Mask>0 & rem(Final_Mask - 0.0005,1)==0;

%This then creates the inverse of the overlapPix matrix
overlapPix2 = overlapPix==0;


%This then removes those pixels from the ROI mask
ROImaps = Final_Mask.*overlapPix;


%Finally remove the 0.5 bit from the label
ROImaps(find(ROImaps)) = ROImaps(find(ROImaps)) -0.0005;


    %% Overlay the Mask onto the real image
    
% Two colour image

ROIcolour = label2rgb(ROImaps);


OrgImg = imshow(1.8*A_I3);

hold on
Colourshow = imshow(ROIcolour);
hold off

set(Colourshow, 'AlphaData', OrgImg)


% 
%  % E is the thing you want to overlay on the original.
%  E = A_I3; 
%  imshow(E, 'InitialMag', 'fit')
%  
%  I = im2bw(Final_Mask,0.8).*0.2; 
%  imshow(I, 'InitialMag', 'fit')
%  
%  imshow(E, 'InitialMag', 'fit')
%  % Make a truecolor all-green image.
%  green = cat(3, zeros(size(E)), ones(size(E)), zeros(size(E)));
%  hold on 
%  h = imshow(green); 
%  hold off
%  
%   % Use our influence map as the 
%  % AlphaData for the solid green image.
%  set(h, 'AlphaData', I)
% 
% 


%% Save Final Results    

acceptstr=input('Accept ROIs? y/n  :','s');


    if strcmpi(acceptstr,'y');
        GRABinfo.bwcells = Final_Mask;
    
        GRABinfo.ROImaps = ROImaps;

        for jj=1:size(Area_files,1);

            GRABname=Area_files(jj).name;
            load(fullfile(targetDir,Date,Area,GRABname));
            GRABinfo.bwcells=Final_Mask;
            GRABinfo.ROImaps = ROImaps;
            GRABinfo.Xpil = Xpil;
            GRABinfo.Ypil = Ypil;
            save(fullfile(targetDir,Date,Area,GRABname),'GRABinfo');

        end;
        
    end
    

end
 
 %% Breakcells have a look. bwcells is used by the ROI selection GUI that we already have. ROImap is the file labelled with numbered ROIs
    