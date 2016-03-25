% Semi Automated ROI
%function = [GRABinfo.bwcells, GRABinfo.
clear Xx
clear Yy

%load('/Users/Yves/Documents/MATLAB/20131017/Area08/(20131017_21_00_48)-_Cortex_20131017_Area10_noise_GRABinfo.mat')


A = GRABinfo.AllMeanStack;
normalising_factor = max(max(A));
A = A./normalising_factor;



A_I = A + 30;
A_I = A_I(:,:,1); 


% morphological opening is basically an erosion followed by a dilation
background_full = imopen(A_I,strel('disk',15));

% Display the Background Approximation as a Surface
%set(gca,'ydir','reverse');
A_I2 = A_I - background_full;
A_I3 = imadjust(A_I2);


% %% Remove boutons
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
radius = 17;


roifig=figure('Position',get(0,'ScreenSize')); imshow(A_I3);

ii = 0;

while ishandle(roifig)


% Select your circle
[x, y] = ginput(1);

ii = ii + 1;


if x ~= 0
    
Xx(ii) = round(x);
Yy(ii) = round(y);
viscircles([x, y], radius,'EdgeColor','b','LineWidth',0.5);
end

if ii - length(Xx) > 0    
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

    ii = length(Xx);
    
   

 for j = 1:ii
 B(:,:,j) = A_I(Yy(j)-radius:Yy(j)+radius,Xx(j)-radius:Xx(j)+radius);
 end
% % Corresponding matrix element
% %figure; imshow(B); 
I3 = B;

%%
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



%% Identify Circles
%subplot(1,2,1), imshow(I_BW)
for j = 1:ii
[CentersDark, RadiiDark] = imfindcircles(I3(:,:,j),[5 7], 'Sensitivity',1,'Objectpolarity', 'dark','Edgethreshold',0.001,'Method','Twostage');

if RadiiDark ~= 0
centersDark(j,1) = CentersDark(1,1);
centersDark(j,2) = CentersDark(1,2);
radiiDark(j) = RadiiDark(1);
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

for j = 1:ii
figure, subplot(1,2,2), imshow(cells(:,:,j)); subplot(1,2,1), imshow(I3(:,:,j));
end

%% Create the complete mask

Size_Mask = size(A_I3);
Final_Mask = zeros(Size_Mask);

for j = 1:ii
    
  Final_Mask(Yy(j)-radius:Yy(j)+radius,Xx(j)-radius:Xx(j)+radius)  =  Final_Mask(Yy(j)-radius:Yy(j)+radius,Xx(j)-radius:Xx(j)+radius) + cells(:,:,j);
    
end

    figure, imshow(Final_Mask)
    
    
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

 
 %% Breakcells have a look. bwcells is used by the ROI selection GUI that we already have. ROImap is the file labelled with numbered ROIs
    