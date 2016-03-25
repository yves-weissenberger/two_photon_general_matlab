% Semi Automated ROI

A = GRABinfo.AllMeanStack;
normalising_factor = max(max(A));
A = A./normalising_factor;



A_I = A + 30;
A_I = A_I(:,:,1); 
% 
% %%%%%%% This is only a test
% se = strel('disk',1);
% I2 = imerode(A_I,se); figure, imagesc(I2); colormap('gray')
% A_I = I2;
% %%%%%%%

% morphological opening is basically an erosion followed by a dilation
background = imopen(A_I,strel('disk',15));

% Display the Background Approximation as a Surface
%set(gca,'ydir','reverse');
A_I2 = A_I - background;
A_I3 = imadjust(A_I2);
figure; imshow(A_I3);

xmax = 512;
ymax = 512; %size(GRABinfo.AllMeanStack);
radius = 13;

% Select your circle
point_coordinates = round(ginput(1));
y = point_coordinates(1);
x = point_coordinates(2);

B = A_I(x-radius:x+radius,y-radius:y+radius);
% Corresponding matrix element
%figure; imshow(B);
%point = y*xmax + x;

I3 = B;


I = B + 30;
I = I(:,:,1); 

% morphological opening is basically an erosion followed by a dilation
background = imopen(I,strel('disk',15));

% Display the Background Approximation as a Surface
%set(gca,'ydir','reverse');
I2 = I - background;
I3 = imadjust(I2);
figure; imshow(I3);

% %% K-means
% 
% %a = 5;
% I2 = im2double(I3)*a;
% I3 = reshape(I2,[],1);
% idx = kmeans(I3,2,'Replicates',5);
% idx2 = reshape(idx, size(I2));
% I4 = I2.*(idx2==1);
% I5 = I2.*(idx2==2);
% figure
% imshow(I4);
% % 
% %% Sobel Filtering
% I = I3;
% 
% hy = fspecial('sobel');
% hx = hy';  
% Iy = imfilter(double(I), hy, 'replicate');
% Ix = imfilter(double(I), hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);
% figure, imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
% % test = im2bw(gradmag);
% % figure, imshow(test)
% % bwlabel(test,8)
% 
% %% Other Stuff
% 
% I = I3;
% se = strel('disk', 2);
% Io = imopen(I, se);
% figure, imshow(Io), title('Opening (Io)')
% 
% 
% Ie = imerode(I, se);
% Iobr = imreconstruct(Ie, I);
% figure, imshow(Iobr), title('Opening-by-reconstruction (Iobr)')
% 


%% Identify Circles

%subplot(1,2,1), imshow(I_BW)
[centersDark, radiiDark] = imfindcircles(I3,[5 7], 'Sensitivity',0.90,'Objectpolarity', 'dark','Edgethreshold',0.05,'Method','Twostage');
viscircles(centersDark(1,:), (radiiDark(1)),'EdgeColor','b','LineWidth',0.5);

% %%
% 
% level = graythresh(I3);
% bw = im2bw(I3,level);
% bw = bwareaopen(bw, 50);
% figure, imagesc(bw)

%% Create the corresponding mask
centre_row = centersDark(1,1);
centre_column = centersDark(1,2);

radius = radiiDark(1);

% Select the annulus of supposed bright pixels.
[rr, cc] = meshgrid(1:length(I3));
mask1 = sqrt((rr-centre_row).^2+(cc-centre_column).^2)<=radius + 4;
mask2 = sqrt((rr-centre_row).^2+(cc-centre_column).^2)>=radius;
mask = mask1.*mask2;

figure, imshow(mask)

%% Layover mask and find mean brightness of annulus

after_mask = I3.*mask;
figure, imshow(after_mask);

% Find nonzero elements after the mask
% nonzero_idx = find(after_mask);
% mean_bright = max(after_mask(nonzero_idx));

%  per = imfill(after_mask);
%  figure, imshow(per)

% %% Fill in the black circle
% 
% black_circle = sqrt((rr-centre_row).^2+(cc-centre_column).^2)<=radius + 1;
% 
% black_idx = find(black_circle);
% 
% After_fill = I3;
% After_fill(black_idx) = mean_bright;
% 
% figure, imshow(After_fill);
% 
% %% Integrate with the full picture
% 
% Final_image = I3;
% 
% light_idx = find(per);
% Final_image(light_idx) = max(;
% figure, imshow(Final_image);
% 
% %%
% [centersBright, radiiBright] = imfindcircles(Final_image,[10 20], 'Sensitivity',0.90,'Objectpolarity', 'bright','Edgethreshold',0.2,'Method','Twostage');
% viscircles(centersBright, (radiiBright),'EdgeColor','b','LineWidth',0.5);
% 

%%
close all

a = find(after_mask);
b = mean2(after_mask(a));
c = im2bw(I3, b);
d = after_mask | c;
e = imopen(d, strel('disk',1));
figure, imshow(I3);

[L, num] = bwlabel(e);

obj = zeros(1,num);
for i = 1:num
    
    obj(i) = length(find(L==i));
end

[size_obj, index] =  max(obj(:));

cell = (L==index);
figure, imshow(cell)
