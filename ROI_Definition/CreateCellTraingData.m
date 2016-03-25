
load('/Users/Yves/Desktop/Bouton_Grabinfo_files/Area02/(20130514_15_01_53)-_Area1_73_noise_GRABinfo.mat')


Im = Visible_Boutons(GRABinfo.AllMeanStack,0);

roifig=figure('Position',get(0,'ScreenSize')); imshow(Im);

ii = 0;

radius = 8;

while ishandle(roifig)


% Select your circle
[x, y] = ginput(1);

ii = ii + 1;


if x ~= 0
    
Xx(ii) = round(x);
Yy(ii) = round(y);
viscircles([x, y], (radius-4),'EdgeColor','b','LineWidth',0.2);
end

if ii - length(Xx) > 0    
close gcf
end

x = [];
y = [];

end


  ii = length(Xx);
  
  [ymax xmax] = size(GRABinfo.AllMeanStack);
  
% Check that everything is within bounds
toosmall_x = find(Xx<=radius);
Xx(toosmall_x) = radius+2;

toosmall_y = find(Yy<=radius);
Yy(toosmall_y) = radius+2;

toobig_x = find(Xx>=(xmax - radius));
Xx(toobig_x) = (xmax - radius) - 2;

toobig_y = find(Yy>=(xmax - radius));
Yy(toobig_y) = (ymax - radius) - 2;
    

  A = [];

 for j = 1:ii
 A(:,:,j) = Im(Yy(j)-radius:Yy(j)+radius,Xx(j)-radius:Xx(j)+radius);
 end
 
 
 %% How to create a test matrix.
 
 %This script will subdivide the image into 'sumimages' with size 17x17. It
 %will oversample the image, so that after creating one image, it will move
 %stepSize pixels along and then after finishing that column will move down
 %4 pixels and do the same thing for the next row.
 
 load('/Users/Yves/Desktop/Bouton_Grabinfo_files/Area05/(20130514_18_28_52)-_Area5_933x-425y-38z_tones3_GRABinfo.mat')

Im = Visible_Boutons(GRABinfo.AllMeanStack,0);
ImSize = size(Im);


radius = 8;

stepSize = 4;
BoutonRow = []; 
 
 for jj = 1:floor(ImSize(1) - stepSize*radius)/stepSize

for ii = 1:floor(ImSize(1) - stepSize*radius)/stepSize
    
    
    BoutonRow(:,:,ii + (jj-1)*120) = Im(jj*stepSize:(jj*stepSize) +2*radius,ii*stepSize:(ii*stepSize)+2*radius);
    BoutonIm = Visible_Boutons(BoutonRow(:,:,ii),0);
    %imshow(BoutonIm)
    

    
    %a = input('contains_bouton');
    %BoutonPresent(ii) = a;
    
end


 end

 Training.Images = BoutonRow;
 %% 
 
 %This just suppresses small radius warning messages
 w.id = 'images:imfindcircles:warnForSmallRadius';
 warning('off',w.id)
rmpath('folderthatisnotonpath');
 

 %This will 
 a = find(Group ==1);
 
 plar = 1;
 
 for i = 1:length(a)
     
Images(:,:,i) = Visible_Boutons(Training.Images(:,:,a(i)),0);
     
%pp(:,:,i) = im2bw(Training.Images(:,:,a(i)));
%figure, imshow(pp(:,:,i))

%figure, imshow(imgradient(Images(:,:,i),'Prewitt'));

%figure, imshow(edge(Images(:,:,i),'canny'))


%par(:,:,i) = watershed(Images(:,:,i),8);
%figure, imshow(Training.Images(:,:,a(i)))

 [centersBright, radiiBright] = imfindcircles(Images(:,:,i),[1 5], 'Sensitivity',0.80,'Objectpolarity', 'bright','Edgethreshold',0.05,'Method','Twostage');
 
 if numel(centersBright) ~= 0
     centers(plar,:) = centersBright(1,:);
     radii(plar) = radiiBright(1);
     viscircles(centersBright(1,:), (radiiBright(1)+2),'EdgeColor','b','LineWidth',0.5);
     
     plar = plar+1;
 end

 end
 
 
 

%% Select the annulus of supposed bright pixels.

[rr, cc] = meshgrid(1:size(Images,1));


for j = 1:length(centers)
mask1(:,:,j) = sqrt((rr-centers(j,1)).^2+(cc-centers(j,2)).^2)<=radii(j);
%mask2(:,:,j) = sqrt((rr-centers(j,1)).^2+(cc-centers(j,2)).^2)>=radii(j);
end
mask = mask1;%.*mask2;
 
 

 
 
%% Now need to create the layover mask

 Mask = zeros(size(GRABinfo.AllMeanStack));

%Works, note that 120 is used because floor(ImSize(1) -
%stepSize*radius)/stepSize =120 , where ImSize refers to the original image
%not subimages
 
stepSize = 4;
radius = 6;

 for i = 1:size(mask,3)
     
     Mask((1+floor(a(i)/120))*stepSize:(1+floor(a(i)/120))*stepSize + 2*radius,rem(a(i),120)*stepSize:rem(a(i),120)*stepSize+2*radius) = Mask((1+floor(a(i)/120))*stepSize:(1+floor(a(i)/120))*stepSize + 2*radius,rem(a(i),120)*stepSize:rem(a(i),120)*stepSize+2*radius)+ mask(:,:,i);
     
     
 end
 
 finalMask = im2bw(Mask,0.8);
 
 
 %% Overlay the two
 
 
 
 % E is the thing you want to overlay on the original.
 E = Visible_Boutons(GRABinfo.AllMeanStack,0); 
 imshow(E, 'InitialMag', 'fit')
 
 I = finalMask*0.5; 
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

