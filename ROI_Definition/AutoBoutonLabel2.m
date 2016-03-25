

%% This Part subdivides the image into small subimages each of which will be labelled with bouton present or bouton absent.

%clear all

%load('/Users/Yves/Desktop/Bouton_Grabinfo_files/Area01/(20130514_15_08_27)-_Area1_73_tones1_GRABinfo.mat')
%load('/Users/Yves/Desktop/Bouton_Grabinfo_files/Area04/(20130514_17_29_49)-_Area4_607x-428y-38z_tones1_GRABinfo.mat')
% load('/Users/Yves/Desktop/2p_MGB_imaging/Boutons_preprocessed/20140116/Area09/(20140116_20_29_38)-_MGB20140116_Area9_tones1_GRABinfo.mat')
load('SVMHansel.mat')
%clearvars -except SVMStruct


Img = sessionMean; %Visible_Boutons(GRABinfo.AllMeanStack,0);
ImSize = size(Img);

tic

radius = 8;

stepSize = 4;

%initialise
BoutonRow = zeros(2*radius+1, 2*radius+1,(floor(ImSize(1) - stepSize*radius)/stepSize)^2);
BoutonIm = zeros(size(BoutonRow));


for jj = 1:floor(ImSize(1) - stepSize*radius)/stepSize
    
    for ii = 1:floor(ImSize(1) - stepSize*radius)/stepSize
        
        
        BoutonRow(:,:,ii + (jj-1)*120) = Img(jj*stepSize:(jj*stepSize) +2*radius,ii*stepSize:(ii*stepSize)+2*radius);
        BoutonIm(:,:,ii + (jj-1)*120) = Visible_Boutons(BoutonRow(:,:,ii + (jj-1)*120),0);

       
        
    end
    
    
end

%PUNY ATTEMPT AT PARALLELISATION
% maxxs  = (max(max(BouonRow)));
% maxMatrix = repmat(maxxs,size(BoutonRow,1),size(BoutonRow,2));
% BoutonImg = BoutonRow./maxMartix;

Training.Images = BoutonRow;


%% This for loop Performs the HoG

Im = BoutonIm(1:13,1:13,:);


for kk=1:size(Im,3)
    
    
    %once this actually works, reset the radius to 8
    
    B = [];
    
    B = Im(:,:,kk);
    
    %Note that this is a custom script. I don't really understand how they
    %normalise their histograms so I just won't do that.
    
    %B = imresize(B,4);
    
    % Corresponding matrix element
    %figure; imshow(B);
    
    
    
    %% Calculate the Gradients
    
    %not sure why using conv2 and atan2, think about it..
    
    grX =  [-1,0,1];
    grY = grX';
    
    %not clear on why, when I try to convolve the two, have gradient along the
    %edges by the left
    
    EdgeX = [];
    EdgeY = [];
    
    EdgeX = conv2(B,grX,'same');
    EdgeY = conv2(B,grY,'same');
    
    
    %This code is just made much more ugly, note that use unclipped EdgeX and Y
    %for calculation of the angle but then clip it and repad it for calculating
    %the magnitude of the gradient
    
    %gradient orientation, think angle (not sure if need to reverse X and Y.
    gr_rad = atan2(EdgeX, EdgeY);
    
    
    %Clip the zero padded edges
    EdgeX = EdgeX(2:end-1,2:end-1);
    EdgeY = EdgeY(2:end-1,2:end-1);
    
    
    
    %gradient magnitude think hypothenuse length of triangle with sides EdgeY and
    %EdgeX
    gr_mag = sqrt(EdgeX.^2 + EdgeY.^2);
    
    gr_mag= padarray(gr_mag,[1 1]);
    % gr_mag(:,end+1) = 0;
    % gr_mag(end+1,:) = 0;
    
    
    
    gr_deg = 180*gr_rad/pi + 90;
    
    %figure, imshow(gr_mag);
    
    
    % might be worth trying to visualise the output of this.
    %quiver(reshape(grX,,grY,)
    
    %based on my
    
    
    
    %% Time to do the binning into cells
    %The image is the subdivided into small connected regions called cells and
    %compiling a histogram of gradient directions/edge orientations within each
    %cell. Can obviously choose different sizes for each cell
    
    
    cellSize = 2;
    
    %Binwidth for angles in degrees
    orientBinW = 30;
    
    %Here 12 is the number of bins into which to around the gradient
    gr_degRound = round(gr_deg/orientBinW)*orientBinW;
    gr_radRound = pi*gr_degRound/180 ;
    
    orientBinN = 360/orientBinW;
    %gives the number of columns of cells per image. This is not a
    %generalisation, because it surely relies on the size of the cells.
    numCells = floor(length(gr_degRound)/cellSize)^2;
    numRows = sqrt(numCells);
    numColumns = sqrt(numCells);
    
    %Just initialise to sped
    orientCells = zeros(cellSize, cellSize, numCells);
    magnCells = orientCells;
    
    
    %Ok, so this creates your cells
    for jj = 0:(numRows-1)
        
        for ii = 0:(numColumns-1)
            
            orientCells(:,:,(numColumns*jj + ii +1 )) = gr_degRound(1+cellSize*jj: cellSize + cellSize*jj, 1+ii*cellSize: cellSize+ii*cellSize);
            magnCells(:,:,(numColumns*jj + ii +1 )) = gr_mag(1+cellSize*jj: cellSize + cellSize*jj, 1+ii*cellSize: cellSize+ii*cellSize);
        end
        
    end
    
    
    %at cellSize=2, you get an output where size(,3)=36, because another wont
    %fit. Think that the image is 13^2 large and each cell has 4 pixels, so
    %have 42 potential  images, however, since we are dealing with square
    %'subimages', then having 7, instead of 6, rows and columns would require
    %fitting in 49 squares
    
    %The output of the above loops is the orientation of the pixel gradient for
    %each pixel in the cell and its magnitude. It is organised so that the
    %first two dimensions are the dimensions of the cell and the third
    %dimension is the number of cells in one 'sample' image.
    
    
    %% want you want to do now is create an array, of length of the number of bins of your angles
    %and fill it in with the sum of the elements that correspond to the value
    %in the magnCells bin
    
    %First find what values are actually present in the orientation of each
    %thing
    
    
    %orientHist = zeros(49,13);
    
    %for all the cells in one image
    for ii= 1:size(orientCells,3)
        
        temp = magnCells(:,:,ii);
        
        %for all orientations in one cell
        for jj = 1:orientBinN + 1               %ERROR Here. Jwp 25th jan 2014
            
            
            
            indices = find(orientCells(:,:,ii) ==(-120 + jj*30));
            
            %convert indices into the correct for selecting both rows and
            %columns
            %[Row, Column] = ind2sub([size(orientCells)], indices);
            
            
            orientHist(ii,jj,kk) = sum(temp(indices));
            
        end
        
    end
    
    
    %The output of this
    
    
end

sizeOutp = size(orientHist);

Histf = reshape(orientHist, sizeOutp(1)*sizeOutp(2), sizeOutp(3))';

Group = svmclassify(SVMStruct,Histf);



%% This will find the circles (putative boutons) in those subimages the classifier has deemed to contain boutons
a = find(Group ==1);

plar = 1;

warning('off','all');
warning;



for i = 1:length(a)
    
    Images(:,:,i) = BoutonIm(:,:,a(i));
    
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
        %viscircles(centersBright(1,:), (radiiBright(1)+2),'EdgeColor','b','LineWidth',0.5);
        
        plar = plar+1;
        
        RegGr(:,:,i) = regiongrowing(Images(:,:,i),round(centersBright(1,2)),round(centersBright(1,1)),0.07);
        subimMask(:,:,i) = imdilate(RegGr(:,:,i),strel('disk',2));
        
        
    else
        
        RegGr(:,:,i) = zeros(size(Images(:,:,1)));
        subimMask(:,:,i) = zeros(size(Images(:,:,1)));
    end
    
    
  
end

warning('on','all');
warning

mask = subimMask;





%% Now need to create the layover mask

Mask = zeros(size(Img));

%Works, note that 120 is used because floor(ImSize(1) -
%stepSize*radius)/stepSize = 120 , where ImSize refers to the original image
%not subimages

stepSize = 4;
radius = 8;

for i = 1:size(mask,3)
    
    k = Img((1+floor(a(i)/120))*stepSize:(1+floor(a(i)/120))*stepSize + 2*radius,1+rem(a(i),120)*stepSize:1+rem(a(i),120)*stepSize+2*radius).*mask(:,:,i);
    idxSubIm = find(k);
    
    labelMean = mean2(k(idxSubIm));
    
    subimMean = mean2(Img((1+floor(a(i)/120))*stepSize:(1+floor(a(i)/120))*stepSize + 2*radius,1+rem(a(i),120)*stepSize:1+rem(a(i),120)*stepSize+2*radius));
    
    %NOTE I HAVE PUT -0.1 as an arbitrary threshold here!
    if (subimMean - labelMean) >= -0.1
        
        mask(:,:,i) = zeros(size(mask(:,:,1)));
        
    end
    
    Mask((1+floor(a(i)/120))*stepSize:(1+floor(a(i)/120))*stepSize + 2*radius,1+rem(a(i),120)*stepSize:1+rem(a(i),120)*stepSize+2*radius) = Mask((1+floor(a(i)/120))*stepSize:(1+floor(a(i)/120))*stepSize + 2*radius,1+rem(a(i),120)*stepSize:1+rem(a(i),120)*stepSize+2*radius)+ mask(:,:,i);
    
    
end

%% See if the selected area is brighter than the surround and discard otherwise

%Mask(find(Mask==1)) = 0;
finalMask = im2bw(Mask,0.8);
ROImask = bwlabel(finalMask);


%% Overlay the two

figure, 

% E is the thing you want to overlay on the original.
E = Visible_Boutons(Img,0);
imshow(E, 'InitialMag', 'fit')

I = finalMask*0.2;
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


toc
