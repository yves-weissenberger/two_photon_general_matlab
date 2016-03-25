%% Select a subpart of your image

clear all

%load('/Users/Yves/Desktop/Bouton_Grabinfo_files/Area01/(20130514_15_01_53)-_Area1_73_noise_GRABinfo.mat')
load('/Users/Yves/Documents/MATLAB/Training_Cells.mat')
Im = Training.Images;

%Im = BoutonRow(1:13,1:13,:);


for kk=1:size(Im,3)


%once this actually works, reset the radius to 8
radius = 13;

B = [];

B = Im(:,:,kk);

%Note that this is a custom script. I don't really understand how they
%normalise their histograms so I just won't do that.
B = Visible_Boutons(B);

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


%% This is just a very short control script for mixing up the labels to see how the chance classifier performs
%idx = randperm(length(Training.Boutonlabel));
%Training.Boutonlabel = Training.Boutonlabel(idx);

%%

SVMStruct = svmtrain(Histf,Training.BoutonLabel);

%Group = svmclassify(SVMStruct,Histf);