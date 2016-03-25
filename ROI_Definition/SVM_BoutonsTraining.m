%% Select a subpart of your image


load('/Users/Yves/Documents/MATLAB/Training.mat')



Im = Training.BoutonRow;

for kk=1:size(Im,3)


%once this actually works, reset the radius to 8
radius = 8;

B = Im(:,:,kk);

%Note that this is a custom script. I don't really understand how they
%normalise their histograms so I just won't do that.
B = Visible_Boutons(B);

% Corresponding matrix element
%figure; imshow(B);



%% Calculate the Gradients

%not sure why using conv2 and atan2, think about it..

grX =  [-1,0,1];
grY = grX';

%not clear on why, when I try to convolve the two, have gradient along the
%edges by the left
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


cellSize = 4;

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


%need to automate the size of these things
for jj = 0:(numRows - 1)

for ii = 0:(numColumns - 1)
    
    orientCells(:,:,(numColumns*jj + ii +1 )) = gr_degRound(1+cellSize*jj: cellSize + cellSize*jj, 1+ii*cellSize: cellSize+ii*cellSize);
    magnCells(:,:,(numColumns*jj + ii +1 )) = gr_mag(1+cellSize*jj: cellSize + cellSize*jj, 1+ii*cellSize: cellSize+ii*cellSize);
end

end


%The output of the above loops is the orientation of the pixel gradient for
%each pixel in the cell and its magnitude. It is organised so that the
%first two dimensions are the dimensions of the cell and the third
%dimension is the number of cells in one 'sample' image.
    
    
%% want you want to do now is create an array, of length of the number of bins of your angles
%and fill it in with the sum of the elements that correspond to the value
%in the magnCells bin

%First find what values are actually present in the orientation of each
%thing


orientHist = zeros(49,13);

for ii= 1:size(orientCells,3)
    
    %have to add 1 for 0
    for jj = 1:orientBinN + 1    %ERROR Here. Jwp 25th jan 2014

        
        
        indices = find(orientCells(:,:,ii) ==(-120 + jj*30));
        
        orientHist(ii,jj,kk) = sum(magnCells(indices));
        
    end 
    
end


%The output of this


end

sizeOutp = size(orientHist);

Histf = reshape(orientHist, sizeOutp(1)*sizeOutp(2), sizeOutp(3))';


SVMStruct = svmtrain(Histf,Training.BoutonLabel);

Group = svmclassify(SVMStruct,Histf)