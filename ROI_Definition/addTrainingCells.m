function addTraining(targetDir, Area, label)


%eg addTraining('/Users/Yves/Desktop/Bouton_Grabinfo_files','Area03',0)

%has been trained on 20130514 Area 4
clearvars -except targetDir Area label


%dir returns all files in a given directory. f = fullfile('myfolder','mysubfolder','myfile.m')
%returns a string combining all the inputs
Area_files = dir([fullfile(targetDir,Area) '/*.mat']);

%Load image to be processed
GRABname=Area_files(1).name;
load(fullfile(targetDir,Area,GRABname));

%Deafults accepting output of the function to no
acceptstr='n';



%% Preprocessing.

A = GRABinfo.AllMeanStack;


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
Im = A_I3;

roifig=figure('Position',get(0,'ScreenSize')); imshow(Im);

ii = 0;

radius = 13;

%% Select the Cells

while ishandle(roifig)


% Select your circle
[x, y] = ginput(1);

ii = ii + 1;


if x ~= 0
    
Xx(ii) = round(x);
Yy(ii) = round(y);
viscircles([x, y], (radius),'EdgeColor','b','LineWidth',0.2);
end

if ii - length(Xx) > 0    
close gcf
end

x = [];
y = [];

end


%%

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
 
%% Add to previous training data


load('Training_Cells.mat')

if label==1
    
    labels = ones(length(A),1);
    
elseif label==0
    
    labels = zeros(length(A),1);
    
end

test1 = cat(3,Training.Images,A);
test2 = cat(2,Training.BoutonLabel,labels');

idx = randperm(length(test2));

Training.BoutonLabel = test2(idx);

Training.Images = test1(:,:,idx);

save('Training_Cells','Training')





