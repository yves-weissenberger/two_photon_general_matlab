function I2 = Visible_Boutons(I,show)

num = nargin;

if nargin == 1
    
    show = 0;
    
end

% Area_files = dir([fullfile(targetDir,Date,Area) '/*.mat']);
% 
% load('/Users/Yves/Desktop/Boutons_preprocessed/20140116/Area02/(20140116_15_16_17)-_MGB20140116_Area02_tones1_GRABinfo.mat')

%I = GRABinfo.AllMeanStack;
normalising_factor = max(max(I));
I = I./normalising_factor;

I2 = imadjust(I);

if show == 1
figure, imshow(I2)
end

end