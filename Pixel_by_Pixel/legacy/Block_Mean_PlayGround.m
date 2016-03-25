%% Calculate Selectivity



%inter-stimulus-interval
ISI = 45;

%calculate the cumulative sum
A = cumsum(uint32(Pixeltraces),3);


%Find the cumulative sum at the presentation of each stimulus
B = A(:,:,1:45:end);

%For memory safekeeping
clear A


%calculate the response sum in each window by subtracting the cumulative sum
%from the previous window (45 frames prior)
C = (cat(3,B,zeros(size(B,1),size(B,2),1))- cat(3,zeros(size(B,1),size(B,2),1),B));

%convert to a mean
C = C./ISI;

%Load the trial sequence
load('tones_outDat.mat')
stim_seq = outDat.trialOrder;

%Sound sequence is repeated 3x so need to repeat the matrix 3x
allTrialOrder = repmat(stim_seq,3,1);


%Ok, so I think stimuli 1-4 are the same frequency but different levels, so
%should get them all into multiples of 4...
freqOrder = ceil(allTrialOrder./4);
%freqOrder = mod(allTrialOrder,25);
%%

%initialise matrix with each pixel in the first two dimensions, the
%stimulus_index in the 3rd and the repetitions in the fourth
%D = zeros(size(B,1),size(B,2),100,9,'uint32');

for stim_idx=1:25
    
    D(:,:,stim_idx,:) = C(:,:,freqOrder==stim_idx);

end



E = mean(D,4);

%% Find the mean and max of E

mean_val = mean(E,3);

[max_val,max_idx] = max(E,[],3);


%This is the sound selectivity index
F = (max_val-mean_val)./(max_val+mean_val);
F = F./(max(max(F)));
F2 = histeq(F);


%% Now create the plots

%One color for each frequency
Color_list = jet(25);
%clrs = reshape(Color_list,25,1,3);

%calculate the range of intensities and scale between 0 and 1
top_val = max(max(F));
bot_val = min(min(F));




%%

final_Img1 = zeros(512,32);
final_Img2 = zeros(512,32);
final_Img3 = zeros(512,32);


for ii=1:25
    
    %[x,y] = ind2sub([512,32],find(max_idx==ii));
    idxs = find(max_idx==ii);
    final_Img1(idxs) = Color_list(ii,1);
    final_Img2(idxs) = Color_list(ii,2);
    final_Img3(idxs) = Color_list(ii,3);
end

final_Img = cat(3,final_Img1,final_Img2,final_Img3);
final_Img = final_Img.*repmat(F,1,1,3);

%%

enhanced_meanStack= Visible_Boutons(mean(Pixeltraces,3),0);

figure('Units','normalized')
%image(final_Img)

imshow(cat(2,repmat(enhanced_meanStack,1,1,3),zeros(512,10,3),final_Img))

%background = zeros(size(max_idx));


%h = imshow(background); hold on

%imshow(max_idx,'Colormap',Color_list); hold off

%opacity = 0;
%set(h, 'AlphaData', opacity*F)


