%function dF=cal c_dF(input)

%Johannes 2012
% 
% % calculate F/F0
% input=input-min(input);
% 
% tmp=prctile(input,[10 70],2);
% 
% 
% low_tc=tmp(1);
% high_tc=tmp(2);
% 
% ind=input>low_tc & input<high_tc;
% 
% 
% F0=median(input(ind),2);
% 
% dF=(input-F0)/F0;
% 


%% First Calculate a Time dependent deltaF

%set the time window for time dependent F0 calculation to 5s

%This is tf0 in Frames
t_f0= 2*375;

testarray = Pixeltraces;

%Now lets do it first for a 1-D array.
trace = squeeze(testarray(1,1,:));

%This is the f0 to subtract from each time point, want to interpolate this
%with the same value 375 times each
f0_vals = squeeze(min(reshape(trace,1,[],t_f0),[],3));


f2 = reshape(repmat(f0_vals,t_f0,1),1,40500);



%% This does it in parallel


testarray = Pixeltraces;
t_f0= 375;

tic
[d1,d2,d3] = size(testarray);

Pixels_f0 = squeeze(min(reshape(testarray,d1,d2,t_f0,d3./t_f0),[],3));

temp1 = repmat(Pixels_f0,[1,t_f0,1]);
f0_mtx = reshape(temp1,d1,d2,[]);

clear temp1
%cat(2,testarray,f0_mtx)

toc
%%
%This is the control plot
 iii = 423; jjj = 14;
 plot(squeeze(Pixeltraces(1,23,:))); hold on; plot(squeeze(f0_mtx(1,23,:)),'r')

 
 %% look at variations in df/f
 


%% Final Parallel Version with pseudo smoothing -Code 5


testarray = Pixeltraces;

t_f01 = 25;
t_f02= 750;


tic
[d1,d2,d3] = size(testarray);

%Pixels_f0 = squeeze(min(reshape(testarray,d1,d2,t_f0,d3./t_f0),[],3));

%This works
interm1 = permute(squeeze(reshape(testarray,d1,d2,t_f02,d3./t_f02)),[1,2,4,3]);
clear testarray

interm2 = reshape(interm1,d1,d2,d3./t_f02,t_f02/t_f01,t_f01);
%clear interm1
interm3 = min(mean(interm2,5),[],4);
%clear interm2

%this is the old one
%temp1 = repmat(Pixels_f0,[1,t_f0,1]);

temp1 = repmat(interm3,[1,t_f02,1]);
%clear interm3
f0_mtx = reshape(temp1,d1,d2,[]);
%clear temp1
%cat(2,testarray,f0_mtx)

toc
 
 
%% Test to see that it (=Code 5) works

for i = 1:10000
    ii = randi(size(interm3,1)); jj = randi(size(interm3,2)); kk = randi(size(interm3,2));
    test1 = squeeze(interm1(ii,jj,kk,:));
    test2 = reshape(test1, t_f02./t_f01,t_f01);
    testFinal(i) = min(mean(test2,2)) - interm3(ii,jj,kk)~=0;
end

if sum(testFinal)==0
    'The test was a success'
else
    'you suck'
end


%%


 mean_vals = mean(f0_mtx,3);
 
 bar3(mean_vals)
 
%% Finally should be able to perform relatively fast filtering by creating one long time matrix that is zero padded




