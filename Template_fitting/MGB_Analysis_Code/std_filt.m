

% from scipy.ndimage.filters import uniform_filter
% 
% def window_stdev(arr, radius):
%     c1 = uniform_filter(arr, radius*2, mode='constant', origin=-radius)
%     c2 = uniform_filter(arr*arr, radius*2, mode='constant', origin=-radius)
%     return ((c2 - c1*c1)**.5)[:-radius*2+1,:-radius*2+1]


kernel = ones(3,3,3);


convn(


%% Notes

%important that this works to see that the standard deviation of a random
%variable is given by E[(X - E[X])^2] which is mathematically the same as
%E[X^2] - E[X]^2 which means can use a uniform filter