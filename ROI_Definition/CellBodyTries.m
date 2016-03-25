

for j = 1:size(I3,3)
% CentreCircle_big = sqrt((rr-CircCentre_small).^2+(cc-CircCentre_small).^2)<=3;
% aa = I3(:,:,j).*CentreCircle_big;

% afind = find(aa);
% medaa = median(aa(afind));

b(:,:,j) = regiongrowing(I3(:,:,j),13, 13,.050);

% if length(find(b(:,:,j))) >= 150
%     
%     [CentersDark, RadiiDark] = imfindcircles(I3(:,:,j),[2 8], 'Sensitivity',1,'Objectpolarity', 'dark','Edgethreshold',0.001,'Method','Twostage');
% 
% if RadiiDark ~= 0
% centersDark(j,1) = CentersDark(1,1);
% centersDark(j,2) = CentersDark(1,2);
% radiiDark(j) = RadiiDark(1);
% end


% % Create the corresponding mask
% [m, n, o] = size(I3);
% 
% % Select the annulus of supposed bright pixels.
% [rr, cc] = meshgrid(1:m);


    
%b(:,:,j) = sqrt((rr-centersDark(j,1)).^2+(cc-centersDark(j,2)).^2)<=radiiDark(j) + 4;



% end



end


warning('on','all');
warning;
    
figure, imshow(b)
figure, imshow(I3(:,:,j))





