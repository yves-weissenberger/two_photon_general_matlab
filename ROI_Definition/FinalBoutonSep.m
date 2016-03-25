

agg = finalMask.*Img;
agg2 = agg;

numROIs = max(max(ROImask));

for ii = 1:numROIs
    
    ROIpix = find(ROImask==ii);
    
    if length(ROIpix)>=70
        
        ROImean = mean2(agg2(ROIpix));
        ROImax = max(max(agg2(ROIpix)));
        
        Thresh(ii) = (0.9*ROImean + 1.1*ROImax)/2;
        
        
        agg2(ROIpix) = im2bw(agg2(ROIpix),Thresh(ii));
        
    end
    
end

agg3 = bwlabel(agg2);

ROIcolour = label2rgb(agg3,'jet',[0.2,0.2,0.2],'shuffle');


OrgImg = imshow(1.6*Visible_Boutons(GRABinfo.AllMeanStack,0));

hold on
Colourshow = imshow(ROIcolour);
hold off

set(Colourshow, 'AlphaData', OrgImg)


%Overlay(Img,agg2,0.2);