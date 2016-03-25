 %This will 
 a = find(Group ==1);
 
 plar = 1;
 
 for i = 1:50
     
Images(:,:,i) = Visible_Boutons(Training.Images(:,:,a(i)),1);
     
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
    % viscircles(centersBright(1,:), (radiiBright(1)+2),'EdgeColor','b','LineWidth',0.5);
     
     plar = plar+1;
 end

    RegGr(:,:,i) = regiongrowing(Images(:,:,i),round(centersBright(1,2)),round(centersBright(1,1)),0.10)
    subimMask = imdilate(RegGr(:,:,i),strel('disk',2))
    
    %figure, imshow

 end