function [Xpil, Ypil] = ICNeuropilMask(GRABinfo)

ROICentres = reshape([(GRABinfo.ROIprops(:).Centroid)],2,length([(GRABinfo.ROIprops(:).Centroid)])/2)';


lengthEveryThing = max(max(GRABinfo.ROImaps));
radius = 13;

%% Create the Neuropil Mask

%as output basically want an extra GRABinfo thing that has the [y x]
%coordinates of the ROImasks for each cell

Ypil = zeros((2*diameter)^2,lengthEveryThing);
Xpil = zeros((2*diameter)^2,lengthEveryThing);

Xxpil = ROICentres(:,1);
Yypil = ROICentres(:,2);

% Check that everything is within bounds
toosmall_xpil = find(Xx<=2*radius);
Xxpil(toosmall_xpil) = 2*radius+1;

toosmall_ypil = find(Yy<=2*radius);
Yypil(toosmall_ypil) = 2*radius+1;

toobig_xpil = find(Xx>=(xmax - 2*radius));
Xxpil(toobig_xpil) = (xmax - 2*radius) - 1;

toobig_ypil = find(Yy>=(xmax - 2*radius));
Yypil(toobig_ypil) = (ymax - 2*radius) - 1;


%This could go wrong, if x and y are zero for a given lookup.
for jj = 1:lengthEveryThing
    
    
  empty_Mask = zeros(size(Final_Mask));
  
  cellpill = imdilate(GRABinfo.ROImaps(Yypil(jj)-2*radius:Yypil(jj)+2*radius,Xxpil(jj)-2*radius:Xxpil(jj)+2*radius),strel('disk',4));
  cellpill = imfill(cellpill,'holes');
  neuropil(:,:,jj) = ones(size(cellpill)) - cellpill;
  
  %This allows easy lookup of the x and y positions of the neuropil mask
  empty_Mask(Yypil(jj)-2*radius:Yypil(jj)+2*radius,Xxpil(jj)-2*radius:Xxpil(jj)+2*radius) = neuropil(:,:,jj);
  
  [y,x] = find(empty_Mask==1);

  Ypil(1:length(y),jj) = y;
  Xpil(1:length(x),jj) = x;
  
end

