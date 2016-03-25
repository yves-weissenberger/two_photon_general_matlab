function [reliableSounds,reliableCorrMat,ClustTree,order,clustStore] = Clusters_from_CorrMat(corrMat,reliableOnly,thresh,toPlot)


if reliableOnly==true
    
    
    %% This will select then only the reliable responses
    
    %The method for doing this is to look for population responses to one
    %stimulus that have a correlation above a certain value
    
    %First select the correlations between responses to the same stimulus
    D = diag(corrMat);
    
    
    %find the indices of the sounds that are reliable
    reliableSounds = find(D>=thresh);
    %find the indices of the sounds that are unreliable
    unreliableSounds = find(D<=thresh);
    
    
    %This contains the number of reliable sounds at each threshold
    numReliableSounds = length(reliableSounds);

    %reduce the correlation matrix to only those pairs of sounds that are
    %reliable
    reliableCorrMat = corrMat([reliableSounds],[reliableSounds]);
    
    %This returns returns a matrix Zreliable that encodes a tree of hierarchical 
    %clusters of the rows of the real matrix reliableCorrMat. Z is a 
    %(m ? 1)-by-3 matrix, where m is the number of reliable frequencies.
    %The first two columns contain the nodes that are linked, the third
    %contains the linkage distance. 
    ClustTree = linkage(reliableCorrMat,'complete');
    
    %Create the dendrogram by magic. orderR puts the order of stimuli in
    %the correlation matrix to the right place
    [Hreliable,Treliable,order] = dendrogram(ClustTree,0);
    
    %This closes the dendrogram plot if not plotting
%     if toPlot==0
%         close(gcf)
%     end
    
   
        
    %% This just goes a bit of plotting
        %Make some variables pretty for plotting
        freqOrderR = ceil(reliableSounds./4);
        freqLevelR = rem(reliableSounds,4)./10;
        clustStore = reliableCorrMat(order,order);

    if toPlot==true

        %Do the actual plots
        figure('Units','Normalized'),
        subplot(1,2,1)
        pcolor(flipud(reliableCorrMat'));

        ax = gca;
        ax.XTickLabel = flipud(freqOrderR+freqLevelR);
        ax.XTick = 1:length(freqOrderR);
        ax.YTick = 1:length(freqOrderR);
        ax.YTickLabel = freqOrderR+freqLevelR;
        axis square
        subplot(1,2,2)
        pcolor(flipud(reliableCorrMat(order,order)'))
        ax = gca;
        ax.XTickLabel = flipud(freqOrderR(order)+freqLevelR(order));
        ax.XTick = 1:length(freqOrderR);
        ax.YTick = 1:length(freqOrderR);
        ax.YTickLabel = freqOrderR(order)+freqLevelR(order);
        ax.FontSize = 7.5;
        ax.FontWeight = 'bold';
        colorbar()
        axis square
    end
    
    
elseif reliableOnly==0
    
    %%
    
    
    %plot the images based on all responses
    ClustTree = linkage(corrMat,'single');
    [H,T,OUTPERM] = dendrogram(ClustTree,0);
    
    order = OUTPERM;
    %order = fliplr(OUTPERM);
    clustStore = corrMat(order,order);
    
    if toPlot==true
        figure,
        %imshow(corrMat,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
        pcolor(corrMat')
        colorbar()
        figure,
        %imshow(clustMeanStore,[0, maxColourMap],'InitialMagnification',500,'Colormap',jet(255))
        pcolor(clustStore')
        colorbar()
    end
    

end




end