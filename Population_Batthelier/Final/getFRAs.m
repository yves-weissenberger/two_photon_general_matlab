function [fraMeanStore,fraStdStore,pVals,fraMeanAll,fraStdAll] = getFRAs(popResponse)


numCells = size(popResponse,1);
numSign = 0;

signIdxs = zeros(numCells,1);

for cellIdx=1:numCells
    pVals(cellIdx)=anova1(squeeze(popResponse(cellIdx,:,:)),[],'off');
    
    FRA=flipud(reshape(mean(popResponse(cellIdx,:,:),2),4,25));
    fraMeanStore(cellIdx,:,:) = FRA;
    fraMeanAll(cellIdx,:)=FRA(:);


    fraStd = flipud(reshape(std(popResponse(cellIdx,:,:),[],2),4,25));
    fraStdStore(cellIdx,:,:) = fraStd;
    fraStdAll(cellIdx,:) = fraStd(:);

    
    for k=1:4
        smFRA(k,:)=smooth(FRA(k,:),3);
    end
        smFRAstore(cellIdx,:,:) = smFRA;
        
    if pVals<0.005
        signIdxs(cellIdx) = 1;
        numSign = numSign+1;
    end
    
end


end
