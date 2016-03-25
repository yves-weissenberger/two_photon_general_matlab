function [corrMat] = Pop_corr_mat(popResponse)


%% Calculate the correlation matrix between responses.

%Initialise some global parameters
numCells = size(popResponse,1);
numFreqs = size(popResponse,3);




for freq1Idx=1:numFreqs
    
    for freq2Idx=1:numFreqs
        
        %Each row of respFreqn is a 
        respFreq1 = squeeze(popResponse(:,:,freq1Idx))';
        respFreq2 = squeeze(popResponse(:,:,freq2Idx))';
        
        %This gives you the single trial correlation matrix; ie element
        %(m,n) of the matrix contains the correlation between population
        %responses to the mth presentation of respFreq1 and the nth
        %presentation of the respFreq2
        trlCorrMat = ((respFreq1*respFreq2') - numCells*mean(respFreq1,2)*mean(respFreq2,2)')./((numCells-1)*std(respFreq1,0,2)*std(respFreq2,0,2)');
        
        %This gives you the correlation matrix, ie the correlation between
        %responses to each of the stimuli
        corrMat(freq1Idx,freq2Idx) = median(trlCorrMat(:));
        
        
        
    end
    
end

end