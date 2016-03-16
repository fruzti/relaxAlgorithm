function J = nfEvalCostWall(micFreqData,srcFreqData, micPos,srcPos,l,N,fs,wall)
%function J = nfEvalCostWall(micFreqData,srcFreqData, micPos,srcPos,l,N,fs,wall)
%----------------------------------------------------
% micFreqData : matrix containing the FFT of the data in the array
% srcFreqData : matrix containing the FFT of the source signal
% micPos : position of microphones
% srcPos : position of source
% l : number of harmonics
% N : number of samples
% fs : sampling frequency
% wall : wall [type and location] vector

    a = nfGen_a(micPos, srcPos, l, N, fs,1);
    
    K = size(a,2);
    
    rflxSrcs = genWallSrc(srcPos,walls);
    numWalls = size(rflxSrcs,1);
    tmp = zeros(size(a,1),numWalls+1);
    
    for k = 1 : K
        tmp(:,1) = tmp(:,1) + conj(a(:,k)) .* micFreqData(:,k);
        for r = 2:numWalls
            aT = nfGen_a(micPos, rflxSrcs(r,:), l, N, fs,1);
            tmp(:,r) = tmp(:,r) + conj(aT(:,k)) .* ...
                micFreqData(:,k);
        end
    end
    
    tmpR(1) = sum(conj(srcFreqData) .* tmp(:,1));
    
    for r = 2:numWalls
        tmpR(r) = sum(conj(srcFreqData) .* tmp(:,r));
    end
    
    
    J = norm(tmpR)^2;
end