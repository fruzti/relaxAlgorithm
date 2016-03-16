function J = nfEvalCostWall(micFreqData,srcFreqData, micPos,srcPos,l,N,fs,walls)
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
    
    rflxSrcs = genWallSrc(srcPos,walls)';
    numWalls = size(rflxSrcs,2);
%     tmp = zeros(size(a,1),numWalls+1);
    tmp = 0;
    for r = 1:numWalls
        a = a + nfGen_a(micPos, rflxSrcs(:,r), l, N, fs,1);
    end
    
    for k = 1 : K
        tmp = tmp + conj(a(:,k)) .* micFreqData(:,k);
%         tmp(:,1) = tmp(:,1) + conj(a(:,k)) .* micFreqData(:,k);
%         for r = 1:numWalls
%             aT = nfGen_a(micPos, rflxSrcs(:,r), l, N, fs,1);
%             tmp(:,r+1) = tmp(:,r+1) + conj(aT(:,k)) .* ...
%                 micFreqData(:,k);
%         end
    end
    
%     tmpR(1) = sum(conj(srcFreqData) .* tmp(:,1));
    
%     for r = 1:numWalls
%         tmpR(r+1) = sum(conj(srcFreqData) .* tmp(:,r));
%     end

    tmpR = sum(conj(srcFreqData) .* tmp);
    
    J = norm(tmpR)^2;
end