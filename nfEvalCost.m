function J = nfEvalCost(micFreqData,srcFreqData, micPos,srcPos,l,N,fs)
% function J = nfEvalCost(micFreqData,srcFreqData, micPos,srcPos,l,N,fs)
%----------------------------------------------------
% micFreqData : matrix containing the FFT of the data in the array
% srcFreqData : matrix containing the FFT of the source signal
% micPos : position of microphones
% srcPos : position of source
% l : number of harmonics
% N : number of samples
% fs : sampling frequency

    a = nfGen_a(micPos, srcPos, l, N, fs,1);
    
    tmp = 0;
    
    K = size(a,2);
    
    for k = 1 : K
        tmp = tmp + conj(a(:,k)) .* micFreqData(:,k);
    end
    
    tmp = sum(conj(srcFreqData) .* tmp);
    
    J = norm(tmp)^2;
end