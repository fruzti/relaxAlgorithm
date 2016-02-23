function a = nfGen_a(micPos,srcPos,l,N,fs,flag)
% function a = nfGen_a(micPos,srcPos)
% ------------------------------------
% micPos : position of microphones
% srcPos : position of the source
% l : number of harmonics
% N : number of samples
% fs : sampling frequency

    if nargin < 6
        flag = 0;
    end

    % get attenuation and delay information
    [~, kBeta, kTau] = getMic2SrcParams(micPos, srcPos);

    % delay in number of samples
    eta = kTau*fs;
    
    % number of microphones
    K = length(kBeta);
    
    if flag
%         kBeta = 1./kBeta;
        kBeta = ones(K,1);
    end
    
    a = zeros(N,K);
    for k = 1:K
        a(:,k) = kBeta(k) .* nfGen_d(eta(k), l, N);
    end
    
end