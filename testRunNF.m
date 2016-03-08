% Script to test estability of the ML estimate for NF source localization
% SNR - Test
SNR = -10:3:20;
numIt = 3;

e = zeros(numIt, length(SNR));

for s = 1:length(SNR)
    for it = 1:numIt
        nPow =  mean(var(micTimeData))*10^(-SNR(s)/10);
        n = sqrt(nPow) * randn(size(micTimeData));
        micTimeDataN = micTimeData + n;
        
        [estX,estY,~,~] = nfSequentialMLE_TOA_DOA(micTimeDataN,srcTimeData,...
            srcFreqData,K,micPos,1,N,fs);
        
        e(it,s) = norm(srcPos(1:2) - [estX estY]);
        
    end
    s
end

med_e = median(e);
% std_e = std(e,0,1);
plot(SNR,med_e), xlabel('SNR [dB]'), ylabel('|| xs - xsEst ||')
title('Estimation Error')
%%
% errorbar(SNR,mu_e,std_e)
% Micro - Test
numIt = 3;


posStd = [0.01:0.01:0.1];
e = zeros(numIt, length(posStd));

for s = 1:length(posStd)
    for it = 1:numIt
        nPow =  mean(var(micTimeData))*10^(-60/10);
        n = sqrt(nPow) * randn(size(micTimeData));
        micTimeDataN = micTimeData + n;
        
        nMic = posStd(s) * randn(size(micPos));
        micPosN =  micPos + nMic;
        
        [estX,estY,~,~] = nfSequentialMLE_TOA_DOA(micTimeDataN,srcTimeData,...
            srcFreqData,K,micPosN,1,N,fs);
        
        e(it,s) = norm(srcPos(1:2) - [estX estY]);
        
    end
    s
end
%%
med_e = median(e);
plot(posStd,med_e), xlabel('std [cm]'), ylabel('|| xs - xsEst ||')
title('Estimation Error')