function delaySignal = nfGenDelayDataNonPeriodic(srcTimeData, micPos,...
    srcPos,fs)
    N = length(srcTimeData);
    K = size(micPos,1);
    [~,~,delay] = getMic2SrcParams(micPos',srcPos');
    delaySignal = zeros(N,K);
    for k = 1:K
        delaySample = round(delay(k)*fs);
        if mod(delaySample,2) == 0
            delaySample = delaySample + 1;
        end
        extSignal = [srcTimeData; zeros(delaySample,1)];
        tmp = nfGenDelayData(extSignal, micPos(k,1:2)',...
            srcPos(1:2)',fs);
        delaySignal(:,k) = tmp(1:N);
    end

end