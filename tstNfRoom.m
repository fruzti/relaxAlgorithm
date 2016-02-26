% Test for estimating the virtual sources in a room.

run experimentalSetup.m

% Simulation Parameters
zPos = 1.26;
srcPos = [ladLoud' zPos];
srcPos(:,1) = srcPos(:,1) + 4;
micPos = [arrayElems; ones(1,size(arrayElems,2))*zPos]';
micPos(:,1) = micPos(:,1) + 4;
refMics = [singleRef pairRef]';
refMics(:,1) = refMics(:,1) + 4;

arrayPos = mean(micPos(7:12,:));

clearvars -except arrayPos srcPos micPos refMics

vrtSrcPos = genSrcsFromWalls(srcPos(1:2), [4 6]);

% MLS sequence
srcTimeData = mls(7,5); 
% srcFreqData1 = [zeros(10,1); exp(1i*2*pi*randn(53,1))];
% srcFreqData = [srcFreqData1; 0 ; flipdim(conj(srcFreqData1),1)];
% srcTimeData = applyIFFT(srcFreqData);
% srcTimeData = repmat(srcTimeData,5,1);
% Num of samples after the signal reaches the last array
N1 = 127*2+1;
% N1 = 381;

fileNameRIRs = 'ISM_RIRs.mat';
AuData = ISM_AudioData(fileNameRIRs, srcTimeData);
srcTimeData = srcTimeData(N1:end);
micTimeData = AuData(N1:end,1:18); N = size(micTimeData,1); fs = 10e3;
K = size(micPos,1);

micFreqData = getFreqMicData(micTimeData, N, K);
srcFreqData = applyFFT(srcTimeData, N);
L = 4;
%%

warning off
[estX,estY,estBeta,J] = nfSequentialMLE_TOA_DOA(micTimeData,srcTimeData,...
    srcFreqData,K,micPos,L,N,fs);

% xGrid = 0:0.01:4;
% yGrid = 0:0.01:0.5;
% xGrid = -1:0.1:8;
% yGrid = 0:0.1:12;
xGrid = -3:0.1:8;
yGrid = -3:0.1:12;

figure, imagesc(xGrid,yGrid,J{1}'), set(gca,'YDir','normal')
title('Cost Function for Source Location')
xlabel('X-Axis'), ylabel('Y-Axis')
disp('Estimate')
estimates = [estX estY];
[~,or] = sort(estimates(:,1));
disp(estimates(or,:))
disp('True')
[~,or] = sort(vrtSrcPos(:,1));
disp(vrtSrcPos(or,:))
%%
srcList = runMethodTest(micFreqData, srcTimeData, srcFreqData,...
    micPos, fs, N, K, xGrid, yGrid);
% RELAX-Like
% srcList = runMethodTest2(micFreqData, srcTimeData, srcFreqData,...
%     micPos, fs, N, K, xGrid, yGrid);
%%
yDelay = [];
yDelayF = [];
finalList = [[0.9859 0.3984]; srcList];
for l = 1:(size(srcList,1)+1)
    [~, beta, delay] = getMic2SrcParams(micPos(:,1:2)', finalList(l,:)');
    tmp = [];
    for k = 1:K
        tmp(:,k) = beta(k)*delayPeriodicSignal(srcTimeData,round(delay(k)*fs));
        tmpF(:,k) = beta(k)*delayPeriodicSignal(srcTimeData,(delay(k)*fs));
    end
    yDelay(:,l) = tmp(:);
    yDelayF(:,l) = tmpF(:);
end
bESt = (yDelay(:,2:end)\(micTimeData(:)-yDelay(:,1)));
indx = (abs(bESt) > 0.01) .* sign(bESt) ; 
estY = yDelay(:,1) + yDelay(:,2:end)*indx;
estYF = yDelayF(:,1) + yDelayF(:,2:end)*indx;

figure, plot(estY), hold on, plot(micTimeData(:),'--r')
figure, plot(estYF), hold on, plot(micTimeData(:),'--r')

SegSNR = segsnr(micTimeData(:),estY,fs)
SegSNRF = segsnr(micTimeData(:),estYF,fs)
%%
yDelay = []; yDelayF = [];
cS = [finalList(1,:); c];
for l = 1:(size(srcList,1)+1)
    [~, beta, delay] = getMic2SrcParams(refMics', finalList(l,:)');
    tmp = []; tmpF = [];
    for k = 1:size(refMics,1)
        tmp(:,k) = beta(k)*delayPeriodicSignal(srcTimeData,round(delay(k)*fs));
        tmpF(:,k) = beta(k)*delayPeriodicSignal(srcTimeData,(delay(k)*fs));
    end
    yDelay(:,l) = tmp(:);
    yDelayF(:,l) = tmpF(:);
end
% estY = sum(yDelay,2); estYF = sum(yDelayF,2);
estY = yDelay(:,1) + yDelay(:,2:end)*indx;
estYF = yDelayF(:,1) + yDelayF(:,2:end)*indx;
refSignal = AuData(N1:end,19:end);
figure, plot(estY), hold on, plot(refSignal(:),'--r')
figure, plot(estYF), hold on, plot(refSignal(:),'--r')

SegSNR = segsnr(refSignal(:),estY,fs)
SegSNRF = segsnr(refSignal(:),estYF,fs)
%%
[~, beta, delay] = getMic2SrcParams(
%%
% TO FIX
nSrcList = unique(srcList,'rows');

A = zeros(N*K,size(nSrcList,1) + 1);

tmp = nfGen_a(micPos(:,1:2)', [0.9859 0.3984]',floor(N/2),N,fs);
A(:,1) = tmp(:) .* repmat(srcFreqData,K,1);

for l = 1:size(nSrcList,1)
    tmp = nfGen_a(micPos(:,1:2)',nSrcList(l,:)',floor(N/2),N,fs);
    A(:,l+1) = tmp(:) .* repmat(srcFreqData,K,1);
end

estS = A\micFreqData(:);