% Test for estimating the virtual sources in a room.

run experimentalSetup.m

% Simulation Parameters
zPos = 1.26;
srcPos = [ladLoud' zPos];
srcPos(:,1) = srcPos(:,1) + 4;
micPos = [arrayElems; ones(1,size(arrayElems,2))*zPos]';
micPos(:,1) = micPos(:,1) + 4;

arrayPos = mean(micPos(7:12,:));

clearvars -except arrayPos srcPos micPos

vrtSrcPos = genSrcsFromWalls(srcPos(1:2), [4 6]);

% MLS sequence
srcTimeData = mls(7,5); 
% Num of samples after the signal reaches the last array
N1 = 127*2+1;

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