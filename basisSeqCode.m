%% Generatin Microphone Data
lenSeq = 20e3;                      % Length of sequence
fs = 96e3;
srcTimeData = randn(lenSeq,1);

micTimeData = zeros(K,N);           % Init variable for microphone data

for k = 1:K
    micTimeData(k,:) = filteR(h(k,:),1,srcTimeData);
end

micNewTimeData = micTimeData;       % Assignment of data
%% Processing each available frame by Matched Filter

lenFrame = 501;                     % Length of frame
offsetSample = 4e3;                 % All microphones received a signal

numFrames = floor(lenSeq-offsetSample)/lenFrame;  % Number of frames to process
estX = zeros(numFrames,1); estY = zeros(numFrames,1);
tJ = cell(numFrames,1);
Jtotal = 0;
for frame = 1:numFrames
    startSample = offsetSample + (frame-1)*lenFrame;
    micData = micNewTimeData(:, startSample + 1 : startSample + lenFrame)';
    srcTimeDataT = srcTimeData(startSample + 1 : startSample + lenFrame);
    srcFreqData = getFreqMicData(srcTimeDataT,lenFrame,1);
    [estX(frame),estY(frame),~,J] = nfSequentialMLE_TOA_DOA(micData,srcTimeDataT,...
        srcFreqData,K,micPos,1,lenFrame,fs);
    tJ{frame} = J{1};
    Jtotal = Jtotal + J{1};
end
%%
Jtotal5 = Jtotal./max(Jtotal(:));
figure, imagesc(xGrid,yGrid,(Jtotal5)'), set(gca,'YDir','normal')
%% Coarse Estimate
xGrid = -3:0.1:12;             
yGrid = -3:0.1:12;
[~,indxEst] = findMaximizer(Jtotal4);
estPos = [xGrid(indxEst(1)) yGrid(indxEst(2))];
%% Refined Estimate
xF = estPos(1)-0.5:0.005:estPos(1)+0.5;
yF = estPos(2)-0.5:0.005:estPos(2)+0.5;
[estPos, e, xF, yF] = extensiveSearch(micNewTimeData,srcTimeData,micPos,...
    fs,xF,yF);
figure, mesh(xF,yF,e'), set(gca,'YDir','normal')
%% Update Data
delayData = nfGenDelayDataNonPeriodic(srcTimeData,micPos,[estPos 1.26],fs);
micNewTimeData = micNewTimeData - delayData';