%% Generatin Microphone Data
lenSeq = 35e3;  % Length of sequence
N = lenSeq; 
fs = 96e3;
% srcTimeData = randn(lenSeq,1);
srcTimeData = [mls(15,1); 0]; lenSeq = length(srcTimeData); N = lenSeq;
K = 18;
kVec = [randperm(6,K/3) randperm(6,K/3)+6 randperm(6,K/3)+12];
% kVec = randperm(18,K);

micTimeData = zeros(K,N);           % Init variable for microphone data

for k = 1:K
    micTimeData(k,:) = filter(h(kVec(k),:),1,srcTimeData);
end
micPosK = micPos(kVec,:);
%%
micNewTimeData = micTimeData;       % Assignment of data
lenFrame = 1001;                     % Length of frame
offsetSample = 4e3;                 % All microphones received a signal

numFrames = floor((lenSeq-offsetSample)/lenFrame);  % Number of frames to process
xGrid = -3:0.1:9;             
yGrid = -3:0.1:12;

micNewTimeData = micTimeData;       % Assignment of data
maskMax = ones(length(xGrid),length(yGrid));

numWallsAndSource = 5;
estPos(1,:) = srcPos(1:2);
%% Main Loop
for l = 1 : numWallsAndSource
    if l > 1
    str2Disp = strcat('Estimation of Source #: ',num2str(l));
    disp(str2Disp)
    %% Processing each available frame by Matched Filter
    estX = zeros(numFrames,1); estY = zeros(numFrames,1);
    tJ = cell(numFrames,1);
    Jtotal = 0;
    disp('----------------------------')
    disp('Starting Coarse Estimate...')
    disp('----------------------------')
    for frame = 1:numFrames
        startSample = offsetSample + (frame-1)*lenFrame;
        micData = micNewTimeData(:, startSample + 1 : startSample + lenFrame)';
        srcTimeDataT = srcTimeData(startSample + 1 : startSample + lenFrame);
        srcFreqData = getFreqMicData(srcTimeDataT,lenFrame,1);
        [estX(frame),estY(frame),~,J] = nfSequentialMLE_TOA_DOA(micData,srcTimeDataT,...
            srcFreqData,K,micPosK,1,lenFrame,fs);
        tJ{frame} = J{1};
        Jtotal = Jtotal + J{1};
    end
    estVrtPos{l} = [estX estY];
    Jt{l} = Jtotal./max(Jtotal(:));
    [~,indxEst] = findMaximizer(Jt{l},maskMax);
    estPos(l,:) = [xGrid(indxEst(1)) yGrid(indxEst(2))];
    cEstPos(l,:) = estPos(l,:);
    disp('Done Coarse Estimate.')
    disp(cEstPos(l,:))
    disp('---------------------')
    %% Refined Estimate
    disp('----------------------------')
    disp('Starting Refinement.........')
    disp('----------------------------')
    if l > 1
        sector = findSectorOfPoint(estPos(1,:),estPos(l,:));
        if (sector == 1 || sector == 3)
            xF = estPos(l,1)-0.5:0.005:estPos(l,1)+0.5;
            yF = estPos(1,2)-0.01:0.005:estPos(1,2)+0.01;
        elseif (sector == 2 || sector == 4)
            yF = estPos(l,2)-0.5:0.005:estPos(l,2)+0.8; 
            xF = estPos(1,1)-0.01:0.005:estPos(1,1)+0.01;
        end
    else
        xF = estPos(l,1)-0.5:0.005:estPos(l,1)+0.5;
        yF = estPos(l,2)-0.5:0.005:estPos(l,2)+0.5; 
    end
%     xF = estPos(l,1)-0.5:0.005:estPos(l,1)+0.5;
%     yF = estPos(l,2)-0.5:0.005:estPos(l,2)+0.5; 
    [estPos(l,:), e, xF, yF] = extensiveSearch(micNewTimeData,srcTimeData,micPosK,...
        fs,xF,yF);
    disp('Done Refinement')
    end
    %% Fixing Mask
    if l > 1
        switch (sector)
            case 1
                maskMax(xGrid > 1.5,:) = 0;
            case 2
                maskMax(:,yGrid > 0.7) = 0;
            case 3
                maskMax(xGrid < 0.4,:) = 0;
            case 4
                maskMax(:,yGrid < 0) = 0;
        end
    else
        maskMax = zeros(length(xGrid),length(yGrid));
        maskMax(xGrid > 0.3 & xGrid < 1, :) = 1;
        maskMax(:, yGrid > -0.2 & yGrid < 0.9) = 1;
        maskMax(xGrid > 0.3 & xGrid < 8,yGrid > -0.1 & yGrid < 9) = 0;
    end
    %% Update Data
    delayData = nfGenDelayDataNonPeriodic(srcTimeData,micPosK,[estPos(l,:) 1.26],fs);
    micNewTimeData = micNewTimeData - delayData';
    disp('Current Estimate')
    disp('----------------------------')
    disp(estPos(l,:))
    disp('----------------------------')
end
plotRoom(micPos, srcPos,vrtPos', estPos')
disp('Microphone Used')
disp(kVec)
disp('Array')
disp(ceil(kVec/6))