%% Generatin Microphone Data
lenSeq = 20e3;                      % Length of sequence
fs = 96e3;
srcTimeData = randn(lenSeq,1);

micTimeData = zeros(K,N);           % Init variable for microphone data

for k = 1:K
    micTimeData(k,:) = filter(h(k,:),1,srcTimeData);
end
%%
micNewTimeData = micTimeData;       % Assignment of data
lenFrame = 501;                     % Length of frame
offsetSample = 4e3;                 % All microphones received a signal

numFrames = floor(lenSeq-offsetSample)/lenFrame;  % Number of frames to process
xGrid = -3:0.1:9;             
yGrid = -3:0.1:12;

micNewTimeData = micTimeData;       % Assignment of data
maskMax = ones(length(xGrid),length(yGrid));

numWallsAndSource = 5;
%% Main Loop
for l = 1 : numWallsAndSource
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
            srcFreqData,K,micPos,1,lenFrame,fs);
        tJ{frame} = J{1};
        Jtotal = Jtotal + J{1};
    end
    estVrtPos{l} = [estX estY];
    Jt{l} = Jtotal./max(Jtotal(:));
    [~,indxEst] = findMaximizer(Jt{l},maskMax);
    estPos(l,:) = [xGrid(indxEst(1)) yGrid(indxEst(2))];
    disp('Done Coarse Estimate.')
    disp('---------------------')
    %% Refined Estimate
    disp('----------------------------')
    disp('Starting Refinement.........')
    disp('----------------------------')
    xF = estPos(l,1)-0.5:0.005:estPos(l,1)+0.5;
    yF = estPos(l,2)-0.5:0.005:estPos(l,2)+0.5; 
    [estPos(l,:), e, xF, yF] = extensiveSearch(micNewTimeData,srcTimeData,micPos,...
        fs,xF,yF);
    disp('Done Refinement')
    %% Fixing Mask
    if l > 1
        phi = pi + atan2(estPos(l,2)-estPos(1,2),estPos(l,1)-estPos(1,1));
        if (phi < 5*pi/4 && phi > 3*pi/4)
            maskMax(xGrid > 1.5,:) = 0;
        elseif (phi > 5*pi/4 && phi < 7*pi/4)
            maskMax(:,yGrid > 0.7) = 0;
        elseif (phi > 7*pi/4 || phi < pi/4)
            maskMax(xGrid < 0.4,:) = 0;
        elseif (phi > pi/4 && phi < 3*pi/4)
            maskMax(:,yGrid < 0) = 0;
        end
    else
        maskMax(xGrid > 0.3 & xGrid < 8,yGrid > -0.1 & yGrid < 8) = 0;
    end
    %% Update Data
    delayData = nfGenDelayDataNonPeriodic(srcTimeData,micPos,[estPos(l,:) 1.26],fs);
    micNewTimeData = micNewTimeData - delayData';
    disp('Current Estimate')
    disp('----------------------------')
    disp(estPos(l,:))
    disp('----------------------------')
end

plotRoom(micPos, srcPos,vrtPos, estPos')