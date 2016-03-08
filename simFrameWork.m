% Simulation Parameters
radius = 0.06;
numOfMics = 3;
numOfLoudSpeakers = 4;
K = numOfLoudSpeakers*numOfMics;
fs = 96e3;
% Size of the room
roomDim = [4.5 6 2.5];
% Reverberation time
T60 = 0.4;
% Horizontal plane
zPos = 1.26;
% Number of Iterations
numIt = 1;
% Length of the Signal
lenSeq = round(T60*fs/2);
if (mod(lenSeq,2) == 1); lenSeq = lenSeq+1; end
% Processing size
lenFrame = 501;
offSetSample = 4e3;
numFrames = floor((lenSeq-offSetSample)/lenFrame);
% Source to estimate
numOfSources = 5; % Four walls and source
%% Iterative Solution
for it = 1:numIt
    %% Creating Data 
    % Position of the arrays
    xCenterPos =  min([0.5 1.8 3.4 4.3] + 0.1*randn(1,4), roomDim(1) - 0.3);
    xCenterPos = max(xCenterPos,0.3);
    yCenterPos = min([3.5 4.3 3.4 2.3] + 0.3*randn(1,4), roomDim(2) - 0.3);
    yCenterPos = max(yCenterPos, 0.3);

    % Source Position
    xSrcPos = (roomDim(1)- 0.5) * rand(1) + 0.5;
    ySrcPos = (2 - 1.5) * rand(1) + 1;
    srcPos = [xSrcPos ySrcPos zPos];
    % Original Grid
    xLlim = -(srcPos(1) + 1); xUlim = 2*(roomDim(1) - srcPos(1)) + 1 + srcPos(1);
    yLlim = -(srcPos(2) + 1); yUlim = 2*(roomDim(2) - srcPos(2)) + 1 + srcPos(2);
    gridInX = xLlim : 0.1 : xUlim;
    gridInY = yLlim : 0.1 : yUlim;
    [xGrid, yGrid] = meshgrid(gridInX, gridInY);
    xPoints = xGrid(:); yPoints = yGrid(:);
    
    % Source Signal
    srcSignal = randn(lenSeq,1);

    % Creating arrays
    L = cell(numOfLoudSpeakers,1);
    micTemporalData = nan(lenSeq,K);
    micPos = nan(K,2); arrayPos = nan(K,2);
    
    %% Creating the simulated data
    for l = 1:numOfLoudSpeakers
        centerPos = [xCenterPos(l) yCenterPos(l) zPos];
        arrayPos(l,:) = centerPos(1:2);
        L{l} = LoudSpeaker(centerPos, numOfMics, radius);
        L{l} = L{l}.genNewSim(srcSignal, srcPos, roomDim, T60,fs);
        indx = (1:numOfMics) + (l-1)*numOfMics;
        micTemporalData(:,indx) = L{l}.micSignals';
        micPos(indx,:) = L{l}.micPos(:,1:2);
    end
    micNewTimeData = micTemporalData;
    vrtSrcPos = genSrcsFromWalls(srcPos(1:2),roomDim(1:2));
    plotRoom(micPos,srcPos(1:2),vrtSrcPos,numOfMics);
    %% Finding the reflections
    Jsrc = cell(numOfSources,1);
    estPos = nan(numOfSources,2);
    
    for l = 1:numOfSources
        %% Coarse Estimate
        Jtotal = 0;
        disp(strcat('Estimation Source #',num2str(l)))
        disp('Init Coarse Estimate')
        disp('---------------------')
        for frame = 1:numFrames
            startSample =  offSetSample + (frame-1)*lenFrame + 1;
            endSample = startSample + lenFrame - 1;
            micData = micNewTimeData(startSample:endSample,:);
            srcData = srcSignal(startSample:endSample);
            srcFreqData = getFreqMicData(srcData,lenFrame,1);
            micFreqData = getFreqMicData(micData,lenFrame,K);
            J = nfEstML(micFreqData,...
                srcFreqData,micPos,fs, xPoints, yPoints);
            Jtotal = Jtotal + J;
        end
        Jsrc{l} = Jtotal./max(Jtotal(:));
        [~,indxMax] = max(Jsrc{l});
        estPos(l,:) = [xPoints(indxMax) yPoints(indxMax)];
        disp(estPos(l,:))
        disp('---------------------')
        %% Refined Estimate
        disp('Init Refined Estimate')
        disp('---------------------')
        if l > 1
            sector = findSectorOfPoint(estPos(1,:), estPos(l,:));
            if (sector == 1 || sector == 3) % Grid in x-axis
                xF = estPos(l,1) - 0.5 : 0.005 : estPos(l,1) + 0.5;
                yF = estPos(1,2) - 0.01 : 0.005 : estPos(1,2) + 0.01;
            elseif (sector == 2 || sector == 4) % Grid in y-axis
                yF = estPos(l,2) - 0.5 : 0.005 : estPos(l,2) + 0.8;
                xF = estPos(1,1) - 0.01 : 0.005 : estPos(1,1) + 0.01;
            end
        else
            xF = estPos(l,1) - 0.5 : 0.005 : estPos(l,1) + 0.5;
            yF = estPos(l,2) - 0.5 : 0.005 : estPos(l,2) + 0.5;
        end
        tmpMicPos = [micPos zPos*ones(K,1)];
        estPos(l,:) = extensiveSearch(micNewTimeData', srcSignal, tmpMicPos, fs, xF, yF);
        disp(estPos(l,:))
        disp('---------------------')
        %% Reducing Feasible set
        if l > 1
            switch (sector)
                case 1
                    mask = ( mask .* (xGrid < mXUlim) ) ~= 0;
                case 2
                    mask = ( mask .* (yGrid < mYUlim) ) ~= 0;
                case 3
                    mask = ( mask .* (xGrid > mXLlim) ) ~= 0;
                case 4
                    mask = ( mask .* (yGrid > mYLlim) ) ~= 0;
            end
            xPoints = xGrid(mask);
            yPoints = yGrid(mask);
        else
            mXLlim = estPos(1,1) - 0.3; mXUlim = estPos(1,1) + 0.3;
            mYLlim = estPos(1,2) - 0.3; mYUlim = estPos(1,2) + 0.3;
            maxX = 2*max(max(micPos(1,:)) - estPos(1,1),0) + mXUlim; 
            maxY = 2*max(max(micPos(2,:)) - estPos(1,2),0) + mYUlim; 
            minX = 2*min(min(micPos(1,:)) - estPos(1,1),0) + mXLlim; 
            minY = 2*min(min(micPos(2,:)) - estPos(1,2),0) + mYLlim; 
            nXgrid = xGrid .* (xGrid >= mXLlim) .* ...
                (xGrid <= mXUlim) .* ~((yGrid < maxY) & (yGrid >= minY) );
            nYgrid = yGrid .* (yGrid >= mYLlim) .* ...
                (yGrid <= mYUlim) .* ~((xGrid < maxX) & (xGrid >= minX) );
            mask = (nXgrid + nYgrid) ~= 0;
            xPoints = xGrid(mask);
            yPoints = yGrid(mask);
        end
        %% Update Data
        delayData = nfGenDelayDataNonPeriodic(srcSignal, micPos,...
            estPos(l,:), fs);
        micNewTimeData = micNewTimeData - delayData;
    end
    
end
vrtSrcPos = genSrcsFromWalls(srcPos,roomDim(1:2));
plotRoom(micPos, srcPos,vrtSrcPos,numOfMics, estPos')