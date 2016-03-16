function [estPos, estMaxErr, estMeanErr] = roomEstSimulation(radius, numOfMics, numOfLoudSpeakers,...
    K, fs, roomDim, T60, zPos, lenSeq, lenFrame,offSetSample, numFrames,...
    numOfSources, simDim)
%% Creating Data 
    % Position of the arrays
    xBase = 0.3 + (0:(numOfLoudSpeakers-1))/numOfLoudSpeakers * roomDim(1);
    yBase = 2 + 2.5*rand(1,numOfLoudSpeakers);
    xCenterPos =  min(xBase + 0.1*randn(1,numOfLoudSpeakers), roomDim(1) - 0.3);
    xCenterPos = max(xCenterPos,0.3);
    yCenterPos = min(yBase + 0.3*randn(1,numOfLoudSpeakers), roomDim(2) - 0.3);
    yCenterPos = max(yCenterPos, 0.3);

    % Source Position
    xSrcPos = (roomDim(1)- 0.5) * rand(1) + 0.5;
    ySrcPos = (2 - 1.5) * rand(1) + 1;
    srcPos = [xSrcPos ySrcPos zPos];
    vrtSrcPos = genSrcsFromWalls(srcPos,roomDim(1:2));
    
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
        L{l} = LoudSpeaker(centerPos, numOfMics, radius, simDim);
        L{l} = L{l}.genNewSim(srcSignal, srcPos, roomDim, T60,fs);
        indx = (1:numOfMics) + (l-1)*numOfMics;
        micTemporalData(:,indx) = L{l}.micSignals';
        micPos(indx,:) = L{l}.micPos(:,1:2);
    end
    micNewTimeData = micTemporalData;
    %% Finding the reflections
    Jsrc = cell(numOfSources,1);
    estPos = nan(numOfSources,2);
    
    % Speed ups the simulation (Source location given a priori)
    estPos(1,:) = srcPos(1:2);

    for l = 1:numOfSources
        if l > 1 % Assumes source location known
            %% Coarse Estimate
            Jtotal = 0;
            for frame = 1:numFrames
                startSample =  offSetSample + (frame-1)*lenFrame + 1;
                endSample = startSample + lenFrame - 1;
                micData = micNewTimeData(startSample:endSample,:);
                srcData = srcSignal(startSample:endSample);
                srcFreqData = getFreqMicData(srcData,lenFrame,1);
                micFreqData = getFreqMicData(micData,lenFrame,K);
                if l > 10
                    J = nfEstML(micFreqData,...
                        srcFreqData,micPos,fs, xPoints, yPoints,w);
                else
                    J = nfEstML(micFreqData,...
                        srcFreqData,micPos,fs, xPoints, yPoints);
                end
                Jtotal = Jtotal + J;
            end
            Jsrc{l} = Jtotal./max(Jtotal(:));
            [~,indxMax] = max(Jsrc{l});
            estPos(l,:) = [xPoints(indxMax) yPoints(indxMax)];
            %% Refined Estimate
            if l > 1
                sector = findSectorOfPoint(estPos(1,:), estPos(l,:));
                if (sector == 1 || sector == 3)     % Grid in x-axis
%                     xF = estPos(l,1) - 0.8 : 0.005 : estPos(l,1) + 0.8;
                    xF = estPos(l,1) - 2 : 0.005 : estPos(l,1) + 2;
                    yF = estPos(1,2) - 0.01 : 0.005 : estPos(1,2) + 0.01;
                elseif (sector == 2 || sector == 4) % Grid in y-axis
%                     yF = estPos(l,2) - 0.8 : 0.005 : estPos(l,2) + 0.8;
                    yF = estPos(l,2) - 2 : 0.005 : estPos(l,2) + 2;
                    xF = estPos(1,1) - 0.01 : 0.005 : estPos(1,1) + 0.01;
                end
            else
                xF = estPos(l,1) - 0.5 : 0.005 : estPos(l,1) + 0.5;
                yF = estPos(l,2) - 0.5 : 0.005 : estPos(l,2) + 0.5;
            end
            [XM, YM] = meshgrid(xF,yF);
            xP = XM(:); yP = (YM(:));
%             tmpMicPos = [micPos zPos*ones(K,1)];
            if l > 2 % Includes knowledge of walls (After the 1st reflx is known)
%                 estPos(l,:) = extensiveSearchWalls(micNewTimeData', srcSignal,...
%                     tmpMicPos, fs, xF, yF,w);
                J = nfEstML(micFreqData,...
                    srcFreqData,micPos,fs, xP, yP,w);
%             elseif l == 2 % Removes any biased in the cost function
%                 estPos(l,:) = extensiveSearch(micNewTimeData', srcSignal,...
%                     tmpMicPos, fs, xF, yF);
            else          % Original Cost function 
                J = nfEstML(micFreqData,...
                    srcFreqData,micPos,fs, xP, yP);
%                 estPos(l,:) = extensiveSearch(micNewTimeData', srcSignal,...
%                     tmpMicPos, fs, xF, yF,1);
            end
           
            
            [~,indxMax] = max(J);
            estPos(l,:) = [xP(indxMax) yP(indxMax)];
        end
        %% Reducing Feasible set
        if l > 1
            switch (sector)
                case 1 %right wall
                    mask = ( mask .* (xGrid < mXUlim) ) ~= 0;
                    w(l-1,:) = [1 0.5*(estPos(l,1) + estPos(1,1))];
                case 2 % front wall
                    mask = ( mask .* (yGrid < mYUlim) ) ~= 0;
                    w(l-1,:) = [0 0.5*(estPos(l,2) + estPos(1,2))];
                case 3 % left wall
                    mask = ( mask .* (xGrid > mXLlim) ) ~= 0;
                    w(l-1,:) = [1 0.5*(estPos(l,1) + estPos(1,1))];
                case 4 % rear wall
                    mask = ( mask .* (yGrid > mYLlim) ) ~= 0;
                    w(l-1,:) = [0 0.5*(estPos(l,2) + estPos(1,2))];
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
    %% Computing Error
    or = [3 4 1 2]; % Reflection order in Virtual sources
    
    sortSrcs = nan(4,2);
    for i = 2:numOfSources
        sector = findSectorOfPoint(estPos(1,:),estPos(i,:));
        sortSrcs(or == sector,:) = estPos(i,:);
    end
    e = norms(vrtSrcPos' - sortSrcs,2,2);
    estMaxErr = max(e);
    estMeanErr = mean(e);
end