function [estPos, estMaxErr, estMeanErr] = roomEstSimulationSNR(K, fs,...
    lenFrame,offSetSample, numFrames, srcPos,roomDim,...
    micPos, micNewTimeData, srcSignal, numOfSources, nPow)

    % Original Grid
    xLlim = -(srcPos(1) + 1); xUlim = 2*(roomDim(1) - srcPos(1)) + 1 + srcPos(1);
    yLlim = -(srcPos(2) + 1); yUlim = 2*(roomDim(2) - srcPos(2)) + 1 + srcPos(2);
    gridInX = xLlim : 0.1 : xUlim;
    gridInY = yLlim : 0.1 : yUlim;
    [xGrid, yGrid] = meshgrid(gridInX, gridInY);
    xPoints = xGrid(:); yPoints = yGrid(:);
    
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
                micData = micData + randn(size(micData))*sqrt(nPow);
                srcData = srcSignal(startSample:endSample);
                srcFreqData = getFreqMicData(srcData,lenFrame,1);
                micFreqData = getFreqMicData(micData,lenFrame,K);
                if l > 2
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
                switch (sector)
                    case 1
                        xF = (estPos(l,1) - 1) : 0.005 : ...
                            min(estPos(l,1) + 2,xUlim-1);
                        yF = estPos(1,2) - 0.03 : 0.005 : estPos(1,2) + 0.03;
                    case 2
                        yF = max((estPos(l,2) - 1),maxY) : 0.005 : ...
                            min(estPos(l,2) + 3.3,yUlim-1);
                        xF = estPos(1,1) - 0.01 : 0.005 : estPos(1,1) + 0.01;
                    case 3
                        xF = max(estPos(l,1) - 2,xLlim) : 0.005 : ...
                            min(estPos(l,1) + 2,minX);
                        yF = estPos(1,2) - 0.01 : 0.005 : estPos(1,2) + 0.01;
                    case 4
                        yF = max(estPos(l,2)  - 3,yLlim) : 0.005 : ...
                            min(estPos(l,2) + 1,minY);
                        xF = estPos(1,1) - 0.01 : 0.005 : estPos(1,1) + 0.01;
                end
%                 if (sector == 1 || sector == 3)     % Grid in x-axis
% %                     xF = estPos(l,1) - 0.8 : 0.005 : estPos(l,1) + 0.8;
%                     xF = max(estPos(l,1) - 2,xLlim) : 0.005 : min(estPos(l,1) + 2,xUlim);
%                     yF = estPos(1,2) - 0.01 : 0.005 : estPos(1,2) + 0.01;
%                 elseif (sector == 2 || sector == 4) % Grid in y-axis
% %                     yF = estPos(l,2) - 0.8 : 0.005 : estPos(l,2) + 0.8;
%                     yF = estPos(l,2) - 1 : 0.005 : estPos(l,2) + 3.3;
%                     xF = estPos(1,1) - 0.01 : 0.005 : estPos(1,1) + 0.01;
%                 end
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
            sector = findSectorOfPoint(estPos(1,:), estPos(l,:));
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
            maxX = 2*max(max(micPos(:,1)) - estPos(1,1),0) + mXUlim; 
            maxY = 2*max(max(micPos(:,2)) - estPos(1,2),0) + mYUlim; 
            minX = 2*min(min(micPos(:,1)) - estPos(1,1),0) + mXLlim; 
            minY = 2*min(min(micPos(:,2)) - estPos(1,2),0) + mYLlim; 
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
    vrtSrcPos = genSrcsFromWalls(srcPos,roomDim(1:2));
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