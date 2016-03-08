run experimentalSetup.m

% Simulation Parameters
zPos = 1.26;
srcPos = [ladLoud' zPos];
srcPos(:,1) = srcPos(:,1) + 4;
micPos = [arrayElems; ones(1,size(arrayElems,2))*zPos]';
micPos(:,1) = micPos(:,1) + 4;
refMics = [singleRef pairRef]';
refMics(:,1) = refMics(:,1) + 4;

arrayPos1 = mean(micPos(1:6,:));    arrayPos{1} = arrayPos1;
arrayPos2 = mean(micPos(7:12,:));   arrayPos{2} = arrayPos2;
arrayPos3 = mean(micPos(13:18,:));  arrayPos{3} = arrayPos3;

clearvars -except arrayPos srcPos micPos refMics
%% Generating RIRs
c = 344.8;                  % Speed of sound [m/s]
fs = 96e3;                  % Sampling frequency [Hz]

L = [4.5 6 2.5];            % Room dimensions [x y z]
beta = 0.4;                 % Reverberation time [s]
nsample = floor(beta*fs);             % Number of samples
order = -1;                 % -1 -> maximum reflection order
mtype = 'omnidirectional';  % Type of microphone
dim = 3;                    % Dimension of the room
orientation = [0 0];        % Microphone orientation [azimuth elevation]
hp_filter = 1;              % Enable high-pass filter

setupFilter.L = L;
setupFilter.beta = beta;
setupFilter.nsample = nsample;
setupFilter.order = order;
setupFilter.mtype = mtype;
setupFilter.dim = dim;
setupFilter.hp_filter = hp_filter;
setupFilter.c = c;
setupFilter.fs = fs;
setupFilter.orientation = orientation;

K = size(micPos,1);
for k = 1:K
    h(k,:) = rir_generator(c, fs, micPos(k,:), srcPos, L, beta, nsample,...
        mtype, order, dim,orientation, hp_filter);
    disp(k)
end

%% Generating the microphone data
% order = 7;                       % Order of MLS
% P = 20;                          % Periods
% srcTimeData = mls(order,P); 
N1 = 20e3;
srcTimeData = randn( N1 ,1);
% srcTimeData = repmat(srcTimeData,40,1);
% triaWnd = triang(1002); win = triaWnd(502:end);
% windT = repmat(win,ceil(N1/501),1);
% srcTimeData = srcTimeData.*windT(1:N1);
% srcTimeData = mls(nextpow2(N1),1);
N = length(srcTimeData);
micTimeData = zeros(K,N);
for k = 1:K
    micTimeData(k,:) = filter(h(k,:),1, srcTimeData);
end
lenSeq = 500+1;
offsetSample = 4e3;
clear estX estY
%%
numSec = 30;
for seq = 1:numSec
    startSample = offsetSample + (seq-1)*lenSeq;
    micTimeDataT = micTimeData(:, startSample + 1 : startSample + lenSeq)';
    srcTimeDataT = srcTimeData(startSample + 1 : startSample + lenSeq);
    srcFreqData = getFreqMicData(srcTimeDataT,lenSeq,1);
    [estX(seq),estY(seq),~,J] = nfSequentialMLE_TOA_DOA(micTimeDataT,srcTimeDataT,...
        srcFreqData,K,micPos,1,lenSeq,fs);
    tJ{seq} = J{1};
end
%%
jT = ones(size(tJ{1}));
for seq = 1:numSec
    jT = jT + tJ{seq};
%     Jp{seq} = nfEvalCost(micFreqData,srcFreqData, micPos',...
%                 srcPos,l,N,fs);
    tmpM = tJ{seq};
    xjT(:,seq) = sum(tmpM(:,34:36),2);
    xjT(:,seq) = xjT(:,seq)/max(xjT(:,seq));
    
    yjT(:,seq) = sum(tmpM(39:42,:));
    yjT(:,seq) = yjT(:,seq)/max(yjT(:,seq));
end
jT  = jT'/max(jT(:));
figure, imagesc(xGrid,yGrid,jT), set(gca,'YDir','normal')
%%
plot(xGrid,sum(xjT,2))
figure,plot(yGrid,sum(yjT,2))
figure, plot(xGrid,sum(jT(32:36,:)))
hold on, plot(micPos(:,1),2*ones(18,1),'k*')
figure,
plot(yGrid,(sum(jT(:,39:42),2)))
hold on, plot(micPos(:,2),2*ones(18,1),'k*')
%% Testing Sequential ML Estimator
warning off
[estX,estY,~,J] = nfSequentialMLE_TOA_DOA(micTimeData,srcTimeData,...
    srcFreqData,K,micPos,1,N,fs);
xGrid = -3:0.1:12;             
yGrid = -3:0.1:12;
%%
figure, imagesc(xGrid,yGrid,J{1}'), set(gca,'YDir','normal')
title('Cost Function for Source Location')
xlabel('X-Axis'), ylabel('Y-Axis')
disp('Estimate')
estimates = [estX estY];
disp(estimates)
disp('True')
disp(srcPos(1:2))

%%
[estX, estY, estBeta, J] = nfRelaxMLE_TOA_DOA(micTimeData',...
    srcTimeData, srcFreqData, K, micPos,L,N,fs,setupFilter);