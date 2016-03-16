clear all, 
close all
% Simulation Parameters
radius = 0.06;
numOfMics = 3;
% Size of the room
roomDim = [4.5 6 2.5];

% numOfLoudSpeakers = 1:6 ;
numOfLoudSpeakers = 1:4;

K = numOfLoudSpeakers*numOfMics;
fs = 96e3;

% Reverberation time
T60 = [0.4 0.5 0.6];
% T60 = 0.4;

% Horizontal plane
zPos = 1.26;

% Number of Iterations
numIt = 50;
% Length of the Signal
lenSeq = round(T60*fs*3/2);
if (mod(lenSeq,2) == 1); lenSeq = lenSeq+1; end
% Processing size
lenFrame = 3001;
offSetSample = 4e3;
numFrames = floor((lenSeq-offSetSample)/lenFrame);
% Source to estimate
numOfSources = 5; % Four walls and source

% Simulation Dimension
% simDim = [2 3];
simDim = [3];

estMaxErr = nan([length(numOfLoudSpeakers), length(T60), ...
    length(simDim), numIt]);
estMeanErr = nan([length(numOfLoudSpeakers), length(T60), ...
    length(simDim), numIt]);
estPos = nan([numOfSources, 2, length(numOfLoudSpeakers), length(T60), ...
    length(simDim), numIt]);
%%
% Run Simulation
for nL = 1:length(numOfLoudSpeakers)
% for nL = 
    for t60 = 1:length(T60)
        for nD = 1:length(simDim)
            for it = 1:numIt
                disp(strcat('Running It: ',num2str(it)))
                disp('--------------------------------')
                try
                    [estPos(:,:,nL,t60,nD,it),...
                        estMaxErr(nL,t60,nD,it),...
                        estMeanErr(nL,t60,nD,it)] = roomEstSimulation(...
                        radius, numOfMics, numOfLoudSpeakers(nL),...
        K(nL), fs, roomDim, T60(t60), zPos, lenSeq(t60), lenFrame,...
        offSetSample, numFrames(nL),numOfSources, simDim(nD));
                catch
                    disp('Error')
                end
            end
        end
    end
end

% 12sensors2DT60:04 estMeanErr -> 0.0129 meanMax -> 0.0517 -median = 0

%%
mEstMeanErr = nanmean(estMeanErr,4)
mEstMaxErr =  nanmean(estMaxErr,4);

plot([0.4 0.5 0.6], mEstMaxErr)
