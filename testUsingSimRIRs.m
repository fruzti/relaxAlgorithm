% Test of Sequential ML using room acoustics
run experimentalSetup.m

% Simulation Parameters
zPos = 1.26;
srcPos = [ladLoud' zPos];
srcPos(:,1) = srcPos(:,1) + 4;
micPos = [arrayElems; ones(1,size(arrayElems,2))*zPos]';
micPos(:,1) = micPos(:,1) + 4;

arrayPos2 = mean(micPos(7:12,:));
arrayPos1 = mean(micPos(1:6,:));
arrayPos3 = mean(micPos(13:18,:));

clearvars -except arrayPos1 arrayPos2 arrayPos3 srcPos

N = 201; 
w0 = 2*pi/N; 
fs =  10e3; c = 344.8;

trueDOA1 = atan2(srcPos(2) - arrayPos1(2), srcPos(1) - arrayPos1(1));
trueDOA2 = atan2(srcPos(2) - arrayPos2(2), srcPos(1) - arrayPos2(1));
trueDOA3 = atan2(srcPos(2) - arrayPos3(2), srcPos(1) - arrayPos3(1));

d1 = norm(srcPos - arrayPos1); 
d2 = norm(srcPos - arrayPos2);
d3 = norm(srcPos - arrayPos3);

trueTOA1 = w0*d1*fs/c;
trueTOA2 = w0*d2*fs/c;
trueTOA3 = w0*d3*fs/c;

if N*c/fs < d2
    error('Ambiguity in Range : Folded range')
end

% srcTimeData = mls(16,1);
N1 = 1e3;
srcTimeData = randn(N1,1);
% srcTimeData = mls(7,6); 
N1 = 3*128;
% [srcTimeData, fT] = audioread('voice.wav');
% srcTimeData = resample(srcTimeData(:,1), fs, fT);

% Generate Microphones Signals
fileNameRIRs = 'ISM_RIRs_MLE.mat'; % T60 = 0.4 // fs = 20e3
fs = 20e3;
AuData = ISM_AudioData(fileNameRIRs, srcTimeData);
srcTimeData = srcTimeData(N1:end);
micTimeData = AuData(N1:end,7:12); N = size(micTimeData,1);
w0 = 2*pi/N;
% offset = 100;
% srcTimeData = srcTimeData(offset:N + offset - 1);
% srcFreqData = applyFFT(srcTimeData,N);

% array's info
K = 6;
rho = 0.06;
p = w0*rho*fs/c;    % frequency-radius equivalence

% Only one array (the second - almost in the center)
% micTimeData = AuData(offset:N+offset-1,7:12);
micFreqData = getFreqMicData(micTimeData, N, K);
srcFreqData = applyFFT(srcTimeData, N);

% number of reflections
L = 4; % three main walls and source
%%
% [estDOA, estTOA, ~, J] = sequentialMLE_TOA_DOA(micTimeData,...
%     srcTimeData, srcFreqData, K, p, L, N);
[estDOA, estTOA, ~, J] = sequentialMLE_TOA_DOA(micTimeData,...
    srcTimeData, srcFreqData, K, p, L, N);

estBeta = updateBeta(micTimeData, srcTimeData, estDOA,...
                    estTOA, N, K, L, p);

disp('Estimates')
disp(strcat('------DOAs---','--TOAs---', '-EstBetas---'))
estimates = [mod(rad2deg(estDOA),360), estTOA*c/(w0*fs)];
[~,or] = sort(estimates(:,1));
disp(estimates(or,:))
disp('True')
disp(strcat('------DOAs---','--TOAs---', '-EstBetas---'))
disp([mod(rad2deg(trueDOA2),360) d2])
vrtPos = genSrcsFromWalls(srcPos,[4 6]);
[tVrtPos, rVrtPos] = cart2pol(vrtPos(1,:)-arrayPos2(1), vrtPos(2,:)-arrayPos2(2));
p2 = polar([trueDOA2;tVrtPos'], [d2; rVrtPos'], 'og');
hold on, p1 = polar(deg2rad(estimates(:,1)), estimates(:,2), 'xr');
legend([p1 p2], 'Estimates', 'True')

