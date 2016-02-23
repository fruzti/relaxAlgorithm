% Test of Sequential ML using room acoustics
run experimentalSetup.m

% Simulation Parameters
zPos = 1.26;
srcPos = [ladLoud' zPos];
srcPos(:,1) = srcPos(:,1) + 4;
micPos = [arrayElems; ones(1,size(arrayElems,2))*zPos]';
micPos(:,1) = micPos(:,1) + 4;

arrayPos = mean(micPos(7:12,:));

clearvars -except arrayPos srcPos

N = 201; 
w0 = 2*pi/N; 
fs =  10e3; c = 344.8;

trueDOA = atan2(srcPos(2) - arrayPos(2), srcPos(1) - arrayPos(1));
d = norm(srcPos - arrayPos);
trueTOA = w0*d*fs/c;

if N*c/fs < d
    error('Ambiguity in Range : Folded range')
end

% srcTimeData = mls(16,1);
N1 = 1e3;
srcTimeData = randn(N1,1);
% [srcTimeData, fT] = audioread('voice.wav');
% srcTimeData = resample(srcTimeData(:,1), fs, fT);

% Generate Microphones Signals
fileNameRIRs = 'ISM_RIRs.mat'; % T60 = 0.6
AuData = ISM_AudioData(fileNameRIRs, srcTimeData);

offset = 100;
srcTimeData = srcTimeData(offset:N + offset - 1);
srcFreqData = applyFFT(srcTimeData,N);

% array's info
K = 6;
rho = 0.06;
p = w0*rho*fs/c;    % frequency-radius equivalence

% Only one array (the second - almost in the center)
micTimeData = AuData(offset:N+offset-1,7:12);

% number of relfections
L = 4; % three main walls and source
%%
% [estDOA, estTOA, ~, J] = sequentialMLE_TOA_DOA(micTimeData,...
%     srcTimeData, srcFreqData, K, p, L, N);
[estDOA, estTOA, ~, J] = sequentialMLE_TOA_DOA(micTimeData,...
    srcTimeData, srcFreqData, K, p, L, N);

estBeta = updateBeta(micTimeData, srcTimeData, estDOA,...
                    estTOA, N, K, L, p);
%%
disp('Estimates')
disp(strcat('------DOAs---','--TOAs---', '-EstBetas---'))
estimates = [mod(rad2deg(estDOA),360), estTOA*c/(w0*fs)];
[~,or] = sort(estimates(:,1));
disp(estimates(or,:))
disp('True')
disp(strcat('------DOAs---','--TOAs---', '-EstBetas---'))
disp([mod(rad2deg(trueDOA),360) trueTOA*c/(w0*fs)])

p1 = polar(deg2rad(estimates(:,1)), estimates(:,2), 'xr');
hold on, p2 = polar(trueDOA, d, 'og');
legend([p1 p2], 'Estimates', 'True')

