% Test of Sequential ML algorithm from a single virtual source
clear all, close all

P = 80;
fs = 20e3;
f0 = 5e3;
N = (P/f0 * fs) + 1;         % Hardcoded simulation parameters
w0 = 2*pi/N;
c = 344.8;

K = 6;
rho = 0.06;
p = w0*rho*fs/c;    % frequency-radius equivalence

% range2source
dSrc = 1.5; %[m]
etaSrc = w0*dSrc*fs/c;
% doa of the source
DOAsrc = pi/3;
% source position w.r.t origin
srcPos = dSrc * [cos(DOAsrc); sin(DOAsrc)];
betaSrc = 1/(4*pi*dSrc);

% infinite wall description
dist2wall = 2.5; %[m]
n = [-1; 0]; % normal to the wall
% position of the virtual source
rfxPos = srcPos + 2*dist2wall * n;
% doa and range of virtual source
DOArfx = atan2(rfxPos(2), rfxPos(1));
dRfx = norm(rfxPos);
etaRfx = w0*dRfx*fs/c;
betaRfx = 1/(4*pi*dRfx);

dist2wall = 0.5;
n = [1; 0]; % normal to the wall
% position of 2nd virtual source
rfxPos2 = srcPos + 2*dist2wall * n;
% doa and toa of virtual source
DOArfx2 = atan2(rfxPos2(2), rfxPos2(1));
dRfx2 = norm(rfxPos2);
etaRfx2 = w0*dRfx2*fs/c;
betaRfx2 = 1/(4*pi*dRfx2);


dist2wall = 3;
n = [0; -1]; % normal to the wall
% position of 2nd virtual source
rfxPos3 = srcPos + 2*dist2wall * n;
% doa and toa of virtual source
DOArfx3 = atan2(rfxPos3(2), rfxPos3(1));
dRfx3 = norm(rfxPos3);
etaRfx3 = w0*dRfx3*fs/c;
betaRfx3 = 1/(4*pi*dRfx3);


trueDOA = [DOAsrc; DOArfx; DOArfx2; DOArfx3];
trueTOA = [etaSrc; etaRfx; etaRfx2; etaRfx3]; 
trueBeta = [betaSrc; betaRfx; betaRfx2; betaRfx3];
L = length(trueTOA);

warning off

[micFreqData, srcFreqData, trueDOA, trueTOA, ...
    micTimeData, srcTimeData] = genTstMicData(K, p,L, trueDOA, trueTOA,...
    P,fs,f0,trueBeta);

[estDOA, estTOA, ~, J] = sequentialMLE_TOA_DOA(micTimeData,...
    srcTimeData, srcFreqData, K, p, L, N);

estBeta = updateBeta(micTimeData, srcTimeData, estDOA,...
                    estTOA, N, K, L, p);
%%
disp('Estimates')
disp(strcat('------DOAs---','--TOAs---', '-EstBetas---'))
estimates = [mod(rad2deg(estDOA),360), estTOA*c/(w0*fs), estBeta];
[~,or] = sort(estimates(:,1));
disp(estimates(or,:))
disp('True')
disp(strcat('------DOAs---','--TOAs---', '----Betas---'))
trueValues = [mod(rad2deg(trueDOA),360) trueTOA*c/(w0*fs) trueBeta];
[~,or] = sort(trueValues(:,1));
disp(trueValues(or,:))

p1 = polar(deg2rad(trueValues(:,1)), trueValues(:,2),'og');
hold on, p2 = polar(deg2rad(estimates(:,1)), estimates(:,2), 'xr');
legend([p1 p2],'True Values', 'Estimates')
title('Estimation of ISM')