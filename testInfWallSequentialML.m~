% Test of Sequential ML algorithm from a single virtual source
clear all, close all

P = 150;
fs = 48e3;
f0 = 10e3;
N = (P/f0 * fs) + 1;         % Hardcoded simulation parameters
w0 = 2*pi/N;
c = 344.8;

K = 6;
rho = 0.12;

% range2source
dSrc = 1.5; %[m]
etaSrc = w0*dSrc*fs/c;
% doa of the source
DOAsrc = pi/3;
% source position w.r.t origin
srcPos = dSrc * [cos(DOAsrc); sin(DOAsrc)];
% infinite wall description
dist2wall = 2.5; %[m]
n = [-1; 0]; % normal to the wall
% position of the virtual source
rfxPos = srcPos + 2*dist2wall * n;
% doa and range of virtual source
DOArfx = atan2(rfxPos(2), rfxPos(1));
dRfx = norm(rfxPos);
etaRfx = w0*dRfx*fs/c;

dist2wall = 0.5;
n = [1; 0]; % normal to the wall
% position of 2nd virtual source
rfxPos2 = srcPos + 2*dist2wall * n;
% doa and toa of virtual source
DOArfx2 = atan2(rfxPos2(2), rfxPos2(1));
dRfx2 = norm(rfxPos2);
etaRfx2 = w0*dRfx2*fs/c;


dist2wall = 3;
n = [0; -1]; % normal to the wall
% position of 2nd virtual source
rfxPos3 = srcPos + 2*dist2wall * n;
% doa and toa of virtual source
DOArfx3 = atan2(rfxPos3(2), rfxPos3(1));
dRfx3 = norm(rfxPos3);
etaRfx3 = w0*dRfx3*fs/c;


% trueDOA = [DOAsrc; DOArfx; DOArfx2; DOArfx3];
% trueTOA = [etaSrc; etaRfx; etaRfx2; etaRfx3]; 
trueDOA = [DOAsrc; DOArfx];
trueTOA = [etaSrc; etaRfx];
L = length(trueTOA);
%%
[micFreqData, srcFreqData, trueDOA, trueTOA, ...
    micTimeData, srcTimeData] = genTstMicData(K, rho,L, trueDOA, trueTOA,...
    P,fs);

[estDOA, estTOA, estBeta] = sequentialMLE_TOA_DOA(micTimeData,...
    srcTimeData, srcFreqData, K, rho, L, N);

estBeta = updateBeta(micTimeData, srcTimeData, estDOA,...
                    estTOA, N, K, L, rho);
%%
disp('Estimates')
disp(strcat('------DOAs---','--TOAs---', '-EstBetas---'))
disp([rad2deg(estDOA), estTOA*c/(w0*fs), estBeta])
disp('True')
disp(strcat('------DOAs---','--TOAs---'))
disp([rad2deg(trueDOA) trueTOA*c/(w0*fs)])