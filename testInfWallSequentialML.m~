% Test of Sequential ML algorithm from a single virtual source
clear all, close all

P = 20;
fs = 48e3;

N = (P/1e3 * fs) + 1;         % Hardcoded simulation parameters
w0 = 2*pi/N;
c = 344.8;

K = 6;
rho = 0.06;

% range2source
dSrc = 1.5; %[m]
etaSrc = w0*dSrc*fs/c;
% doa of the source
DOAsrc = pi/2;
% source position w.r.t origin
srcPos = dSrc * [cos(DOAsrc); sin(DOAsrc)];
% infinite wall description
dist2wall = 1; %[m]
n = [-1; 0]; % normal to the wall
% position of the virtual source
rfxPos = srcPos + dist2wall * n;
% doa and range of virtual source
DOArfx = atan2(rfxPos(2), rfxPos(1));
dRfx = norm(rfxPos);
etaRfx = w0*dRfx*fs/c;

dist2wall = 2.5;
n = [1; 0]; % normal to the wall
% position of 2nd virtual source
rfxPos2 = srcPos + dist2wall + n;
% doa and toa of virtual source
DOArfx2 = atan2(rfxPos2(2), rfxPos2(1));
dRfx2 = norm(rfxPos2);
etaRfx2 = w0*dRfx2*fs/c;

trueDOA = [DOAsrc; DOArfx; DOArfx2];
trueTOA = [etaSrc; etaRfx; etaRfx2]; L = length(trueTOA);

[micFreqData, srcFreqData, trueDOA, trueTOA, ...
    micTimeData, srcTimeData] = genTstMicData(K, rho,L, trueDOA, trueTOA,...
    P,fs);

[estDOA, estTOA, ~] = sequentialMLE_TOA_DOA(micTimeData,...
    srcTimeData, srcFreqData, K, rho, L, N);

estBeta = updateBeta(micTimeData, srcTimeData, estDOA,...
                    estTOA, N, K, L, rho);
 
disp('Estimates')
disp(strcat('------DOAs---','--TOAs---', '-EstBetas---'))
disp([rad2deg(estDOA), estTOA*c/(w0*fs), estBeta])
disp('True')
disp(strcat('------DOAs---','--TOAs---'))
disp([rad2deg(trueDOA) trueTOA*c/(w0*fs)])