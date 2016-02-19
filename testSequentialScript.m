% Test of Sequential ML algorithm estimating L sources, with L known a-priori
clear all, close all

K = 6;
rho = 0.06;
L = 3;

% Creates K virtual signals coming from L random waves impinging in a UCA
% with radius rho.
trueDOA = 2*pi/4 * (0:(L-1))';
trueTOA = 2*pi/4 * (0:(L-1))';
[micFreqData, srcFreqData, trueDOA, trueTOA, ...
    micTimeData, srcTimeData] = genTstMicData(K, rho,L, trueDOA, trueTOA);

N = length(srcTimeData);

% RELAX : BEGIN
[estDOA, estTOA, estBeta] = sequentialMLE_TOA_DOA(micTimeData,...
    srcTimeData, srcFreqData, K, rho, L, N);

% Gain Correction
S = zeros(N*K,L);
for l = 1:L
    tmp = genDelayData(srcTimeData,estDOA(l), estTOA(l), K, rho);
    S(:,l) = tmp(:);
end

hatBeta = S\micTimeData(:);

disp('Estimates')
disp(strcat('------DOAs---','--TOAs---', '-EstBetas-','-FinalBetas---'))
disp([estDOA, estTOA, estBeta hatBeta])