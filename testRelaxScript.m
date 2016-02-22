% Test of RELAX algorithm estimating L sources, with L known a-priori
% [It has issues with convergence by errors in estimation]
clear all, close all

K = 6;
rho = 0.06;
L = 4;

% Creates K virtual signals coming from L random waves impinging in a UCA
% with radius rho.
trueDOA = 2*pi/4 * (0:(L-1))';
trueTOA = 2*pi/4 * (0:(L-1))';
[micFreqData, srcFreqData, trueDOA, trueTOA, ...
    micTimeData, srcTimeData] = genTstMicData(K, rho,L, trueDOA, trueTOA);

N = length(srcTimeData);

% RELAX : BEGIN
%%
[estDOA, estTOA, estBeta] = estRELAX_TOA_DOA(micTimeData, srcTimeData, ...
    srcFreqData, K, rho, L, N);
%%
disp('Estimates')
disp(strcat('------DOAs---','---TOAs---', '---Betas---'))
disp([estDOA, estTOA, estBeta])