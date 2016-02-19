% Test of functions for evaluation of cost function 
clear all, close all

K = 6;
rho = 0.06;
L = 4;

% Creates K virtual signals coming from L random waves impinging in a UCA
% with radius rho.
trueDOA = 2*pi/4 * (0:(L-1))';
trueTOA = 2*pi/4 * (0:(L-1))';
[micFreqData, srcFreqData, trueDOA, trueTOA,micTimeData,srcTimeData] =...
    genTstMicData(K, rho,L, trueDOA, trueTOA);

% Computes the MLE for the DOA and TOA based on NLS
[estDOA, estTOA, estBeta, J] = estML_TOA_DOA(micFreqData, srcFreqData, K, rho);

imagesc(J), xlabel('Range [Eta]'), ylabel('DOA [Theta]')

% Estimates results
disp('Estimate'); 
disp(strcat('DOA [Theta]: ',num2str(estDOA)))
disp(strcat('Range [Eta]: ',num2str(estTOA)))
disp(strcat('Amplitude [Beta]: ', num2str(estBeta)))

