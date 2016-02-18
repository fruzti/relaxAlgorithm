% Test of functions for evaluation of cost function 
clear all, close all

K = 6;
rho = 0.06;
L = 3;

% Creates K virtual signals coming from L random waves impinging in a UCA
% with radius rho.
[micFreqData, srcData, trueDOA, trueTOA,micTimeData,srcTimeData] =...
    genTstMicData(K, rho,L);

% Computes the MLE for the DOA and TOA based on NLS
[estDOA, estTOA, estBeta, J] = estML_TOA_DOA(micFreqData, srcData, K, rho);

imagesc(J), xlabel('Range [Eta]'), ylabel('DOA [Theta]')

% Estimates results
disp('Estimate'); 
disp(strcat('DOA [Theta]: ',num2str(deg2rad(estDOA))))
disp(strcat('Range [Eta]: ',num2str(deg2rad(estTOA))))
disp(strcat('Amplitude [Beta]: ', num2str(estBeta)))

