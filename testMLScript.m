% Test of functions for evaluation of cost function 
% Single sinusoid
clear all, close all

K = 6;
rho = 0.06;
% 
% % sinosoid frequency
% f = 1e3;
% w = 2*pi*f;
% % number of samples
% P = 5;
% fs = 10e3;
% samples = ceil(fs * P/f);
% 
% t = (0:samples)/fs;
% % original signal
% s = sin(w*t).'; N = length(s); l = floor(N/2);
% a = applyFFT(s,N);
% 
% % simulation parameters
% theta = pi/3; theta2 = 2*pi/3; theta3 = pi;
% eta = pi/3; eta2 = pi/2; eta3 = 4*pi/3;
% 
% % microphone data
% micData1 = genDelayData(s,theta,eta,K,rho);
% micData2 = genDelayData(s,theta2,eta2,K,rho);
% micData3 = genDelayData(s,theta3,eta3,K,rho);
% 
% micData = micData1 + micData2 + micData3;
% 
% micFreqData = zeros(size(micData));
% for k = 1:K
%     micFreqData(:,k) = N*applyFFT(micData(:,k),N); 
% end

[micFreqData, srcData, trueDOA, trueTOA] = genTstMicData(K, rho);

[estDOA, estTOA, estBeta,J] = estML_TOA_DOA(micFreqData, srcData, K, rho);

imagesc(J), xlabel('Range [Eta]'), ylabel('DOA [Theta]')

disp('Estimate'); 
disp(strcat('DOA [Theta]: ',num2str(estDOA)))
disp(strcat('Range [Eta]: ',num2str(estTOA)))
disp(strcat('Amplitude [Beta]: ', num2str(estBeta)))

