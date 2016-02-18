% Test of functions for evaluation of cost function 
% Single sinusoid
clear all, close all
% sinosoid frequency
f = 1e3;
w = 2*pi*f;
% number of samples
P = 5;
fs = 10e3;
samples = ceil(fs * P/f);

t = (0:samples)/fs;
% original signal
s = sin(w*t); N = length(s); l = floor(N/2);
a = applyFFT(s,N);

% simulation parameters
theta = pi;
eta = pi/2;
K = 6;
rho = 0.06;

% microphone data
micData = genDelayData(s,theta,eta,K,rho);
micFreqData = zeros(size(micData));
for k = 1:K
    micFreqData(:,k) = N*applyFFT(micData(:,k),N); 
end

% grid for cost function
Theta = deg2rad(0:359); Eta = deg2rad(0:359);
J = zeros(length(Theta),length(Eta));
for i = 1:length(Theta)
    for j = 1:length(Eta)
        J(i,j) = evalCost(micFreqData, Theta(i), a, Eta(j), K, l, rho);
    end
end

mesh(J), xlabel('Range [Eta]'), ylabel('DOA [Theta]')
[maxR,maxRindx] = max(J);
[~,etaEst] = max(maxR);
thetaEst = maxRindx(etaEst);
disp('Estimate'); 
disp(strcat('Range [Eta]: ',num2str(etaEst)))
disp(strcat('DOA [Theta]: ',num2str(thetaEst)))
