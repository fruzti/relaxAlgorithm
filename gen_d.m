function d = gen_d(theta,l,rho,K,k)
% d : delay matrix
% theta : direction of arrival
% l : number of harmonics
% rho : radius of UCA
% K : number of microphones
% k : microphone element

    eta = -rho * cos(theta - 2*pi*(k-1)/K);
    d = gen_f(eta,l);
    
end