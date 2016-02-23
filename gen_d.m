function d = gen_d(theta,l,K,k,p)
% function d = gen_d(theta,l,K,k,p)
%----------------------------------
% d : delay matrix
% theta : direction of arrival
% l : number of harmonics
% K : number of microphones
% k : microphone element
% p : frequency proportional to radius (p = w0*r*fs/c)

    eta = -p * cos(theta - 2*pi*(k-1)/K);
    d = gen_f(eta,l);
    
end