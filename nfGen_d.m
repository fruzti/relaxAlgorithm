function d = nfGen_d(eta,l,N)
% function d = nfGen_d(eta,l,N)
% ---------------------
% d : delay vector
% eta : delays in samples
% l : number of harmonics
% N : number of samples

    w0 = 2*pi/N;
    d = exp(-1i*(-l:l)*eta*w0).';
    
end