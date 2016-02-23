function n = genPseudoRndNoise(fl,fu,N,fs)
% function n = genPseudoRndNoise(N)
% ------------------------------------
% n : pseudo random noise sequence 
% fl : lower frequency
% fu : upper frequency
% N : length of sequence

%     w0 = 2*pi/N;
    
    lu = round((fu/fs)*N);
    ll = round((fl/fs)*N);
    
%     l = round(N/2);
    
    alpha = zeros(N,1);
    indx = [-lu:-ll ll:lu] + round(N/2);
    numH = length(indx);
    alpha(indx) = exp(1i*2*pi*rand(numH,1));
    
    n = ifft(ifftshift(alpha),'symmetric');

end