function J = evalCost(micData,theta,a,eta,K,l,p)
% function J = evalCost(micData,theta,a,eta,K,l,p)
%----------------------------------------------------
% micData : matrix containing the FFT of the data in the array
% theta : angle of arrival to evaluate
% a : fft of sources
% eta : range to evaluate
% K : number of microphones
% l : number of harmonics
% p : frequency-radius equivalence

    tmp = 0;
    
    for k = 1:K
        tmp = tmp + ( conj(gen_d(theta,l,K,k,p)) .* micData(:,k) );
    end
    
    tmp = conj(a) .* tmp;
    
    f = gen_f(eta,l);
    
    J = abs(f'*tmp)^2;
    
end