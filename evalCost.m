function J = evalCost(micData,theta,a,eta,M,l,rho)
% micData : matrix containing the FFT of the data in the array
% theta : angle of arrival to evaluate
% a : fft of sources
% M : number of microphones
% l : number of harmonics
% rho : radius of array
% eta : range to evaluate

    tmp = 0;
    
    for k = 1:M
        tmp = tmp + conj(gen_d(theta,l,rho,M,k)) .* micData(:,k);
    end
    
    tmp = conj(a) .* tmp;
    
    f = gen_f(eta,l);
    
    J = abs(f'*tmp)^2;
    
end