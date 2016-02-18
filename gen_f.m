function f = gen_f(eta,l)
% f : shifted dtft vector
% eta : fundamental frequency
% l : number of harmonics

    f = exp(1i*eta*(l:-l)).';

end