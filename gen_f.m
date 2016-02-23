function f = gen_f(eta,l)
% function f = gen_f(eta,l)
% --------------------------
% f : shifted dtft vector
% eta : fundamental frequency
% l : number of harmonics

    f = exp(1i*eta*(l:-1:-l)).';

end