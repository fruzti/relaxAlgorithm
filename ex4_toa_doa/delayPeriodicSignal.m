%% delayPeriodicSignal
% Delays a periodic signal with a number of samples.
% 
%% Syntax:
%# delayedSignal = delayPeriodicSignal(signal,delay)
%
%% Description:
% Delays a periodic signal with a number of samples. The delay does not
% have to be an integer number, but can also be fractional. The signal x is
% decomposed into
%
% $$x(n) = \sum_{k=-N/2}{N/2}A[k]\exp(j2\pi kn/N)$$
%
% where N are the number of samples. When N is even, A[k] is given by
%
% $$A[k] = 
%     \begin{cases}
%         X[N/2]/2 & k = -N/2\\
%         X[k+N] & -N/2 < k < 0\\
%         X[k] & 0 <= k < N/2
%         X[N/2]/2 & k = N/2\\
%     \end{cases}
% $$
%
% where X[k] is the DFT of x, and given by
%
% $$ X[k] = \sum_{n=0}^{N-1} x(n)\exp(-j2\pi kn/N)$$
%
% for k = 0,1,...,N-1. When N is uneven, A[k] is given by
%
% $$A[k] = 
%     \begin{cases}
%         X[k+N] & -(N-1)/2 <= k < 0\\
%         X[k] & 0 <= k <= (N-1)/2
%     \end{cases}
% $$
%
% Using this decomposition of x, the signal is then delayed by the desired
% number of samples. Note that the DFT decomposition of the signal does not
% produce a useful result for fractional delays since a real-valued
% signal ends up being complex-valued when shifted with a fractional number
% of samples.
%
% * signal: The periodic signal to be delayed
% * delay: The delay in samples. Does not have to be an integer
% * delayedSignal: The delayed signal
%
%% Examples:
% nData = 100
% time = (0:nData-1)';
% signal = sin(2*pi*3*time/nData);
% delay = 3.5;
% delayedSignal = delayPeriodicSignal(signal,delay);
%
%% See also:
% 
%
function delayedSignal = delayPeriodicSignal(signal,delay)
    nData = length(signal);
    nDataIsEven = mod(nData,2) == 0;
    dftSignal = fft(signal(:));
    % Create the symmetric representation of the signal
    symmetricDftSignal = dft2SymmetricDft(dftSignal);
    % Shift the signal in the frequency domain
    if nDataIsEven
        dftIndices = (-nData/2:nData/2)';
    else
        dftIndices = (-(nData-1)/2:(nData-1)/2)';
    end
    symmetricDftDelayedSignal = ...
        symmetricDftSignal.*exp(-1i*2*pi*dftIndices*delay/nData);
    % Transform the signal back to the time domain
    dftDelayedSignal = symmetricDft2Dft(symmetricDftDelayedSignal,nData);
    delayedSignal = ifft(dftDelayedSignal);
end