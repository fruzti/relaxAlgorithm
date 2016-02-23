%% dft2SymmetricDft
% Computes the symmetric DFT representation from the DFT of a real-valued
% signal.
% 
%% Syntax:
%# symmetricDftData = dft2SymmetricDft(dftData)
%
%% Description:
%Computes the symmetric DFT representation given by
%
% $$x(n) = \sum_{k=-N/2}{N/2}A[k]\exp(j2\pi kn/N)$$
%
% where N are the number of samples, from the DFT representation given by
%
% $$x(n) = \sum_{k=0}{N-1}X[k]\exp(j2\pi kn/N)$$
%
% When $N$ is uneven, the $(N-1)/2$ last data point in the vector X is set
% to be the $(N-1)/2$ first points in the vector A. The first $(N+1)/2$
% points in X are the $(N+1)/2$ last points in A. When $N$ is even, the
% same is done as above to the vector X with the $N/2$'th element removed.
% This element (which is real) is first divided by two and then set as the 
% first and last element of A. Consequently, the vector A will always have
% an uneven number of samples.
% 
% * dftData: The DFT of a real-valued signal
% * symmetricDftData: The symmetric DFT representation of a real-valued
%   signal
%
%% Examples:
% dftData = fft(dataVector);
% symmetricDftData = dft2SymmetricDft(dftData);
%
%% See also:
% symmetricDft2Dft
%
function symmetricDftData = dft2SymmetricDft(dftData)
    nData = length(dftData);
    nDataIsEven = mod(nData,2) == 0;
    % Create the symmetric representation of the signal
    if nDataIsEven
        symmetricDftData = [dftData(nData/2+1)/2;dftData(nData/2+2:end);...
            dftData(1:nData/2);dftData(nData/2+1)/2];
    else
        symmetricDftData = ...
            [dftData((nData+1)/2+1:end);dftData(1:(nData-1)/2+1)];
    end