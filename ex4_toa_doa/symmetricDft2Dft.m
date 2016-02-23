%% symmetricDft2Dft
% Computes the DFT representation from the symmetric DFT representation of a 
% real-valued signal.
% 
%% Syntax:
%# dftData = symmetricDft2Dft(symmetricDftData,nData)
%
%% Description:
% Computes the DFT representation from the symmetric DFT representation of a 
% real-valued signal. The symmetric DFT representation is given by
%
% $$x(n) = \sum_{k=-N/2}{N/2}A[k]\exp(j2\pi kn/N)$$
%
% where N are the number of samples. The DFT representation is given by
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
% * symmetricDftData: The symmetric DFT representation of a real-valued
%   signal
% * nData: Number of acquired data points and the number of samples in the
%   DFT representation
% * dftData: The DFT of a real-valued signal
%
%% Examples:
% dftData = fft(dataVector);
% nData = length(dataVector)
% symmetricDftData = dft2SymmetricDft(dftData);
% dftData = symmetricDft2Dft(symmetricDftData,nData);
% dataVector = ifft(dftData);
% 
%
%% See also:
% dft2SymmetricDft
%
function dftData = symmetricDft2Dft(symmetricDftData,nData)
    nDataIsEven = mod(nData,2) == 0;
    % Transform back to the DFT representation
    if nDataIsEven
        dftData = [symmetricDftData(nData/2+1:end-1);...
            symmetricDftData(1)+symmetricDftData(end);...
            symmetricDftData(2:nData/2);];
    else
        dftData = [symmetricDftData((nData+1)/2:end);...
            symmetricDftData(1:(nData-1)/2)];
    end