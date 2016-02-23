%% myIfft
% Zero pads in the frequency domain so that time-domain interpolation is
% achieved.
% 
%% Syntax:
%# interpolatedData = myIfft(dftData,nIdft)
%
%% Description:
% Performs interpolation in the time-domain by zero padding a DFT vector in
% the "middle". The standard ifft performs the zero padding in the end of a
% DFT vector, but this gives an asymmetric frequency response and produces
% a complex valued time domain signal even if the DFT vector was computed
% from a real-valued signal.
%
% * dftData: A vector containing the dft of some data.
% * nIdft: (optional) The number of samples in the interpolated data
%   vector. The default is the length of the dftData vector.
% * interpolatedData: The interpolated data vector.
%
%% Examples: 
% % Increase the sampling rate by a factor of ten for a simple sinusoid
% freq = 10;
% amp = 2;
% phase = 0.1;
% nData = 100;
% data = amp*sin(2*pi*freq*(0:nData-1)'/nData+phase)
% nIdft = 10*nData;
% interpolatedData = myIfft(fft(data),nSubIntervals);
%
%% See also:
% ifft fft
%
function interpolatedData = myIfft(dftData,nIdft)
    [nDft,nSignals] = size(dftData);
    if nargin == 1
        nIdft = nDft;
    end
    if nIdft == nDft
        % If no interpolation is required, return the normal ifft
        interpolatedData = ifft(dftData);
    else
        % Perform interpolation
        if mod(nDft,2)==0
            % nDft is an even number
            lowerDft = dftData((0:nDft/2)+1,:);
            lowerDft(end,:) = lowerDft(end,:)/2;
            upperDft = [lowerDft(end,:);...
                dftData((nDft/2+1:nDft-1)+1,:)];
            zeroVector = zeros(nIdft-nDft-1,nSignals);
        else
            % nDft is an uneven number
            lowerDft = dftData((0:floor(nDft/2))+1,:);
            upperDft = ...
                dftData((floor(nDft/2)+1:nDft-1)+1,:);
            zeroVector = zeros(nIdft-nDft,nSignals);
        end
        interpolatedData = ifft([lowerDft;zeroVector;upperDft])*nIdft/nDft;
    end
end

