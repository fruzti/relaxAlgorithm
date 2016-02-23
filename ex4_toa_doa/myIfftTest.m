function tests = myIfftTest()
    tests = functiontests(localfunctions);
end

% Validate the produced test
function testInterpolatedData(testCase)
    freqs = [10;15];
    amp = [2;1.2];
    phase = [0.1;-pi];
    nSubIntervals = 3;
    nDataArray = [100,101];
    nNDataArray = length(nDataArray);
    for iNData = 1:nNDataArray
        index = (0:nDataArray(iNData)-1)'/nDataArray(iNData);
        data = [amp(1)*sin(2*pi*freqs(1)*index+phase(1)),...
            amp(2)*sin(2*pi*freqs(2)*index+phase(2))];
        dftData = fft(data);
        interpolatedDataIndex = ...
            (0:nDataArray(iNData)*nSubIntervals-1)'/...
            (nDataArray(iNData)*nSubIntervals);
        expInterPolatedData = ...
            [amp(1)*sin(2*pi*freqs(1)*interpolatedDataIndex+phase(1)),...
            amp(2)*sin(2*pi*freqs(2)*interpolatedDataIndex+phase(2))];
        nIdft = nDataArray(iNData)*nSubIntervals;
        actInterPolatedData = myIfft(dftData, nIdft);
        keyboard
        testCase.assertEqual(actInterPolatedData,expInterPolatedData,...
            'absTol',1e-10);
    end
end
