function tests = symmetricDft2DftTest()
    tests = functiontests(localfunctions);
end

% Test for an even number of samples
function testTransformationEvenSamples(testCase)
    nDataList = [2,4,10,12,50,100,200];
    nNDataList = length(nDataList);
    for iNData = 1:nNDataList
        nData = nDataList(iNData);
        signal = randn(nData,1);
        dftSignal = fft(signal);
        symmetricDftSignal = dft2SymmetricDft(dftSignal);
        actDftSignal = symmetricDft2Dft(symmetricDftSignal,nData);
        testCase.assertEqual(actDftSignal,dftSignal)
    end
end

% Test for an uneven number of samples
function testTransformationUnevenSamples(testCase)
    nDataList = [1,3,11,13,51,101,201];
    nNDataList = length(nDataList);
    for iNData = 1:nNDataList
        nData = nDataList(iNData);
        signal = randn(nData,1);
        dftSignal = fft(signal);
        symmetricDftSignal = dft2SymmetricDft(dftSignal);
        actDftSignal = symmetricDft2Dft(symmetricDftSignal,nData);
        testCase.assertEqual(actDftSignal,dftSignal)
    end
end