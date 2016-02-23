function tests = dft2SymmetricDftTest()
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
        % The middle element of symmetricDftSignal must be real-valued and
        % equal to the first element of dftSignal
        testCase.assertTrue(isreal(symmetricDftSignal(nData/2+1)));
        testCase.assertEqual(symmetricDftSignal(nData/2+1),...
            dftSignal(1));
        % The first and last elements of symmetricDftSignal must be real-valued
        % and equal half of the (nData/2+1)'th element of dftSignal
        testCase.assertTrue(isreal(symmetricDftSignal(1)));
        testCase.assertTrue(isreal(symmetricDftSignal(end)));
        testCase.assertEqual(symmetricDftSignal(1),...
            symmetricDftSignal(end));
        testCase.assertEqual(2*symmetricDftSignal(1),...
            dftSignal(nData/2+1));
        % Finally, symmetricDftSignal must be conjugate symmetric
        testCase.assertEqual(symmetricDftSignal,...
            conj(flipud(symmetricDftSignal)));
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
        % The middle element of symmetricDftSignal must be real-valued and
        % equal to the first element of dftSignal
        testCase.assertTrue(isreal(symmetricDftSignal((nData+1)/2)));
        testCase.assertEqual(symmetricDftSignal((nData+1)/2),...
            dftSignal(1));
        % Finally, symmetricDftSignal must be conjugate symmetric
        testCase.assertEqual(symmetricDftSignal,...
            conj(flipud(symmetricDftSignal)));
    end
end