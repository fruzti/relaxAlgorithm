function tests = delayPeriodicSignalTest()
    tests = functiontests(localfunctions);
end

% test delays for integer valued delays and an even number of samples
function testIntegerDelaysEvenSamples(testCase)
    nData = 100;
    signal = randn(nData,1);
    delayList = [-99,-49,-8,0,2,5,99];
    nDelayList = length(delayList);
    for iDelay = 1:nDelayList
        delayedSignal = delayPeriodicSignal(signal,delayList(iDelay));
        if delayList(iDelay)>0
            inverseDelayedSignal = [delayedSignal(delayList(iDelay)+1:end);...
                delayedSignal(1:delayList(iDelay))];
        else
            inverseDelayedSignal = [delayedSignal(end+delayList(iDelay)+1:end);...
                delayedSignal(1:end+delayList(iDelay))];
        end
        testCase.assertEqual(signal,...
            inverseDelayedSignal,'AbsTol',1e-12);
    end
end

% test delays for integer valued delays and an uneven number of samples
function testIntegerDelaysUnevenSamples(testCase)
    nData = 101;
    signal = randn(nData,1);
    delayList = [-99,-49,-8,0,2,5,99];
    nDelayList = length(delayList);
    for iDelay = 1:nDelayList
        delayedSignal = delayPeriodicSignal(signal,delayList(iDelay));
        if delayList(iDelay)>0
            inverseDelayedSignal = [delayedSignal(delayList(iDelay)+1:end);...
                delayedSignal(1:delayList(iDelay))];
        else
            inverseDelayedSignal = [delayedSignal(end+delayList(iDelay)+1:end);...
                delayedSignal(1:end+delayList(iDelay))];
        end
        testCase.assertEqual(signal,...
            inverseDelayedSignal,'AbsTol',1e-12);
    end
end

% test delays for fractional delays and an even number of samples
function testFractionalDelaysEvenSamples(testCase)
    nData = 100;
    nOversampling = 100;
    nUpData = nData*nOversampling;
    time = (0:nUpData-1)';
    refSignal = real(exp(1i*2*pi*time*(0:nData/2)/nUpData)*...
        [randn(nData/2,2)*[1;1i];randn(1,1)]);
    refSignalExtended = [refSignal;refSignal;refSignal];
    signal = refSignal(1:nOversampling:end);
    delayList = [-99.02,-49.5,-8.01,0.83,2.3,5.41,99.01];
    nDelayList = length(delayList);
    for iDelay = 1:nDelayList
        delay = delayList(iDelay);
        delayOversampled = delay*nOversampling;
        delayedSignalOversampled = ...
            refSignalExtended((1:nUpData)+nUpData-delayOversampled);
        % Down-sample the signals
        expDelayedSignal = delayedSignalOversampled(1:nOversampling:end);
        actDelayedSignal = delayPeriodicSignal(signal,delay);
        testCase.assertEqual(actDelayedSignal,...
            expDelayedSignal,'AbsTol',1e-12);
    end
end

% test delays for fractional delays and an uneven number of samples
function testFractionalDelaysUnevenSamples(testCase)
    nData = 101;
    nOversampling = 100;
    nUpData = nData*nOversampling;
    time = (0:nUpData-1)';
    refSignal = real(exp(1i*2*pi*time*(0:(nData-1)/2)/nUpData)*...
        randn((nData-1)/2+1,2)*[1;1i]);
    refSignalExtended = [refSignal;refSignal;refSignal];
    signal = refSignal(1:nOversampling:end);
    delayList = [-99.02,-49.5,-8.01,0.83,2.3,5.41,99.01];
    nDelayList = length(delayList);
    for iDelay = 1:nDelayList
        delay = delayList(iDelay);
        delayOversampled = delay*nOversampling;
        delayedSignalOversampled = ...
            refSignalExtended((1:nUpData)+nUpData-delayOversampled);
        % Down-sample the signals
        expDelayedSignal = delayedSignalOversampled(1:nOversampling:end);
        actDelayedSignal = delayPeriodicSignal(signal,delay);
        testCase.assertEqual(actDelayedSignal,...
            expDelayedSignal,'AbsTol',1e-12);
    end
end