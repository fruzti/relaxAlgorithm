function tests = fibonacciTest
    tests = functiontests(localfunctions);
end

% Validate the produced numbers
function testReturnedNo(testCase)
    validIndexArray = {...
        0,...
        1,...
        2,...
        -1,...
        -7,...
        10,...
        nan,...
        -7:7,...
        70,...
        10*ones(3,3),...
        };
    expFibonacciNoArray = {...
        0,...
        1,...
        1,...
        1,...
        13,...
        55,...
        nan,...
        [13,-8,5,-3,2,-1,1,0,1,1,2,3,5,8,13],...
        190392490709135,...
        55*ones(3,3),...
        };
    invalidIndexArray = {...
        inf,...
        -inf,...
        3*1i,...
        1.1,...
        [-10,8,nan,inf,-inf,7,-2.01],...
        -71,...
        72,...
        };
    indexArray = [validIndexArray, invalidIndexArray];
    nExpFibonacciNoArray = length(expFibonacciNoArray);
    nIndexArray = length(indexArray);
    for iExpNo = 1:nIndexArray
        if iExpNo <= nExpFibonacciNoArray
            actNo = fibonacci(indexArray{iExpNo});
            testCase.assertEqual(actNo,expFibonacciNoArray{iExpNo});
        else
            testCase.assertError(@()fibonacci(indexArray{iExpNo}),...
                'fibonacci:argChk');
        end
    end
end