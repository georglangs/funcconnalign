function results = run()
% runs all test cases in the fcalign package

    suite = matlab.unittest.TestSuite.fromPackage('fcalign.test');
    results = suite.run();
end
