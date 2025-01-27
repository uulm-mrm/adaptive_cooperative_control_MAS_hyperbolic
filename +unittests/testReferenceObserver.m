function [tests] = testReferenceObserver()
tests = functiontests(localfunctions);
end % testMultiAgent()

function testConstructorNoFailure(tc)
adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

failureHappened = false;
try
	arf = ReferenceObserver(adjacencyMatrix, 2, 6, 6, 2, 2);
catch
	failureHappened = true;
end

tc.verifyFalse(failureHappened);
end % testConstructorNoFailure()

function testSimulateNoFailure(tc)
tDomain = quantity.Domain("t", linspace(0, 10, 2001));
adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
S_r = [0 1; 0 0];
w = quantity.Discrete([(tDomain.grid + 2), ones(tDomain.n, 1)], tDomain);

arf = ReferenceObserver(adjacencyMatrix, 2, 6, 6, 2, 2);

failureHappened = false;
try
	simData = arf.simulate(w, S_r);
catch
	failureHappened = true;
end

tc.verifyFalse(failureHappened);
end % testSimulateNoFailure()

function testPlotObservationNoFailure(tc)
storeCurrentFigures = findobj(0, 'type', 'figure');
set(groot,'defaultFigureVisible','off')
tDomain = quantity.Domain("t", linspace(0, 10, 2001));
adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
S_r = [0 1; 0 0];
w = quantity.Discrete([(tDomain.grid + 2), ones(tDomain.n, 1)], tDomain);

arf = ReferenceObserver(adjacencyMatrix, 2, 6, 6, 2, 2);
simData = arf.simulate(w, S_r);

failureHappened = false;
try
	arf.plotObservation(simData);
catch
	failureHappened = true;
end
set(groot,'defaultFigureVisible','on') 
delete(setdiff(findobj(0, 'type', 'figure'), storeCurrentFigures));
tc.verifyFalse(failureHappened);
end % testPlotObservationNoFailure()