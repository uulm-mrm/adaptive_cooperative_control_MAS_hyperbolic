function [tests] = testMultiAgent()
tests = functiontests(localfunctions);
end % testMultiAgent()

function testConstructorNoFailure(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([2+z, -(0.5+exp(-z))]), zDomain, "name", "Lambda");
A = quantity.Symbolic([-1, z/4; -1, -3*z/2], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("controlOutput", "C", quantity.Symbolic([1 0; 0 1 ], zDomain, "name", "C"));

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

failureHappened = false;
try
	MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1);
catch
	failureHappened = true;
end

tc.verifyFalse(failureHappened);
end % testConstructorNoFailure()

function testSimulateNoFailure(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([2+z, -(0.5+exp(-z))]), zDomain, "name", "Lambda");
A = quantity.Symbolic([-1, z/4; -1, -3*z/2], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("controlOutput", "C", quantity.Symbolic([1 0; 0 1 ], zDomain, "name", "C"));

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1);

simulationSetting = {'t', linspace(0, 1, 201)};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic(ones(2, 1), zDomain),...
    "agent3.x", quantity.Symbolic([2; -2], zDomain),...
    "agent4.x", quantity.Symbolic([1; 0], zDomain)};
[simData, myStateSpace] = mas.simulate(simulationSetting, ic);

failureHappened = false;
try
	mas.simulate(simulationSetting, ic);
catch
	failureHappened = true;
end

tc.verifyFalse(failureHappened);
end % testSimulateNoFailure

function testPlotStateNoFailure(tc)
storeCurrentFigures = findobj(0, 'type', 'figure');
set(groot,'defaultFigureVisible','off')
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([2+z, -(0.5+exp(-z))]), zDomain, "name", "Lambda");
A = quantity.Symbolic([-1, z/4; -1, -3*z/2], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("controlOutput", "C", quantity.Symbolic([1 0; 0 1 ], zDomain, "name", "C"));

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1);

simulationSetting = {'t', linspace(0, 1, 201)};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic(ones(2, 1), zDomain),...
    "agent3.x", quantity.Symbolic([2; -2], zDomain),...
    "agent4.x", quantity.Symbolic([1; 0], zDomain)};
[simData, myStateSpace] = mas.simulate(simulationSetting, ic);

failureHappened = false;
try
    mas.plotState(simData);
catch
	failureHappened = true;
end
set(groot,'defaultFigureVisible','on') 
delete(setdiff(findobj(0, 'type', 'figure'), storeCurrentFigures));
tc.verifyFalse(failureHappened);
end % testPlotStateNoFailure()

function testPlotOutputNoFailure(tc)
storeCurrentFigures = findobj(0, 'type', 'figure');
set(groot,'defaultFigureVisible','off')
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([2+z, -(0.5+exp(-z))]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, z/4; -1, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("controlOutput", "C", quantity.Symbolic([1 0; 0 1 ], zDomain, "name", "C"));

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1);

simulationSetting = {'t', linspace(0, 1, 201)};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic(ones(2, 1), zDomain),...
    "agent3.x", quantity.Symbolic([2; -2], zDomain),...
    "agent4.x", quantity.Symbolic([1; 0], zDomain)};
[simData, myStateSpace] = mas.simulate(simulationSetting, ic);

failureHappened = false;
try
    mas.plotState(simData);
catch
	failureHappened = true;
end
set(groot,'defaultFigureVisible','on') 
delete(setdiff(findobj(0, 'type', 'figure'), storeCurrentFigures));
tc.verifyFalse(failureHappened);
end % testPlotOutputNoFailure()

function testPI_vSolvesRegulatorEquations(tc)
    zDomain = quantity.Domain("z", linspace(0, 1, 201));
    syms z t;
    Lambda = quantity.Symbolic(diag([1+z, -1]), zDomain, "name", "Lambda");
    A = quantity.Symbolic([0, 1+z; 1, 0], zDomain, "name", "A");
    Q0 = 1;
    Q1 = 0;

    output = model.Output("plant.controlOutput", "C1", [1 1]);

    adjacencyMatrix = [0 1; 1 0];
    G1 = [0, 1; 1, 0]*quantity.Symbolic(1+z, zDomain, "name", "G1");
    G2 = [1, 0];
    G3 = [0, 1];
    G4 = [0, 0];

    mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
        "G1", G1, "G2", G2, "G3", G3, "G4", G4);

    S = [0 1; -1 0];
    P = eye(2);
    mas.setBacksteppingKernel();
    mas.setRegulatorEquationDistCoefficients(size(S,1));
    PI_v = mas.solveRegulatorEquationsDisturbance(S, P);
    shouldBeZero = mas.Lambda*PI_v.diff()+mas.A0Backstepped*mas.E1.'*PI_v.at(0)-PI_v*S + mas.G1Backstepped*P;
    tc.verifyEqual(shouldBeZero.median(), zeros(size(PI_v)), "AbsTol", 10^-3);
    tc.verifyEqual((mas.Q0*mas.E1.'-mas.E2.')*PI_v.subs("z", 0), - G2*P, "AbsTol", 10^-3);
    tc.verifyEqual(mas.outputBackstepped.out(PI_v), -G4*P, "AbsTol", 10^-3);
end % testPI_vSolvesRegulatorEquations

function testPI_wSolvesRegulatorEquations(tc)
    S_r = [0 1; -1 0];

    zDomain = quantity.Domain("z", linspace(0, 1, 201));
    syms z t;
    Lambda = quantity.Symbolic(diag([2, -1]), zDomain, "name", "Lambda");
    A = quantity.Symbolic([0, 1; 1, 0], zDomain, "name", "A");
    Q0 = 1;
    Q1 = 1;

    output = model.Output("controlOutput", "C0", [1 0]);

    adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

    mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1);
    mas.setBacksteppingKernel();
    mas.setRegulatorEquationCoefficients(size(S_r,1));
    P_r = [1 0];
    PI_w = mas.solveRegulatorEquationsReference(S_r, P_r);
    
       testVal = Lambda * diff(PI_w, "z") + mas.A0Backstepped()*mas.E1.'*subs(PI_w, "z", 0) - PI_w * S_r;
       tc.verifyEqual(max(median(testVal), [], 'all'), 0, "AbsTol", 1e-4);
       tc.verifyEqual(max(median(mas.E2.'*subs(PI_w, "z", 0) - mas.Q0*mas.E1.'*subs(PI_w, "z", 0)), [], 'all'), 0, "AbsTol", 1e-4);
       tc.verifyEqual(max(median(mas.outputBackstepped.out(PI_w) - P_r), [], 'all'), 0, "AbsTol", 1e-4);
end % testPI_wSolvesRegulatorEquations

function testReferenceTracking(tc)
    tDomain = quantity.Domain("t", linspace(0, 10, 2001));
    adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
    S_r = [0 1; -1 0];
    w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);

    zDomain = quantity.Domain("z", linspace(0, 1, 201));
    syms z t;
    Lambda = quantity.Symbolic(diag([2, -1]), zDomain, "name", "Lambda");
    A = quantity.Symbolic([0, 1; 1, 0], zDomain, "name", "A");
    Q0 = 1;
    Q1 = 1;

    output = model.Output("controlOutput", "C0", [1 0]);
    
    adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
    mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1);
    mas.setBacksteppingKernel();
    mas.setRegulatorEquationCoefficients(size(S_r,1));
    [systemTarget, controller] = mas.feedforwardControl(S_r, [1 0], [0 1; -1 0], eye(2));
    
    exogenousSignals = {"agent1.reference", w, "agent2.reference", w, "agent3.reference", w, "agent4.reference", w};
    simulationSetting = {'t', tDomain.grid};
    ic = {"agent1.x", quantity.Symbolic([0; 0], zDomain),...
        "agent2.x", quantity.Symbolic([0; 0], zDomain),...
        "agent3.x", quantity.Symbolic([0; 0], zDomain),...
        "agent4.x", quantity.Symbolic([0; 0], zDomain)};
    simDataControl = mas.simulate(simulationSetting, ic, controller, exogenousSignals);
    for indAgent = 1:mas.N
       tc.verifyEqual(simDataControl.("agent"+indAgent).controlOutput.at(10) - [1 0]*w.at(10), 0, "AbsTol", 1e-2);
       tc.verifyEqual(simDataControl.("agent"+indAgent).controlOutput.at(7.5) - [1 0]*w.at(7.5), 0, "AbsTol", 1e-2);
       tc.verifyEqual(simDataControl.("agent"+indAgent).controlOutput.at(8) - [1 0]*w.at(8), 0, "AbsTol", 1e-2);
       tc.verifyEqual(simDataControl.("agent"+indAgent).controlOutput.at(5) - [1 0]*w.at(5), 0, "AbsTol", 1e-2);
    end
end % testReferenceTracking