function [tests] = testMultiAgent()
tests = functiontests(localfunctions);
end % testMultiAgent()

function testConstructorNoFailure(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
tDomain = quantity.Domain("t", linspace(0, 15, 4001));
syms z t;
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 0;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 0];
G3 = [1, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);
failureHappened = false;
try
    distObserver = DisturbanceObserver(mas, S, P);
catch
	failureHappened = true;
end

tc.verifyFalse(failureHappened);
end % testConstructorNoFailure()

function testSimulateNoFailure(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
tDomain = quantity.Domain("t", linspace(0, 15, 4001));
syms z t;
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 0;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 0];
G3 = [1, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);

distObserver = DisturbanceObserver(mas, S, P, "c_", 4);

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain)};
w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
v = quantity.Discrete.zeros([2, 1], tDomain);
exogenousSignals = {"disturbance", [w; v]};

failureHappened = false;
try
    simData = mas.simulate(simulationSetting, ic, struct(), exogenousSignals, distObserver);
catch
	failureHappened = true;
end

tc.verifyFalse(failureHappened);
end % testSimulateNoFailure

function testObserverErrorToZero(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
tDomain = quantity.Domain("t", linspace(0, 15, 4001));
syms z t;
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 0;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 0];
G3 = [1, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);

distObserver = DisturbanceObserver(mas, S, P, "c_", 4);

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain)};
w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
v = quantity.Discrete.zeros([2, 1], tDomain);
exogenousSignals = {"disturbance", [w; v]};
simData = mas.simulate(simulationSetting, ic, struct(), exogenousSignals, distObserver);

tc.verifyEqual(simData.agent1.observer.disturbanceState(1).at(15) - w(1).at(15), 0, "AbsTol", 1e-2);
tc.verifyEqual(simData.agent1.observer.disturbanceState(2).at(15) - w(2).at(15), 0, "AbsTol", 1e-2);
tc.verifyEqual(simData.agent2.observer.disturbanceState(1).at(15), 0, "AbsTol", 1e-2);
tc.verifyEqual(simData.agent2.observer.disturbanceState(2).at(15), 0, "AbsTol", 1e-2);
end % testObserverErrorToZero

function testObserverErrorFollowsTargetErrorSimple(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
tDomain = quantity.Domain("t", linspace(0, 15, 4001));
syms z t;
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 0;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 0];
G3 = [1, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);

distObserver = DisturbanceObserver(mas, S, P, "c_", 4);

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain)};
w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
v = quantity.Discrete.zeros([2, 1], tDomain);
exogenousSignals = {"disturbance", [w; v]};
simData = mas.simulate(simulationSetting, ic, struct(), exogenousSignals, distObserver);

pde2ode = [];
In = eye(size(mas.H, 1));
A = kron(In, S) - kron(mas.H, distObserver.L_d*mas.E1.'*distObserver.Omega.at(0));
B = kron(In, distObserver.L_d);
ev0 = [];
e0 = [];
evOrig = [];

ind0 = 1;
for indAgent = ind0:size(adjacencyMatrix, 1)
   ev0 = [ev0; simData.disturbance((indAgent - 1)*size(P, 1) + 1:indAgent*size(P, 1)).at(0) - simData.("agent"+indAgent).observer.disturbanceState.at(0)];
   e0 = [e0; distObserver.backsteppingKernel.transform(simData.("agent"+indAgent).x.subs("t", 0) - simData.("agent"+indAgent).observer.x.subs("t", 0))];
end
e0Tilde = e0 - kron(mas.H, distObserver.Omega)*ev0;

for indAgent = ind0:size(adjacencyMatrix, 1)
   systemTarget = model.Transport(Lambda, ...
       "bc0", model.Output("plant.bc0", "C0", -mas.E2.'), ...
       "bc1", model.Output("plant.bc1", "C1", Q1*mas.E2.' - mas.E1.', "C", -distObserver.F),...
       "input", model.Input("empty", "B0", zeros(size(G2, 1), 1))); 
   simDataTarget = systemTarget.simulate("t", tDomain.grid, systemTarget.stateName.pde(1), e0Tilde((indAgent - ind0)*mas.n + 1:(indAgent - ind0 + 1)*mas.n));
   pde2ode = [pde2ode; mas.E1.'*simDataTarget.plant.x.subs("z", 0)];
   evOrig = [evOrig; simData.disturbance((indAgent - 1)*size(P, 1) + 1:indAgent*size(P, 1)) - simData.("agent"+indAgent).observer.disturbanceState];
end
ode = ss(A, -B, eye(length(ev0)), []);
evTarget = stateSpace.simulate(ode, [pde2ode.valueDiscrete], tDomain.grid, ev0);

tc.verifyEqual(median(evOrig - quantity.Discrete(evTarget, tDomain)), zeros(length(evOrig), 1), "AbsTol", 1e-2);
end % testObserverErrorFollowsTargetErrorSimple

function testObserverErrorFollowsTargetErrorUnstable2d(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
tDomain = quantity.Domain("t", linspace(0, 15, 4001));
syms z t;
Lambda = quantity.Symbolic(diag([2-z, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 2; -0.5, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 1; 1, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 1];
G3 = [1, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);

distObserver = DisturbanceObserver(mas, S, P, "c_", 4);

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain)};
w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
v = quantity.Discrete.zeros([2, 1], tDomain);
exogenousSignals = {"disturbance", [w; v]};
simData = mas.simulate(simulationSetting, ic, struct(), exogenousSignals, distObserver);

pde2ode = [];
In = eye(size(mas.H, 1));
A = kron(In, S) - kron(mas.H, distObserver.L_d*mas.E1.'*distObserver.Omega.at(0));
B = kron(In, distObserver.L_d);
ev0 = [];
e0 = [];
evOrig = [];

ind0 = 1;
for indAgent = ind0:size(adjacencyMatrix, 1)
   ev0 = [ev0; simData.disturbance((indAgent - 1)*size(P, 1) + 1:indAgent*size(P, 1)).at(0) - simData.("agent"+indAgent).observer.disturbanceState.at(0)];
   e0 = [e0; distObserver.backsteppingKernel.transform(simData.("agent"+indAgent).x.subs("t", 0) - simData.("agent"+indAgent).observer.x.subs("t", 0))];
end
e0Tilde = e0 - kron(mas.H, distObserver.Omega)*ev0;

for indAgent = ind0:size(adjacencyMatrix, 1)
   systemTarget = model.Transport(Lambda, ...
       "bc0", model.Output("plant.bc0", "C0", -mas.E2.'), ...
       "bc1", model.Output("plant.bc1", "C1", Q1*mas.E2.' - mas.E1.', "C", -distObserver.F),...
       "input", model.Input("empty", "B0", zeros(size(G2, 1), 1))); 
   simDataTarget = systemTarget.simulate("t", tDomain.grid, systemTarget.stateName.pde(1), e0Tilde((indAgent - ind0)*mas.n + 1:(indAgent - ind0 + 1)*mas.n));
   pde2ode = [pde2ode; mas.E1.'*simDataTarget.plant.x.subs("z", 0)];
   evOrig = [evOrig; simData.disturbance((indAgent - 1)*size(P, 1) + 1:indAgent*size(P, 1)) - simData.("agent"+indAgent).observer.disturbanceState];
end
ode = ss(A, -B, eye(length(ev0)), []);
evTarget = stateSpace.simulate(ode, [pde2ode.valueDiscrete], tDomain.grid, ev0);

tc.verifyEqual(median(evOrig - quantity.Discrete(evTarget, tDomain)), zeros(length(evOrig), 1), "AbsTol", 1e-2);
end % testObserverErrorFollowsTargetErrorSimple

function testObserverErrorFollowsTargetErrorTimeVariantLambda2d(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 101));
tDomain = quantity.Domain("t", linspace(0, 10, 2001));
syms z t;
Lambda = quantity.Symbolic(diag([3-2*z, -1+z/2]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 2+z/2; -0.5, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 1; 1, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 1];
G3 = [1, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);

distObserver = DisturbanceObserver(mas, S, P, "c_", 4);

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain)};
w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
v = quantity.Discrete.zeros([2, 1], tDomain);
exogenousSignals = {"disturbance", [w; v]};
simData = mas.simulate(simulationSetting, ic, struct(), exogenousSignals, distObserver);

pde2ode = [];
In = eye(size(mas.H, 1));
A = kron(In, S) - kron(mas.H, distObserver.L_d*mas.E1.'*distObserver.Omega.at(0));
B = kron(In, distObserver.L_d);
ev0 = [];
e0 = [];
evOrig = [];

ind0 = 1;
for indAgent = ind0:size(adjacencyMatrix, 1)
   ev0 = [ev0; simData.disturbance((indAgent - 1)*size(P, 1) + 1:indAgent*size(P, 1)).at(0) - simData.("agent"+indAgent).observer.disturbanceState.at(0)];
   e0 = [e0; distObserver.backsteppingKernel.transform(simData.("agent"+indAgent).x.subs("t", 0) - simData.("agent"+indAgent).observer.x.subs("t", 0))];
end
e0Tilde = e0 - kron(mas.H, distObserver.Omega)*ev0;

for indAgent = ind0:size(adjacencyMatrix, 1)
   systemTarget = model.Transport(Lambda, ...
       "bc0", model.Output("plant.bc0", "C0", -mas.E2.'), ...
       "bc1", model.Output("plant.bc1", "C1", Q1*mas.E2.' - mas.E1.', "C", -distObserver.F),...
       "input", model.Input("empty", "B0", zeros(size(G2, 1), 1))); 
   simDataTarget = systemTarget.simulate("t", tDomain.grid, systemTarget.stateName.pde(1), e0Tilde((indAgent - ind0)*mas.n + 1:(indAgent - ind0 + 1)*mas.n));
   pde2ode = [pde2ode; mas.E1.'*simDataTarget.plant.x.subs("z", 0)];
   evOrig = [evOrig; simData.disturbance((indAgent - 1)*size(P, 1) + 1:indAgent*size(P, 1)) - simData.("agent"+indAgent).observer.disturbanceState];
end
ode = ss(A, -B, eye(length(ev0)), []);
evTarget = stateSpace.simulate(ode, [pde2ode.valueDiscrete], tDomain.grid, ev0);

tc.verifyEqual(median(evOrig - quantity.Discrete(evTarget, tDomain)), zeros(length(evOrig), 1), "AbsTol", 1e-2);
end % testObserverErrorFollowsTargetErrorTimeVariantLambda2d
