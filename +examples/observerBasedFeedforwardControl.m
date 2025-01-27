simTime = 30;
zDomain = quantity.Domain("z", linspace(0, 1, 201));
tDomain = quantity.Domain("t", linspace(0, simTime, 5001));
syms z t;
%% Define the System
% Lambda = quantity.Symbolic(diag([2+z, -(0.5+exp(-z))]), zDomain, "name", "Lambda");
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, -1; 0.2, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("controlOutput", "C0", [1 0]);

disturbanceInput = model.Input("disturbance", ...
    "B", [1, 0; 0, 0]*quantity.Symbolic(1 + z, zDomain, "name", "G1"), ...
    "B0", [1, 0], "B1", [1, 0]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
G1 = [0, 1; 1, 0]*quantity.Symbolic(1 + z, zDomain, "name", "G1");
G2 = [0, 0];
G3 = [0, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

%% Define Reference Observer and observe S_r

splitDomain = tDomain.split(simTime/2);
% S_r = quantity.Piecewise({[0 1; -1 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 2; -2 0] * quantity.Discrete.ones(1, splitDomain(2))});
% S_r = quantity.Discrete(S_r);
% w = quantity.Piecewise({quantity.Discrete([sin(splitDomain(1).grid), cos(splitDomain(1).grid)], splitDomain(1)), ...
%     quantity.Discrete([sin(2*splitDomain(2).grid), cos(2*splitDomain(2).grid)], splitDomain(2))});
% w = quantity.Discrete(w);
S_r = [0 1; -1 0];
w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);

arf = ReferenceObserver(adjacencyMatrix, 2, 4, 4, 4, 4);

simDataRef = arf.simulate(w, S_r);

for indAgent = 1:mas.N
    S_rObs.("agent" + indAgent) = simDataRef.("agent" + (indAgent-1)).S;
end

%exogenousSignals = {"agent1.reference", w, "agent2.reference", w, "agent3.reference", w, "agent4.reference", w};
% simDataControl = mas.simulate(simulationSetting, ic, controller, exogenousSignals);
% 
% mas.plotOutput(simDataControl);

%% Design a controller and disturbance observer

S = [0 1; -1 0];
P = eye(2);
P_r = [1 0];

[systemTarget, controller] = mas.feedforwardControl(S_rObs, P_r, S, P);

distObserver = DisturbanceObserver(mas, S, P, "c_", 8);

%% Run the simulation

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent3.x", quantity.Symbolic([0; 0], zDomain),...
    "agent4.x", quantity.Symbolic([0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent3.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent4.observer.x", quantity.Symbolic([0; 0], zDomain)};
v = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
zeroSig = quantity.Discrete.zeros([2, 1], tDomain);
% w = zeroSig;
exogenousSignals = {"disturbance", [zeroSig; zeroSig; zeroSig; v], "agent1.reference", w, "agent2.reference", w, "agent3.reference", w, "agent4.reference", w};
% exogenousSignals = {"disturbance", [zeroSig; zeroSig; zeroSig; v]};

simData = mas.simulate(simulationSetting, ic, controller, exogenousSignals, distObserver);
% exogenousSignals = {"disturbance", [zeroSig; zeroSig; zeroSig; v], "disturbanceState", [zeroSig; zeroSig; zeroSig; v], "agent1.reference", w, "agent2.reference", w, "agent3.reference", w, "agent4.reference", w};
% simData = mas.simulate(simulationSetting, ic, controller, exogenousSignals);

mas.plotOutput(simData);

figure;
ax = subplot(2, 1, 1);
simDomain = simData.agent4.observer.disturbanceState.domain;
plot(simDomain.grid, simData.agent4.observer.disturbanceState(1).valueDiscrete);
hold on; 
plot(simDomain.grid, v(1).valueDiscrete, '--g');
hold off;

ax = subplot(2, 1, 2);
plot(simDomain.grid, simData.agent4.observer.disturbanceState(2).valueDiscrete);
hold on; 
plot(simDomain.grid, v(2).valueDiscrete, '--g');
hold off;

