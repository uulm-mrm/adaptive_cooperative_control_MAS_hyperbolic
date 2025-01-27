tDomain = quantity.Domain("t", linspace(0, 10, 2001));
adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
%adjacencyMatrix = [0 0; 1 0];
% S_r = [0 1; 0 0];
% w = quantity.Discrete([(tDomain.grid + 2), ones(tDomain.n, 1)], tDomain);
S_r = [0 1; -1 0];
w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);

arf = ReferenceObserver(adjacencyMatrix, 2, 10, 10, 10, 10);

simData = arf.simulate(w, S_r);
% arf.plotObservation(simData);


zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
% Lambda = quantity.Symbolic(diag([2+z, -(0.5+exp(-z))]), zDomain, "name", "Lambda");
Lambda = quantity.Symbolic(diag([2, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 1; 1, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("controlOutput", "C0", [1 0]);

% disturbanceInput = model.Input("disturbance", ...
%     "B", [1, 0; 0, 0]*quantity.Symbolic(1 + z, zDomain, "name", "G1"), ...
%     "B0", [1, 0], "B1", [1, 0]);
disturbanceInput = model.Input("disturbance", ...
    "B", [1, 0; 0, 0]*quantity.Symbolic(1 + z, zDomain, "name", "G1"), ...
    "B0", [1, 0], "B1", [1, 0]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1);

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic(ones(2, 1), zDomain),...
    "agent3.x", quantity.Symbolic([2; -2], zDomain),...
    "agent4.x", quantity.Symbolic([1; 0], zDomain)};

for indAgent = 1:mas.N
    S_rObs.("agent" + indAgent) = simData.("agent" + (indAgent-1)).S;
end

[systemTarget, controller] = mas.feedforwardControl(S_rObs, [1 0], [0 1; -1 0], eye(2));
% [systemTarget, feedforwardControlInput, Transformations] = mas.feedforwardControl(S_r, [1 0], [0 1; -1 0], eye(2));

% exogenousSignals = {"reference", w, "disturbance", [w; w; w; w]};
exogenousSignals = {"agent1.reference", w*0, "agent2.reference", w*0, "agent3.reference", w*0, "agent4.reference", w};

simDataControl = mas.simulate(simulationSetting, ic, controller, exogenousSignals);
% simDataControl = mas.simulate(simulationSetting, ic, feedforwardControlInput);

% mas.plotState(simDataControl);
mas.plotOutput(simDataControl);
% simDataTarget.plant.controlOutput.plot();