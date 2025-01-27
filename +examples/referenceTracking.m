zDomain = quantity.Domain("z", linspace(0, 1, 201));
tDomain = quantity.Domain("t", linspace(0, 20, 4001));
syms z t;
Lambda = quantity.Symbolic(diag([2, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 1; 1, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("controlOutput", "C0", [1 0]);

disturbanceInput = model.Input("plant.disturbance", ...
    "B", [1, 0; 0, 0]*quantity.Symbolic(1 + z, zDomain, "name", "G1"), ...
    "B0", [1, 0], "B1", [1, 0]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "disturbanceInput", disturbanceInput, "Q0", Q0, "Q1", Q1);

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic([1; 0], zDomain),...
    "agent3.x", quantity.Symbolic([2; -2], zDomain),...
    "agent4.x", quantity.Symbolic([1; 0], zDomain)};
S_r = [0 1; -1 0];
[systemTarget, feedforwardControlInput, Transformations] = mas.feedforwardControl(S_r, [1 0], [0 1; -1 0], eye(2));


w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
exogenousSignals = {"reference", w};
ci = {};
for indAgent = 1:mas.N
    ci{end + 1} = "agent" + indAgent;
    ci{end + 1} = feedforwardControlInput.("agent" + indAgent);
end

simDataControl = mas.simulate(simulationSetting, ic, feedforwardControlInput, exogenousSignals);

for indAgent = 1:mas.N
    simDataControl.("agent"+indAgent).reference = w;
end

mas.plotState(simDataControl);