zDomain = quantity.Domain("z", linspace(0, 1, 201));
tDomain = quantity.Domain("t", linspace(0, 5, 1001));
syms z t;
Lambda = quantity.Symbolic(diag([2, -1, -2, -3]), zDomain, "name", "Lambda");
% A = quantity.Symbolic([0, 1, 0, 0; 1, 0, 0, 0; 1, 0, 0, -1; -1, 0, 1, 0], zDomain, "name", "A");
A = quantity.Symbolic([0, 0, 0, 1; 1, 0, 0, 0; 0, -1, 0, 1; 1, 0, -1, 0], zDomain, "name", "A");
Q0 = [1; 1; 1];
Q1 = [1 0 1];

% output = model.Output("controlOutput", "C", [1 0 0 0; 0 0 0 1]*quantity.Discrete.ones(1, zDomain));
output = model.Output("controlOutput", "C0", [1 1 0 0; 0 0 0 1]);

% output = model.Output("controlOutput", "C0", [1 1 1]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
% adjacencyMatrix = [0 0; 1 0];
G1 = [0, 1; 0, 2; 0, 0; -1, 0]*quantity.Symbolic(1 + z, zDomain, "name", "G1");
G2 = [0 , 0; 1 0; 0, 0];
G3 = [1, 0];
G4 = [0, 0; 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(4, 1), zDomain),...
    "agent2.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent3.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent4.x", quantity.Symbolic([0; 0; 0; 0], zDomain)};
S_r = [0 1; -1 0];
S = [0 1; -1 0];
P = eye(2);
[systemTarget, controller] = mas.feedforwardControl(S_r, eye(2), S, P);

%% 

w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
v = quantity.Discrete.zeros([2, 1], tDomain);
exogenousSignals = {'disturbanceState', [v; v; v; w], "disturbance", [v; v; v; w]};

simDataControlBoth = mas.simulate(simulationSetting, ic, controller, exogenousSignals);

for indAgent = 1:mas.N
    simDataTarget = systemTarget.simulate("t", tDomain.grid, ...
        systemTarget.stateName.pde, controller.backsteppingKernel.backsteppingTransformation(ic{2*indAgent}));
    simDataControlBoth.("agent"+indAgent).reference = simDataTarget.plant.controlOutput;
end

mas.plotOutput(simDataControlBoth);
