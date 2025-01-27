zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
% Lambda = quantity.Symbolic(diag([2+z, -(0.5+exp(-z))]), zDomain, "name", "Lambda");
Lambda = quantity.Symbolic(diag([2, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 1; 1, 0], zDomain, "name", "A")*0;
Q0 = 1;
Q1 = 1;

G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 0];
G3 = [1, 0];
G4 = [0, 0];

output = model.Output("controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "G1", G1, "G2", G2, "G3", G3, "G4", G4, "Q0", Q0, "Q1", Q1);

simulationSetting = {'t', linspace(0, 5, 1001)};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic(ones(2, 1), zDomain),...
    "agent3.x", quantity.Symbolic([2; -2], zDomain),...
    "agent4.x", quantity.Symbolic([1; 0], zDomain)};
[simData, myStateSpace] = mas.simulate(simulationSetting, ic);

mas.plotState(simData);
% mas.plotOutput(simData);

% [systemTarget, feedforwardControlInput] = mas.feedforwardControl(0, 1, [0 1; 0 0], eye(2));
% 
% simDataTarget = systemTarget.simulate("t", linspace(0, 5, 1001), ...
%     systemTarget.stateName.pde, ic{6});
% 
% simDataTarget.plant.x.plot();
% simDataTarget.plant.controlOutput.plot();
% 
% simDataControl = mas.simulate(simulationSetting, ic, feedforwardControlInput);
% 
% mas.plotState(simDataControl);
% mas.plotOutput(simDataControl);