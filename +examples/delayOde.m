zDomain = quantity.Domain("z", linspace(0, 1, 61));
tDomain = quantity.Domain("t", linspace(0, 25, 5001));
syms z t;

Lambda = quantity.Symbolic(diag([3]), zDomain, "name", "Lambda");
% A = quantity.Symbolic([0, 0.8; 0.4, 0], zDomain, "name", "A");
% Q0 = 2;
% Q1 = 1;
% 
% G1 = [0; 1]*quantity.Symbolic(1+z, zDomain, "name", "G1");
% G2 = [1];
% G3 = [2];
% G4 = [0.4];
% 
% output = model.Output("controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

ode = ss(eye(2), [1; 0], [0, 1], []);

% mas = MultiAgent(adjacencyMatrix, Lambda);
% 
% simulationSetting = {'t', tDomain.grid};
% % simData = mas.simulate(simulationSetting);
% simData = tool.simulate(simulationSetting, mas);