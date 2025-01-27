zDomain = quantity.Domain("z", linspace(0, 1, 61));
syms z t;

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 1; 1, 0], zDomain, "name", "A");
Q0 = -1;
Q1 = 1;

G1 = [0; 0]*quantity.Symbolic(1+z, zDomain, "name", "G1");
G2 = [0];
G3 = [1];
G4 = [0];

controlInput = model.Input("control", "B1", eye(1));
disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
% controlOutput = model.Output("controlOutput", "C0", [1 0]);
% measurementOutput = model.Output("measurement", "C1", [0 1]);
% output = controlOutput + measurementOutput;
output = model.Output("controlOutput", "C0", [1 0]);

agents{1} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent1");
Lambda = quantity.Symbolic(diag([2, -1]), zDomain, "name", "Lambda"); 
agents{2} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent2");
Lambda = quantity.Symbolic(diag([2, -2]), zDomain, "name", "Lambda"); 
agents{3} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent3");
Lambda = quantity.Symbolic(diag([1, -2]), zDomain, "name", "Lambda"); 
agents{4} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent4");

mas = MultiAgentHeterogenous(agents, adjacencyMatrix, "nCoefDist", 5, "nCoefRef", 5);

tDomain = quantity.Domain("t", linspace(0, 15, 800));
simulationSetting = {'t', tDomain.grid};
P = [1 0];
S = [0 1; -1 0];
Pr = [1 0];
Sr = [0 1; 0 0];
% w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain)*0;
w = quantity.Discrete([1+tDomain.grid, ones(length(tDomain.grid), 1)], tDomain);
v = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
aggregatedDisturbance = [P*v*0; P*v; P*v*0; P*v/2];

arf = ReferenceObserver(adjacencyMatrix, 2, 4, 4, 4, 4, 1, 2);
ic = {"agent1.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent3.x", quantity.Symbolic([0; 0], zDomain),...
    "agent4.x", quantity.Symbolic([1; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent3.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent4.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "S", [0 1; -1 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0],...
    "v", [0; 0; 0; 0; 0; 0; 0; 0]};
exogenousInput = {"disturbance", aggregatedDisturbance, "S", S, "Sr", Sr, "w", w, "Pr", Pr, "P", P};
simDataNew = mas.adaptiveControlSimulation(simulationSetting, 2, 1, arf, "initialCondition", ic, "exogenousInput", exogenousInput, "mu_S", 2, ...
	"distCoef", 5, "adaptationDelayObserver", 10, "adaptationDelayController", 10);

for idx = 1:4
	simDataNew.("agent" + idx).controlOutput.plot();
end