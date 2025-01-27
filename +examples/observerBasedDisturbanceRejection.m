simTime = 15;
zDomain = quantity.Domain("z", linspace(0, 1, 101));
tDomain = quantity.Domain("t", linspace(0, simTime, 2501));
syms z t;
adjacencyMatrix = [0 0 0 0; 1 0 0 0; 1 0 0 0; 1 0 0 0];
%% Define the System
Lambda = quantity.Symbolic(diag([1+z/3, 0.9*cos(0.8*z), -1, -3*exp(-z)]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0.8, -1+z/4, 0.2; 1/(z+1), 0, -0.4, 0; 0, 0, 0, -1; -2, z/2, 0, 0], zDomain, "name", "A");
Q0 = [1, 0.2; 0 1];
Q1 = [1 0; 0 0.6];

output = model.Output("controlOutput", "C0", [1 0 0 0; 0 1 0 0]);
G1 = [0, 1; 0, 0; 0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 0; 1 0];
G3 = [1, 0; 0 0];
G4 = [0, 0; 0 1];
% G1 = [0; 1; 0; 0]*quantity.Symbolic(1, zDomain, "name", "G1");
% G2 = [0; 0];
% G3 = [1; 0];
% G4 = [0; 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);
%% Design a controller and disturbance observer

S = [0 1; -1 0];
P = eye(2);
% S = 0;
% P = eye(1);

distOde = ss(S, zeros(size(S,1)), P, []);
v = stateSpace.simulate(distOde, zeros(size(S,1), length(tDomain.grid)), tDomain.grid, [0;1]);
v = quantity.Discrete(v, tDomain);

S_r = [0 1; -1 0];
P_r = eye(2);
[systemTarget, controller] = mas.feedforwardControl(S_r, P_r, S, P);

distObserver = DisturbanceObserver(mas, S, P, "c_", 6);

%% Run the simulation

simulationSetting = {'t', tDomain.grid};
ic = {};
for indAgent = 1:mas.N
    ic{end+1} = "agent"+indAgent+".x";
    ic{end+1} = quantity.Symbolic(zeros(size(Lambda,1),1), zDomain);
    ic{end+1} = "agent"+indAgent+".observer.x";
    ic{end+1} = quantity.Symbolic(zeros(size(Lambda,1),1), zDomain);
end
zeroSig = quantity.Discrete.zeros([size(S,1), 1], tDomain);
% w = zeroSig;
% exogenousSignals = {"disturbance", [zeroSig; zeroSig; zeroSig; v], "agent1.reference", w, "agent2.reference", w, "agent3.reference", w, "agent4.reference", w};
% exogenousSignals = {"disturbance", [v; zeroSig], "agent1.reference", w, "agent2.reference", w};
% exogenousSignals = {"disturbance", [zeroSig; zeroSig; zeroSig; v]};
dAgg = [zeroSig; v; zeroSig; v];
exogenousSignals = {"disturbance", dAgg};
% exogenousSignals = {"disturbance", dAgg, "disturbanceState", dAgg};


simData = mas.simulate(simulationSetting, ic, controller, exogenousSignals, distObserver);
% simData = mas.simulate(simulationSetting, ic, controller, exogenousSignals, struct([]));
% simData = mas.simulate(simulationSetting, ic, struct(), exogenousSignals, distObserver, true);
% exogenousSignals = {"disturbance", [zeroSig; zeroSig; zeroSig; v], "disturbanceState", [zeroSig; zeroSig; zeroSig; v], "agent1.reference", w, "agent2.reference", w, "agent3.reference", w, "agent4.reference", w};
% simData = mas.simulate(simulationSetting, ic, controller, exogenousSignals);

mas.plotOutput(simData);
figure;
for indAgent = 1:mas.N
    for indDim = 1:length(v)
        ax = subplot(mas.N, length(v), (indAgent - 1)*length(v) + indDim);
        simDomain = simData.("agent" + indAgent).observer.disturbanceState.domain;
        plot(simDomain.grid, simData.("agent" + indAgent).observer.disturbanceState(indDim).valueDiscrete);
        hold on; 
        d = dAgg((indAgent - 1)*length(v) + 1 : indAgent * length(v));
        plot(simDomain.grid, d(indDim).valueDiscrete, '--g');
        hold off;
    end
end

