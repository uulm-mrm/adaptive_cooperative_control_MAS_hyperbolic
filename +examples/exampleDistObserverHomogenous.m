zDomain = quantity.Domain("z", linspace(0, 1, 61));
tDomain = quantity.Domain("t", linspace(0, 10, 1001));
splitDomain = tDomain.split(10);
syms z t;
Lambda = quantity.Symbolic(diag([2, 1, -1-z/2, -3]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0.5, -0.3*(z-1)*exp(z), -0.25*exp(z);...
	0.3*exp(z), 0, 0.25*exp(z), 0.5;...
	0.5, 0.25*exp(z), 0, -0.25*exp(z);...
	-0.25*exp(z), 0.3*(z-1)*exp(z), 0.5, 0], zDomain, "name", "A");
Q0 = [-1, 0; 0, 0.5];
Q1 = [1, 0.5; 0, -0.5];

discOne = quantity.Discrete.ones(1, zDomain);
discOnet = quantity.Discrete.ones(1, tDomain);
discz = quantity.Symbolic(z, zDomain);
G1 = [0, 0; discOne, 0; 0, 0; discOne-discz, 0];
G2 = [1, 0; 0.5, 0];
G3 = [0, 1; -0.5, 0];
G4 = [0, 0; 0, -1];

output = model.Output("controlOutput", "C0", [1 0 0 0; 0 0 0 0], "C", [0 0 0 0; 0 -discOne 0 0]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "G1", G1, "G2", G2, "G3", G3, "G4", G4, "Q0", Q0, "Q1", Q1, "diffusiveDisturbance", false, "nCoefRef", 5, "nCoefDist", 10);

%disturbance model simulation
S = [0 pi; -pi 0];
distOde = ss(S, zeros(size(S,1)), eye(size(S)), []);
v = stateSpace.simulate(distOde, zeros(size(S,1), length(tDomain.grid)), tDomain.grid, [0; 1]);
v = quantity.Discrete(v, tDomain);
S = S*discOnet;
P = eye(2);

%reference model simulation
Sr = [0 1; 0 0];
refOde = ss(Sr, zeros(size(Sr,1)), eye(size(Sr)), []);
w = stateSpace.simulate(refOde, zeros(size(Sr,1), length(tDomain.grid)), tDomain.grid, [0; 1]);
w = quantity.Discrete(w, tDomain);
Pr = [1, 0; 1 0];

simulationSetting = {'t', tDomain.grid};
zero = quantity.Discrete.zeros(size(v), tDomain);
ic = {"agent1.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent2.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent3.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent4.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent3.observer.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent4.observer.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "S", [0 1; -1 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0],...
    "v", [0; 0; 0; 0; 0; 0; 0; 0]};

distObs = AdaptiveDisturbanceObserver(mas, length(S), size(P,1), 2, "targetEvOde", [-1, -2], "nCoef", 10);

aggregatedDisturbance = [P*zero; P*v; -P*v/4; P*v/2];
aggregatedDisturbanceState = [zero; v; v/4; v/2];

arf = ReferenceObserver(adjacencyMatrix, 2, 4, 4, 4, 4, 2, 2, 2);

exogenousInput = {"disturbance", aggregatedDisturbance, "S", S, "Sr", Sr, "w", w, "Pr", Pr, "P", P, "disturbanceState", aggregatedDisturbanceState};
simData = tool.simulate(simulationSetting, mas, "initialCondition", ic, "exogenousInput", exogenousInput,...
    "disturbanceObserver", distObs, "referenceObserver", arf, "controller", true, "Pr", Pr, "adaptationDelayObserver", 10, "adaptationDelayController", 10);

plot(simData.agent4.controlOutput);
plot(simData.agent4.observer.v)

%% 

%Export control output
header = {'t', 'w'};
M = [tDomain.grid, simData.agent1.observer.w(1).valueDiscrete];
for indAgent = 1:mas.N
	for indDim = 1:mas.output.lengthOutput
		header = cat(2, header, ("y"+indAgent)+indDim);
		M = cat(2, M, simData.("agent"+indAgent).controlOutput(indDim).valueDiscrete);
	end
end
data.controlOutput = export.dd(...
    'M', M, ...
    'header', header, ...
    'filename', 'controlOutput', ...
    'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
    );

%Export norm ob the observer errors of S and P matrices
header = {'t'};
M = [tDomain.grid];
for indAgent = 1:mas.N
	eS = norm(S - simData.("agent"+indAgent).observer.S);
	eP = norm(P - simData.("agent"+indAgent).observer.P);
	eF = norm(Sr - simData.("agent"+indAgent).observer.Sr);
	eR = norm(Pr - simData.("agent"+indAgent).observer.Pr);
	header = cat(2, header, "eS"+indAgent);
	header = cat(2, header, "eP"+indAgent);
	header = cat(2, header, "eF"+indAgent);
	header = cat(2, header, "eR"+indAgent);
	M = cat(2, M, eS.valueDiscrete);
	M = cat(2, M, eP.valueDiscrete);
	M = cat(2, M, eF.valueDiscrete);
	M = cat(2, M, eR.valueDiscrete);
end
data.obsErrMatr = export.dd(...
    'M', M, ...
    'header', header, ...
    'filename', 'obsErrorSignalMatrices', ...
    'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
    );

%Export observed v and w
header = {'t', 'w1', 'w2'};
M = [tDomain.grid, simData.agent1.observer.w(1).valueDiscrete, simData.agent1.observer.w(2).valueDiscrete];
for indAgent = 1:mas.N
	for indDim = 1:size(Sr, 1)
		header = cat(2, header, ("w"+indAgent)+indDim);
		M = cat(2, M, simData.("agent"+indAgent).observer.w(indDim).valueDiscrete);
	end
	for indDim = 1:size(S, 1)
		header = cat(2, header, ("v"+indAgent)+indDim);
		M = cat(2, M, simData.("agent"+indAgent).observer.v(indDim).valueDiscrete);
	end
end
data.observedSignals = export.dd(...
    'M', M, ...
    'header', header, ...
    'filename', 'observedSignals', ...
    'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
    );


export.Data.exportAll(data);