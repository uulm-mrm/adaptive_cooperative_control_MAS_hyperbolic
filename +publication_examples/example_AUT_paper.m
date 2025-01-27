zDomain = quantity.Domain("z", linspace(0, 1, 61));
tDomain = quantity.Domain("t", linspace(0, 5, 501));
syms z t;

discOne = quantity.Discrete.ones(1, zDomain);
discOnet = quantity.Discrete.ones(1, tDomain);
discz = quantity.Symbolic(z, zDomain);
%Agent 0
Lambda = quantity.Symbolic(diag([2, 1-z/2, -3]), zDomain, "name", "Lambda");
A = quantity.Discrete.zeros([3, 3], zDomain, "name", "A");
Q0 = [0, 1];
Q1 = [0.5; 0];

G1 = [0; discOne; 0];
G2 = [1];
G3 = [0; 0.5];
G4 = [0; 0];

controlInput = model.Input("control", "B1", eye(2));
disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;

output = model.Output("controlOutput", "C0", [1 0 0; 0 0 0], "C", [0 0 0; 0 -discOne 0]);


agents{1} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(1));

%Agent 1
Lambda = quantity.Symbolic(diag([6+2*z, 4+z, -3.5-z/2, -5.75-3*z/2]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0.5, -0.3*(z-1)*exp(-z), -0.25*exp(z);...
	0.3*exp(z), 0, 0.25*exp(z), 0.5;...
	0.5, 0.25*exp(z), 0, -0.25*exp(z);...
	-0.25*exp(z), 0.3*(z-1)*exp(z), 0.5, 0], zDomain, "name", "A");
Q0 = [-0.25, 2; 0, 0.4];
Q1 = [0.3, 0.5; 2, -0.2];

G1 = [0; discOne; 0; discOne-discz];
G2 = [1; 0.5];
G3 = [0.5; 0];
G4 = [0; -1];

controlInput = model.Input("control", "B1", eye(2));
disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
outputDynamic = -sqrt(discz+1);
output = model.Output("controlOutput", "C0", [1 0 0 0; 0 0 0 0], "C", [0 0 0 0; 0 outputDynamic 0 0]);

agents{2} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(2));

%Agent 2
Lambda = quantity.Symbolic(diag([4+z, 2.5+z, -4+z/2]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, -0.5, 0.25*exp(z);...
	-0.25*exp(z), 0, 0.5;...
	0.5, 0.25*exp(z), 0], zDomain, "name", "A");
Q0 = [-0.4, 0.2];
Q1 = [0.1; -0.5];

G1 = [0; discOne; 0];
G2 = [1];
G3 = [-0.5; 0.5];
G4 = [0; -1];

controlInput = model.Input("control", "B1", eye(2));
disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
output = model.Output("controlOutput", "C0", [1 0 0; 0 0 0], "C", [0 0 0; 0 outputDynamic 0]);

agents{3} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(3));

%Agent 3
Lambda = quantity.Symbolic(diag([4+z/3, 2.5-z/2, -4.75-z]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, -0.5, -0.3*(z-1)*exp(z);...
	0.3*exp(z), 0, 0.5;...
	0.1, 0, 0], zDomain, "name", "A");
Q0 = [-0.4, 0];
Q1 = [0.25; -0.5];

G1 = [0; discOne; 0];
G2 = [1];
G3 = [-1; 0.5];
G4 = [0.5; 0];

disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
output = model.Output("controlOutput", "C0", [1 0 0; 0 0 0], "C", [0 0 0; 0 outputDynamic 0]);

agents{4} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(4));

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 2 0 0.5; 0 1 0.5 0];

mas = MultiAgentHeterogenous(agents, adjacencyMatrix, "nCoefDist", 10, "nCoefRef", 10);

%disturbance model simulation
S = 0;
distOde = ss(S, zeros(size(S,1)), eye(size(S)), []);
v = stateSpace.simulate(distOde, zeros(size(S,1), length(tDomain.grid)), tDomain.grid, 1);
v = quantity.Discrete(v, tDomain);
S = S*discOnet;
P = 1;

%reference model simulation
Sr = [0 pi 0 0; -pi 0 0 0; 0 0 0 2*pi; 0 0 -2*pi 0];
refOde = ss(Sr, zeros(size(Sr,1)), eye(size(Sr)), []);
w = stateSpace.simulate(refOde, zeros(size(Sr,1), length(tDomain.grid)), tDomain.grid, [0; 2; 0; 1]);
w = quantity.Discrete(w, tDomain);
Pr = [1, 0 0 0; 0 0 1 0];

simulationSetting = {'t', tDomain.grid};
zero = quantity.Discrete.zeros(size(v), tDomain);
ic = {"agent1.x", quantity.Symbolic([0; 0; 0], zDomain),...
    "agent2.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent3.x", quantity.Symbolic([0; 0; 0], zDomain),...
    "agent4.x", quantity.Symbolic([0; 0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent3.observer.x", quantity.Symbolic([0; 0; 0], zDomain),...
    "agent4.observer.x", quantity.Symbolic([0; 0; 0], zDomain),...
    "S", [0; 1; 1; 1],...
    "v", [0; 0; 0; 0]};

v1 = v;
v2 = -v*2;
v3 = v/2;
aggregatedDisturbance = [P*zero; P*v1; P*v2; P*v3];
aggregatedDisturbanceState = [zero; v1; v2; v3];

arf = ReferenceObserver(adjacencyMatrix, length(Sr), 2, 2, 2, 2, size(Pr, 1), 2, 2);

exogenousInput = {"disturbance", aggregatedDisturbance, "S", S, "Sr", Sr, "w", w, "Pr", Pr, "P", P, "disturbanceState", aggregatedDisturbanceState};
icOpen = {"agent1.x", quantity.Symbolic([1; 1; 1], zDomain),...
    "agent2.x", quantity.Symbolic([1; 1; 1; 1], zDomain),...
    "agent3.x", quantity.Symbolic([1; 1; 1], zDomain),...
    "agent4.x", quantity.Symbolic([1; 1; 1], zDomain)};
simDataOpen = mas.network.simulate(simulationSetting, icOpen);

% Export control output of open loop systems
header = {'t'};
M = [tDomain.grid];
for indAgent = 1:mas.N
	for indDim = 1:2
		header = cat(2, header, ("y"+indAgent)+indDim);
		M = cat(2, M, simDataOpen.("agent"+indAgent).controlOutput(indDim).valueDiscrete);
	end
end
data.controlOutput = export.dd(...
    'M', M, ...
    'header', header, ...
    'filename', 'controlOutput', ...
    'basepath', "C:\Users\xzb84\Documents\tarikspapers\Authors` reply - AUT- Adaptive disturbance observer" ...
    );
export.Data.exportAll(data);


%% 

% Export control output
header = {'t', 'w'};
M = [tDomain.grid, simData.agent1.observer.w(1).valueDiscrete];
for indAgent = 1:mas.N
	for indDim = 1:2
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
	eS = maxNorm(S - simData.("agent"+indAgent).observer.S);
	eP = maxNorm(P - simData.("agent"+indAgent).observer.P);
	eF = maxNorm(Sr - simData.("agent"+indAgent).observer.Sr);
	eR = maxNorm(Pr - simData.("agent"+indAgent).observer.Pr);
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
header = {'t', 'w1', 'w2', 'v2', 'v3'};
M = [tDomain.grid, simData.agent1.observer.w(1).valueDiscrete, simData.agent1.observer.w(3).valueDiscrete, v2.valueDiscrete, v3.valueDiscrete];
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

%Export observer errors for w and v
header = {'t', 'ev1', 'ev2', 'ev3', 'ew1', 'ew2', 'ew3'};

ev1 = v1-simData.agent2.observer.v;
ev2 = v2-simData.agent3.observer.v;
ev3 = v3-simData.agent4.observer.v;

ew1 = maxNorm(w-simData.agent2.observer.w);
ew2 = maxNorm(w-simData.agent3.observer.w);
ew3 = maxNorm(w-simData.agent4.observer.w);

M = [tDomain.grid, ev1.valueDiscrete, ev2.valueDiscrete, ev3.valueDiscrete, ew1.valueDiscrete, ew2.valueDiscrete, ew3.valueDiscrete];
data.obsErrors = export.dd(...
    'M', M, ...
    'header', header, ...
    'filename', 'observerErrors', ...
    'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
    );


%Export some exponential functions to show exponential decay
eigH = eig(mas.H);
minEigH = min(eigH(eigH~=0));
maxEigH = max(eigH)^2;
e1 = 8*exp(-2*minEigH*linspace(0, 5, 501));
e2 = 8*exp(-2*minEigH*linspace(0, 5, 501));
e3 = 80*exp(-2*minEigH*linspace(0, 5, 501));
e4 = 8*exp(-2*minEigH*linspace(0, 5, 501));
header = {'t', 'e1', 'e2', 'e3', 'e4'};
M = [linspace(0, 5, 501).', e1.', e2.', e3.', e4.'];
data.exonentials = export.dd(...
    'M', M, ...
    'header', header, ...
    'filename', 'exponentials', ...
    'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
    );

export.Data.exportAll(data);