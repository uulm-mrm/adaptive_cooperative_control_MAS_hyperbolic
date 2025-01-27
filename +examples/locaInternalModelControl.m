tDomain = quantity.Domain("t", linspace(0, 45, 2001));
zDomain = quantity.Domain("z", linspace(0, 1, 21));
exportData = true;

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 2 0 0.5; 0 0.5 0.5 0];
splitDomain = tDomain.split(20);
S = quantity.Piecewise({[0 0 0; 0 0 pi/2; 0 -pi/2 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 0 0; 0 0 pi/4; 0 -pi/4 0] * quantity.Discrete.ones(1, splitDomain(2))});
S = quantity.Discrete(S);

syms z t;
imc = {4};
agents = {4};

controlInput = model.Input("control", "B1", eye(2));

%Agent1
Lambda = quantity.Symbolic(diag([3+z, 2-z/2, -(1.5+exp(-z))]), zDomain, "name", "Lambda");
% Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0, 0; exp(-z/4)/3, 0, 0; 0 -exp(-z/4)/3, 0], zDomain, "name", "A");
Q0 = [0.4, 0.2];
Q1 = [0.2; 0.3];

dLambda = quantity.Symbolic(diag([0.4, -0.4, 0]), zDomain);
dA = quantity.Symbolic([0, 0, 0.1; -0.2, 0, 0; 0.1 0, 0], zDomain);
dQ0 = [-0.1, 0];
dQ1 = [0; 0.1];

G1 = [0; 0.2; -0.1]*quantity.Symbolic(1+z, zDomain, "name", "G1");
G2 = [0];
G3 = [0.5; 0];
G4 = [0; 0.3];

disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
C = quantity.Discrete([ones(floor(length(zDomain.grid)/2), 1); zeros(floor(1*length(zDomain.grid)/2)+1, 1)]/2, zDomain);
output = model.Output("controlOutput", "C", [C, -C, 0; 0, 0, C/2]*0+[1, 1, 0; -1, 0, 1]/2*0, "C1", [0, 0, 1; 0, 0, 0]/2, "C0", [0.2, 0.4, 0; -0.3, 0.6, 0]);

b = {ones(length(S), 1)*2, ones(length(S), 1)*2};
mu = 4;
L = 4*eye(length(S));

agents{1} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(1));
imc{1} = InternalModelController(adjacencyMatrix, Lambda+dLambda, b, mu, L, "A", A+dA, "output", output, "Q0", Q0+dQ0, "Q1", Q1+dQ1, "nCoef", 6);

%Agent2
Lambda = quantity.Symbolic(diag([2+z, 1+z/2, -(1.5+exp(-z)), -3-z/4]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0, 0, -exp(z/2)/4; exp(-z/4)/5, 0, 0, 0; 0 -exp(-z/4)/5, 0, 0; exp(z/2)/4, 0, 0, 0], zDomain, "name", "A");
Q0 = [0.4, 0.2; 0, -0.2];
Q1 = [0.2, -0.4; 0.2, 0.1];

dLambda = quantity.Symbolic(diag([0, -z/4, 0.4, -1]), zDomain);
dA = quantity.Symbolic([0, 0, 0.5, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, -0.5, 0, 0], zDomain);
dQ0 = [0, 0; 0.4, 0];
dQ1 = [-0.4, 0; 0, 0];

G1 = [0; 0.2; -0.1; 0]*quantity.Symbolic(1+z, zDomain, "name", "G1");
G2 = [0.5; 0];
G3 = [0; 0.3];
G4 = [-0.2; 0.3];

disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
C = quantity.Discrete([ones(floor(length(zDomain.grid)/2), 1); zeros(floor(1*length(zDomain.grid)/2)+1, 1)]/2, zDomain);
output = model.Output("controlOutput", "C", [C, -C, 0, 0; 0, 0, C/2, 0]*0+[1, 0, 0, 0; 0, 0, 1, 0.5]/2, "C1", [0, 0, 1, 0.5; 0, 0, 0, -1]/2, "C0", [0.2, 0.4, 0, 0; -0.3, 0.6, 0, 0]);

b = {ones(length(S), 1)*2, ones(length(S), 1)*2};
mu = 4;
L = 4*eye(3);

agents{2} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(2));
imc{2} = InternalModelController(adjacencyMatrix, Lambda+dLambda, b, mu, L, "A", A+dA, "output", output, "Q0", Q0+dQ0, "Q1", Q1+dQ1, "nCoef", 6);

%Agent3
Lambda = quantity.Symbolic(diag([3-z/2, 1+z/4, -1-exp(-z/2)]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, z/4, 0; exp(-z/4)/4, 0, 0; 0 -exp(-z/4)/4, 0], zDomain, "name", "A");
Q0 = [0, 0.3];
Q1 = [0.2; -0.6];

dLambda = quantity.Symbolic(diag([0, -0.2, 0.6]), zDomain);
dA = quantity.Symbolic([0, -0.2, 0; 0, 0, 0; 0 0.4, 0], zDomain);
dQ0 = [0, 0];
dQ1 = [0.3; -0.1];

G1 = [0.4; -0.2; -0.3]*quantity.Symbolic(1+z, zDomain, "name", "G1");
G2 = [-0.2];
G3 = [0.5; 0];
G4 = [-0.6; 0.3];

disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
C = quantity.Discrete([ones(floor(length(zDomain.grid)/2), 1); zeros(floor(1*length(zDomain.grid)/2)+1, 1)]/2, zDomain);
output = model.Output("controlOutput", "C", [C, -C, 0; 0, 0, C/2]*0+[1, 1, 0; -1, 0 1]/2, "C1", [0, 0, 0; 1, 0, 0]/2, "C0", [0.2, 0.4, 0; -0.3, 0.6, 0]);

b = {ones(length(S), 1)*2, ones(length(S), 1)*2};
mu = 4;
L = 4*eye(3);

agents{3} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(3));
imc{3} = InternalModelController(adjacencyMatrix, Lambda+dLambda, b, mu, L, "A", A+dA, "output", output, "Q0", Q0+dQ0, "Q1", Q1+dQ1, "nCoef", 6);

%Agent4
Lambda = quantity.Symbolic(diag([2+z, 1-exp(-z/4)*z/2, -2-z]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0.5, 0; exp(-z/2)/5, 0, 0; 0.5 -exp(-z/2)/5, 0], zDomain, "name", "A");
Q0 = [0.2, 0.6];
Q1 = [-0.5; 0.1];

dLambda = quantity.Symbolic(diag([-0.5, 0, 0.2]), zDomain);
dA = quantity.Symbolic([0, 0, 0; z/4, 0, 0; -0.6, 0, 0], zDomain);
dQ0 = [0, 0.2];
dQ1 = [0; 0.1];

G1 = [0; 0.2; 0.3]*quantity.Symbolic(1+z, zDomain, "name", "G1");
G2 = [0.3];
G3 = [0.5; -0.5];
G4 = [0; 0.3];

disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
C = quantity.Discrete([ones(floor(length(zDomain.grid)/2), 1); zeros(floor(1*length(zDomain.grid)/2)+1, 1)]/2, zDomain);
output = model.Output("controlOutput", "C", [C, -C, 0; 0, 0, C/2]*0+[1, 1, 0; -1, 0 1]/2, "C1", [0, 0, 1; 0, 0, 0]/2, "C0", [0.2, 0.4, 0; -0.3, 0.6, 0]);

b = {ones(length(S), 1), ones(length(S), 1)};
mu = 4;
L = 4*eye(3);

agents{4} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(4));
imc{4} = InternalModelController(adjacencyMatrix, Lambda+dLambda, b, mu, L, "A", A+dA, "output", output, "Q0", Q0+dQ0, "Q1", Q1+dQ1, "nCoef", 6);

mas = MultiAgentHeterogenous(agents, adjacencyMatrix);

arf = ReferenceObserver(adjacencyMatrix, 4, 4, 4, 4, 4, 2, 4, 4);

simulationSetting = {'t', tDomain.grid};
% disturbance = quantity.Piecewise({quantity.Discrete(sin(1*splitDomain(1).grid), splitDomain(1)), ...
%     quantity.Discrete(sin(2*splitDomain(2).grid), splitDomain(2))});
% disturbance = quantity.Discrete(disturbance);
% disturbance = ones(4, 1)*disturbance;
disturbance = quantity.Piecewise({quantity.Discrete.ones(1, splitDomain(1)), ...
    quantity.Discrete.ones(1, splitDomain(2))});
disturbance = quantity.Discrete(disturbance);
disturbance = [disturbance; disturbance*-1; disturbance/2; disturbance/4];
reference = quantity.Piecewise({quantity.Discrete([2*sin(pi/2*splitDomain(1).grid), sin(pi/2*splitDomain(1).grid+pi/4)], splitDomain(1)), ...
    	quantity.Discrete([sin(pi/4*splitDomain(2).grid), 2*sin(pi/4*splitDomain(2).grid+pi/4)], splitDomain(2))});
	reference = quantity.Discrete(reference);
% reference = quantity.Piecewise({quantity.Discrete(sin(1*splitDomain(1).grid), splitDomain(1)), ...
%     quantity.Discrete(sin(2*splitDomain(2).grid), splitDomain(2))})*[1;1];
% w = quantity.Discrete.ones(1, tDomain);
w = quantity.Piecewise({quantity.Discrete([2*sin(pi/2*splitDomain(1).grid), 2*cos(pi/2*splitDomain(1).grid), sin(pi/2*splitDomain(1).grid+pi/4), cos(pi/2*splitDomain(1).grid+pi/4)], splitDomain(1)), ...
    quantity.Discrete([sin(pi/4*splitDomain(2).grid), cos(pi/4*splitDomain(2).grid), 2*sin(pi/4*splitDomain(2).grid+pi/4), 2*cos(pi/4*splitDomain(2).grid+pi/4)], splitDomain(2))});
F = quantity.Piecewise({[0 pi/2 0 0; -pi/2 0 0 0; 0 0 0 pi/2; 0 0 -pi/2 0] * quantity.Discrete.ones(1, splitDomain(1)), ...
	[0 pi/4 0 0; -pi/4 0 0 0; 0 0 0 pi/4; 0 0 -pi/4 0] * quantity.Discrete.ones(1, splitDomain(2))});
F = quantity.Discrete(F);
Rm = [1 0 0 0; 0 0 1 0];
exogenousInput = {"disturbance", disturbance, "S", S, "F", F, "R", Rm, "w", w};

icAgent1 = quantity.Symbolic([0; 0; 0], zDomain);
icV1 = [0; 0; 0; 0; 0; 0];

ic = {"agent1.x", icAgent1,...
    "agent2.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
    "agent3.x", quantity.Symbolic([0; 0; 0], zDomain),...
    "agent4.x", quantity.Symbolic([0; 0; 0], zDomain),...
    "agent1.S", [0 1 0; 0 0 0; 0 0 0],...
    "agent2.S", [1 0 0; 0 0 0; 0 0 0],...
    "agent3.S", [0 1 0; 0 0 0; 0 0 0],...
    "agent4.S", [0 0 0; 0 0 0; 0 0 0],...
    "agent1.v", icV1,...
    "agent2.v", [0; 0; 0; 0; 0; 0],...
    "agent3.v", [0; 0; 0; 0; 0; 0],...
	"agent4.v", [0; 0; 0; 0; 0; 0]};
c = 3;
R = 2;
Q = eye(6)*2;
adaptationDelay = 0;
simData = tool.simulateLocalInternalModel(simulationSetting, mas, imc, arf, "initialCondition", ic, "exogenousInput", exogenousInput, "c", c, "R", R, "Q", Q, "adaptationDelay", adaptationDelay);

for idx = 1:length(adjacencyMatrix)
	simData.("agent"+idx).controlOutput.plot();
% 	simData.("agent"+idx).v.plot();
% 	simData.("agent"+idx).x.plot();
end

exogenousInput = {"disturbance", disturbance, "S", S, "reference", reference*0};
simDataDist = tool.simulateInternalModel(simulationSetting, mas, imc, "initialCondition", ic, "exogenousInput", exogenousInput, "c", c, "R", R, "Q", Q, "adaptationDelay", adaptationDelay);
exogenousInput = {"disturbance", disturbance*0, "S", S, "reference", reference};
simDataRef = tool.simulateInternalModel(simulationSetting, mas, imc, "initialCondition", ic, "exogenousInput", exogenousInput, "c", c, "R", R, "Q", Q, "adaptationDelay", adaptationDelay);

%% 
basepath = "C:\Users\xzb84\Documents\Überblick\InternesModelReglerAdaptiv\dataLocal";
if exportData
	% Export control output
	header = {'t'};
	M = [tDomain.grid];
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
    	'basepath', basepath ...
    	);
	
	% Export control output with disturbance
	header = {'t'};
	M = [tDomain.grid];
	for indAgent = 1:mas.N
		for indDim = 1:2
			header = cat(2, header, ("y"+indAgent)+indDim);
			M = cat(2, M, simDataDist.("agent"+indAgent).controlOutput(indDim).valueDiscrete);
		end
	end
	data.controlOutputDist = export.dd(...
	    'M', M, ...
	    'header', header, ...
	    'filename', 'controlOutputDist', ...
	    'basepath', basepath ...
	    );
	% Export control output with reference
	header = {'t'};
	M = [tDomain.grid];
	for indAgent = 1:mas.N
		for indDim = 1:2
			header = cat(2, header, ("y"+indAgent)+indDim);
			M = cat(2, M, simDataRef.("agent"+indAgent).controlOutput(indDim).valueDiscrete);
		end
	end
	data.controlOutputRef = export.dd(...
	    'M', M, ...
	    'header', header, ...
	    'filename', 'controlOutputRef', ...
	    'basepath', basepath ...
	    );
	
	%Export norm ob the observer errors of S 
	header = {'t'};
	M = [tDomain.grid];
	for indAgent = 1:mas.N
		eS = maxNorm(S - simData.("agent"+indAgent).S);
		header = cat(2, header, "eS"+indAgent);
		M = cat(2, M, eS.valueDiscrete);
		er = maxNorm(Rm*w - simData.("agent"+indAgent).r);
		header = cat(2, header, "er"+indAgent);
		M = cat(2, M, er.valueDiscrete);
	end
	data.obsErrMatr = export.dd(...
    	'M', M, ...
    	'header', header, ...
    	'filename', 'obsErrorSignalMatrices', ...
    	'basepath', basepath ...
    	);
	
	% % export the agent states
	tPoints = 101;
	zPoints = 7;		% reduce spatial discretization to avoid memory problems in latex
	%agent 1
	mas.network.agent(1).exportSolutionOnCharacterstic(simData, "agent1.x", ...
		"fileName", "dataA1", "path", basepath, ...
		"silent", false, ...
		"tPoints", tPoints, "zPoints", zPoints, ...
		"initialCondition", icAgent1);
	
	%agent 2
	mas.network.agent(2).exportSolutionOnCharacterstic(simData, "agent2.x", ...
		"fileName", "dataA2", "path", basepath, ...
		"silent", false, ...
		"tPoints", tPoints, "zPoints", zPoints);
	
	%agent 3
	mas.network.agent(3).exportSolutionOnCharacterstic(simData, "agent3.x", ...
		"fileName", "dataA3", "path", basepath, ...
		"silent", false, ...
		"tPoints", tPoints, "zPoints", zPoints);
	
	%agent 4
	mas.network.agent(4).exportSolutionOnCharacterstic(simData, "agent4.x", ...
		"fileName", "dataA4", "path", basepath, ...
		"silent", false, ...
		"tPoints", tPoints, "zPoints", zPoints);
	
	% export the reference signal
	data.reference = export.dd(...
            	'M', [tDomain.grid, reference(1).valueDiscrete, reference(2).valueDiscrete], ...
            	'header', {'t', 'r1', 'r2'}, ...
            	'filename', 'reference', ...
            	'basepath', basepath ...
            	);
	
	export.Data.exportAll(data);
end
