tDomain = quantity.Domain("t", linspace(0, 10, 501));
zDomain = quantity.Domain("z", linspace(0, 1, 21));
exportData = false;

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 2 0 0.5; 0 0.5 0.5 0];
splitDomain = tDomain.split(15);
S = quantity.Piecewise({[0 0 0; 0 0 pi/2; 0 -pi/2 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 0 0; 0 0 pi/4; 0 -pi/4 0] * quantity.Discrete.ones(1, splitDomain(2))});
S = quantity.Discrete(S);
% S = quantity.Discrete.zeros(1, tDomain);

syms z t;
discOne = quantity.Discrete.ones(1, zDomain);

N = 3;

if N == 3
	Lambda = quantity.Symbolic(diag([3+z, 1.5+z/2, -2+0.5*z]), zDomain, "name", "Lambda");
	A = quantity.Symbolic([0, -0.5, 0.25*exp(z);...
		-0.25*exp(z), 0, 0.2;...
		0.2, 0.25*exp(z), 0], zDomain, "name", "A");
	Q0 = [-0.4, 0.2];
	Q1 = [0.1; -0.3];
	
	G1 = [0; discOne/4; 0];
	G2 = [.6];
	G3 = [-0.5; 0.5];
	G4 = [0; -.4];
	output = model.Output("controlOutput", "C", [.8, -.2, .6; 0, 0, 0]*quantity.Discrete(zDomain.grid, zDomain)^2, "C0", [0.8, 0.4, 0; 0, .6, 0], "C1", [0, 0, .5; 0, 0, 1]);

	dLambda = quantity.Symbolic(diag([0, 0.5, -0.5]), zDomain);
	dA = quantity.Symbolic([0, -0.2, 0.1;...
		-0.25, 0, -0.2;...
		0, 0.2, 0], zDomain);
	dQ0 = [-0.1, -0.2];
	dQ1 = [0; 0.1];
elseif N == 2
	Lambda = quantity.Symbolic(diag([2+z/2, -3+z/4]), zDomain, "name", "Lambda");
	A = quantity.Symbolic([0, z; 1-z, 0], zDomain, "name", "A")*.5;
	Q0 = .5;
	Q1 = .5;
	
	G1 = [0; 1]*quantity.Symbolic(1-z/2, zDomain, "name", "G1");
	G2 = [0.3];
	G3 = [0.8];
	G4 = [-0.4];
	output = model.Output("controlOutput", "C", ([.5, -.2]*quantity.Discrete(zDomain.grid, zDomain)^2 + [.2, .4]), "C0", [1, 0], "C1", [0, 1]);
end

% Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
% A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
% Q0 = 0;
% Q1 = 0;
% 
% G1 = [0; 0]*quantity.Symbolic(1+z, zDomain, "name", "G1");
% G2 = [0];
% G3 = [0];
% G4 = [0];

% output = model.Output("controlOutput", "C0", [0, 1]);
% output = model.Output("controlOutput", "C", [0 1]*quantity.Discrete(zDomain.grid, zDomain));

% b = {[1;1;1]};
% b = {1};
mu = 4;
L = 4*eye(length(S));
if N == 3
	b = {ones(length(S), 1)*2, ones(length(S), 1)*2};
	reference = quantity.Piecewise({quantity.Discrete([2*sin(pi/2*splitDomain(1).grid), sin(pi/2*splitDomain(1).grid+pi/4)], splitDomain(1)), ...
    	quantity.Discrete([sin(pi/4*splitDomain(2).grid), 2*sin(pi/4*splitDomain(2).grid+pi/4)], splitDomain(2))});
	reference = quantity.Discrete(reference);
	icAgent1 = quantity.Symbolic([0; 0; 0], zDomain);
	icV1 = [0; 0; 0; 0; 0; 0];
	ic = {"agent1.x", icAgent1,...
    	"agent2.x", quantity.Symbolic([0; 0; 0], zDomain),...
    	"agent3.x", quantity.Symbolic([0; 0; 0], zDomain),...
    	"agent4.x", quantity.Symbolic([0; 0; 0], zDomain),...
    	"agent1.S", S.at(0),...
    	"agent2.S", [0 0 0; 0 0 0; 0 0 0],...
    	"agent3.S", [0 0 0; 0 0 0; 0 0 0],...
    	"agent4.S", [0 0 0; 0 0 0; 0 0 0],...
    	"agent1.v", icV1,...
    	"agent2.v", [0; 0; 0; 0; 0; 0],...
    	"agent3.v", [0; 0; 0; 0; 0; 0],...
		"agent4.v", [0; 0; 0; 0; 0; 0]};
elseif N == 2
	b = {ones(length(S), 1)*2};
	reference = quantity.Piecewise({quantity.Discrete([sin(pi*splitDomain(1).grid)], splitDomain(1)), ...
    	quantity.Discrete([sin(2*pi*splitDomain(2).grid)], splitDomain(2))});
	reference = quantity.Discrete(reference);
	icAgent1 = quantity.Symbolic([0; 0], zDomain);
	icV1 = [0; 0; 1];
	ic = {"agent1.x", icAgent1,...
    	"agent2.x", quantity.Symbolic([0; 2], zDomain),...
    	"agent3.x", quantity.Symbolic([1; 0], zDomain),...
    	"agent4.x", quantity.Symbolic([2; 1], zDomain),...
    	"agent2.S", [0 0 0; 0 0 pi/2; 0 -pi/2 0]*0,...
    	"agent3.S", [0 0 0; 0 0 pi; 0 -pi 0],...
    	"agent4.S", [0 0 0; 0 0 pi/2; 0 -pi/2 0]*0,...
    	"agent1.v", icV1,...
    	"agent2.v", [0; 0; 0],...
    	"agent3.v", [0; 2; 0],...
		"agent4.v", [0; 0; 0]};
end
imc = InternalModelController(adjacencyMatrix, Lambda+dLambda, b, mu, L, "A", A+dA, "output", output, "Q0", Q0+dQ0, "Q1", Q1+dQ1, "nCoef", 6);

% SObs = imc.reconstructInternalModel(S);

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

simulationSetting = {'t', tDomain.grid};
disturbance = quantity.Piecewise({quantity.Discrete.ones(1, splitDomain(1)), ...
    quantity.Discrete.ones(1, splitDomain(2))});
disturbance = quantity.Discrete(disturbance);
disturbance = [disturbance; disturbance*-1; disturbance/2; disturbance/4];
exogenousInput = {"disturbance", disturbance, "S", S, "reference", reference};


c = 3;
R = 2;
Q = eye(length(S)*output.lengthOutput)*2;
if exportData
	adaptationDelay = 0;
else
	adaptationDelay = 10;
end
simData = tool.simulateInternalModel(simulationSetting, mas, imc, "initialCondition", ic, "exogenousInput", exogenousInput, "c", c, "R", R, "Q", Q, "adaptationDelay", adaptationDelay);

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
basepath = "C:\Users\xzb84\Documents\Überblick\InternesModelReglerAdaptiv\data";
if exportData
	% Export control output with disturbance and reference signal
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
	end
	data.obsErrMatr = export.dd(...
	    'M', M, ...
	    'header', header, ...
	    'filename', 'obsErrorSignalMatrices', ...
	    'basepath', basepath ...
	    );
	
	% export the agent states
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
	% 
	% export the reference signal
	data.reference = export.dd(...
	            'M', [tDomain.grid, reference(1).valueDiscrete, reference(2).valueDiscrete], ...
	            'header', {'t', 'r1', 'r2'}, ...
	            'filename', 'reference', ...
	            'basepath', basepath ...
	            );
	% 
	export.Data.exportAll(data);
end