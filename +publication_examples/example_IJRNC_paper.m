tDomain = quantity.Domain("t", linspace(0, 20, 8001));
zDomain = quantity.Domain("z", linspace(0, 1, 21));
exportData = true;

adjacencyMatrix = [0 0 0 0; 1.4 0 0 0.7; 0 0.9 0 0; 0 1.2 1 0];
splitDomain = tDomain.split(60);
S0 = 0;
S1 = [0 pi/2.5; -pi/2.5, 0];
S2 = [0 pi/5; -pi/5, 0];
S3 = [0 pi/7; -pi/7, 0];
S = blkdiag(S0, S1, S2);
charpol = charpoly(S);
S = [[zeros(1, length(S)-1); eye(length(S)-1)], -flip(charpol(2:end)).'];
S = quantity.Discrete.ones(1, tDomain)*S;

syms z t;
discOne = quantity.Discrete.ones(1, zDomain);

Lambda = quantity.Symbolic(diag([3+z, 1.5+z/2, -2+0.5*z, -3+2*z/3]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, -0.5, 0.25*exp(z), 0;...
	0, 0, 0.2, 0.5*exp(z);...
	0.2, 0.25*exp(z), 0, 0; ...
	0.5*exp(z/2), -0.2*exp(z), 0.5, 0], zDomain, "name", "A");
Q0 = [-0.4, 0.2; -0.2, 0.6];
Q1 = [0.1, 0.4; -0.3, -0.8];

G1 = [0; discOne/4; 0; -discOne/2];
G2 = [.6; -0.4];
G3 = [-0.5; 0.5];
G4 = [0; -.4];

C = [.8, -.2, .6, 0; 0, 0, 0, 0]*quantity.Discrete(zDomain.grid, zDomain)^2;
C0 = [0.8, 0.4, 0, 0; 0, .6, 0, 0.4];
C1 = [0, 0, .5, 0.5; 0, 0, 1, 0];
output = model.Output("controlOutput", "C", C, "C0", C0, "C1", C1);

dLambda = quantity.Symbolic(diag([0, 0.5, -0.5, 0.2]), zDomain);
dA = quantity.Symbolic([0, -0.2, 0.1, 0.3;...
	-0.25, 0, -0.2, 0.1;...
	0, 0.2, 0, 0.1; ...
	0.2, -0.1, 0, 0], zDomain);
dQ0 = [-0.1, -0.2; 0, 0.1];
dQ1 = [0, 0; 0.1, -0.2];
dC = [-.2, 0, .1, .2; 0, -.3, 0, .2];
dC0 = [0.1, 0.1, -.2, 0; -.3, 0, -.1, 0.1];
dC1 = [0, .4, -.1, 0; 0, 0, -.3, .2];

controlInput = model.Input("control", "B1", eye(size(G2, 1)));
disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;

agents{1} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent1");
agents{2} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent2");
disturbanceInput = model.Input("disturbance", "B", G1/2,...
	"B0", G2/2, "B1", G3/2, "D", misc.Gain("disturbance", G4/2, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
output = model.Output("controlOutput", "C", C+dC, "C0", C0+dC0, "C1", C1+dC1);
agents{3} = model.Transport(Lambda+dLambda, "A", A+dA, "input", input, "output", output, "Q0", Q0+dQ0, "Q1", Q1+dQ1, "prefix", "agent3");
disturbanceInput = model.Input("disturbance", "B", G1*2,...
	"B0", G2*2, "B1", G3*2, "D", misc.Gain("disturbance", G4*2, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
output = model.Output("controlOutput", "C", C-dC, "C0", C0-dC0, "C1", C1-dC1);
agents{4} = model.Transport(Lambda-dLambda, "A", A-dA, "input", input, "output", output, "Q0", Q0-dQ0, "Q1", Q1-dQ1, "prefix", "agent4");

mas = MultiAgentHeterogenous(agents, adjacencyMatrix);


mu = 4;
L = 4*eye(length(S));

b = {ones(length(S), 1)*2, ones(length(S), 1)*2};
reference = quantity.Discrete([-2*sin(pi/2.5*splitDomain(1).grid) + 1.5*cos(pi/5*splitDomain(1).grid), 1.8*cos(pi/2.5*splitDomain(1).grid) + sin(pi/5*splitDomain(1).grid)+1], tDomain);
icAgent1 = quantity.Symbolic([0; 0; 0; 0], zDomain);
icV1 = [0; 0; 0; 0; 0; 0];
ic = {"agent1.x", icAgent1,...
	"agent2.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
	"agent3.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
	"agent4.x", quantity.Symbolic([0; 0; 0; 0], zDomain),...
	"agent1.S", S.at(0),...
	"agent2.S", [[zeros(1, length(S)-1); eye(length(S)-1)],zeros(length(S), 1)],...
	"agent3.S", [[zeros(1, length(S)-1); eye(length(S)-1)],zeros(length(S), 1)],...
	"agent4.S", [[zeros(1, length(S)-1); eye(length(S)-1)],zeros(length(S), 1)],...
	"agent1.v", zeros(length(S)*size(C, 1), 1),...
	"agent2.v", zeros(length(S)*size(C, 1), 1),...
	"agent3.v", zeros(length(S)*size(C, 1), 1),...
	"agent4.v", zeros(length(S)*size(C, 1), 1)};
ics = {"agent1.x", quantity.Symbolic([1; 1; 1; 1], zDomain),...
	"agent2.x", quantity.Symbolic([1; 1; 1; 1], zDomain),...
	"agent3.x", quantity.Symbolic([1; 1; 1; 1], zDomain),...
	"agent4.x", quantity.Symbolic([1; 1; 1; 1], zDomain)};
simulationSetting = {'t', tDomain.grid};
simDataOpen = mas.network.simulate(simulationSetting, ics);
imc = InternalModelController(adjacencyMatrix, Lambda, b, mu, L, "A", A, "output", output, "Q0", Q0, "Q1", Q1, "nCoef", 40);



disturbance = quantity.Piecewise({quantity.Discrete.ones(1, splitDomain(1)), ...
    quantity.Discrete.ones(1, splitDomain(2))});
disturbance = quantity.Discrete(disturbance);
disturbance = [disturbance; disturbance*-1; disturbance/2; disturbance/4];
exogenousInput = {"disturbance", disturbance, "S", S, "reference", reference};


c = 5;
R = 2;
Q = eye(length(S)*output.lengthOutput)*2;
if exportData
	adaptationDelay = 150;
else
	adaptationDelay = 50;
end
simData = tool.simulateInternalModel(simulationSetting, mas, imc, "initialCondition", ic, "exogenousInput", exogenousInput, "c", c, "R", R, "Q", Q, "adaptationDelay", adaptationDelay);

for idx = 1:length(adjacencyMatrix)
	simData.("agent"+idx).controlOutput.plot();
end
%% 
output = model.Output("controlOutput", "C", C, "C0", C0, "C1", C1);
mas.singleSys{1}.setBacksteppingKernel();
outputTil = output.backstepping(mas.singleSys{1}.backsteppingKernel, "inverse", true);
invLambda = 1/Lambda;
phi = int(subs(invLambda, "z", "zeta"), "zeta", 0, "z");
eigVals = eig(S.at(0));
Ns = [];
for indEv = 1:length(eigVals)
    psi = expm(eigVals(indEv)*(phi - phi.subs("z", "zeta")));
    M = psi.subs("zeta", 0)*(mas.singleSys{1}.E1+mas.singleSys{1}.E2*Q0)...
        - int(psi*invLambda.subs("z", "zeta")*mas.singleSys{1}.backsteppingKernel.getA0Target().subs("z","zeta"), "zeta", 0, "z");
    N = outputTil.out(M);
    Ns = [Ns, det(N)];
    if det(N) == 0
       observable = false; 
    end
end
%% 
basepath = "C:\Users\xzb84\Documents\tarikspapers\IJRNC cooperative robust hyperbolic\data";
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
	    'basepath', basepath, 'N', 501 ...
	    );
	% Export the internal state with disturbance and reference signal
	header = {'t'};
	M = [tDomain.grid];
	for indAgent = 1:mas.N
		for indDim = 1:3
			header = cat(2, header, ("y"+indAgent)+indDim);
			M = cat(2, M, simData.("agent"+indAgent).v(indDim).valueDiscrete);
		end
	end
	data.internalState = export.dd(...
	    'M', M, ...
	    'header', header, ...
	    'filename', 'internalState', ...
	    'basepath', basepath ...
	    );
	% Export control output with disturbance
% 	header = {'t'};
% 	M = [tDomain.grid];
% 	for indAgent = 1:mas.N
% 		for indDim = 1:2
% 			header = cat(2, header, ("y"+indAgent)+indDim);
% 			M = cat(2, M, simDataDist.("agent"+indAgent).controlOutput(indDim).valueDiscrete);
% 		end
% 	end
% 	data.controlOutputDist = export.dd(...
% 	    'M', M, ...
% 	    'header', header, ...
% 	    'filename', 'controlOutputDist', ...
% 	    'basepath', basepath ...
% 	    );
	% Export control output with reference
% 	header = {'t'};
% 	M = [tDomain.grid];
% 	for indAgent = 1:mas.N
% 		for indDim = 1:2
% 			header = cat(2, header, ("y"+indAgent)+indDim);
% 			M = cat(2, M, simDataRef.("agent"+indAgent).controlOutput(indDim).valueDiscrete);
% 		end
% 	end
% 	data.controlOutputRef = export.dd(...
% 	    'M', M, ...
% 	    'header', header, ...
% 	    'filename', 'controlOutputRef', ...
% 	    'basepath', basepath ...
% 	    );
	
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
	eigH = eig(mas.H);
	minEigH = min(eigH(eigH~=0));
	e1 = 5*exp(-4*minEigH*linspace(0, 10, 1001));
	header = {'t', 'e1'};
	M = [linspace(0, 10, 1001).', e1.'];
	data.exonentials = export.dd(...
	    'M', M, ...
	    'header', header, ...
	    'filename', 'exponentials', ...
	    'basepath', basepath ...
	    );
	
	% Export control gains
	header = {'t'};
	M = [tDomain.grid];
	for indAgent = 1:mas.N
		eK = reshape(max(max(simData.agent1.K(:, :, end) - simData.("agent"+indAgent).K(:, :, :))), [size(simData.("agent"+indAgent).K, 3), 1]);
		
		header = cat(2, header, ("K"+indAgent));
		M = cat(2, M, eK);
	end
	data.controlGains = export.dd(...
	    'M', M, ...
	    'header', header, ...
	    'filename', 'controlGains', ...
	    'basepath', basepath ...
	    );
	export.Data.exportAll(data);
end