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

% G1 = quantity.Discrete([0, 0; 1, 0; 0, cos(z); 1-zDomain.grid, -0.5], zDomain, "name", "G1");
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
% Q0 = [-0.25, 0; 0, 0.4];
% Q1 = [0.3, 0.5; 0, -0.2];
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
% outputDynamic = discOne;
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

% %disturbance model simulation
% S1 = [0 1; -1 0];
% S2 = [0 2; -2 0];
% S = quantity.Piecewise({S1 * quantity.Discrete.ones(1, splitDomain(1)), S2 * quantity.Discrete.ones(1, splitDomain(2))});
% S = quantity.Discrete(S);
% % P = quantity.Piecewise({[1 0; 0 1] * quantity.Discrete.ones(1, splitDomain(1)), [2 0; 0 2] * quantity.Discrete.ones(1, splitDomain(2))});
% P = eye(2);
% distOde1 = ss(S1, zeros(size(S,1)), eye(size(S)), []);
% v1 = stateSpace.simulate(distOde1, zeros(size(S,1), length(splitDomain(1).grid)), splitDomain(1).grid, [0; 1]);
% distOde2 = ss(S2, zeros(size(S,1)), eye(size(S)), []);
% v2 = stateSpace.simulate(distOde2, zeros(size(S,1), length(splitDomain(2).grid)), splitDomain(2).grid, v1(end, :));
% v = quantity.Discrete([v1; v2(2:end, :)], tDomain);
S = 0;
distOde = ss(S, zeros(size(S,1)), eye(size(S)), []);
v = stateSpace.simulate(distOde, zeros(size(S,1), length(tDomain.grid)), tDomain.grid, 1);
v = quantity.Discrete(v, tDomain);
S = S*discOnet;
P = 1;

% %reference model simulation
% Sr1 = [0 1; -1 0];
% Sr2 = [0 2; -2 0];
% Sr = quantity.Piecewise({Sr1 * quantity.Discrete.ones(1, splitDomain(1)), Sr2 * quantity.Discrete.ones(1, splitDomain(2))});
% Sr = quantity.Discrete(Sr);
% % Pr = quantity.Piecewise({[1 0; 0 1] * quantity.Discrete.ones(1, splitDomain(1)), [2 0; 0 2] * quantity.Discrete.ones(1, splitDomain(2))});
% Pr = [1 0; 0 1];
% refOde1 = ss(Sr1, zeros(size(Sr,1)), eye(size(Sr)), []);
% w1 = stateSpace.simulate(distOde1, zeros(size(Sr,1), length(splitDomain(1).grid)), splitDomain(1).grid, [0; 1]);
% refOde2 = ss(Sr2, zeros(size(Sr,1)), eye(size(Sr)), []);
% w2 = stateSpace.simulate(distOde2, zeros(size(Sr,1), length(splitDomain(2).grid)), splitDomain(2).grid, v1(end, :));
% w = quantity.Discrete([w1; w2(2:end, :)], tDomain);
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
% simData = mas.adaptiveControlSimulation(simulationSetting, length(S), size(P, 1), arf, "initialCondition", ic, "exogenousInput", exogenousInput, "mu_S", 2, "mu_P", 2, ...
% 	"distCoef", 10, "adaptationDelayObserver", 0, "adaptationDelayController", 0, "checkControllability", true);
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

% plot(simData.agent2.controlOutput);
% plot(simData.agent4.observer.v)
% toc
% profile viewer
% mas.plotState(simData);
% simData.agent2.observer.v.plot();
% simData.agent2.observer.x.plot();
% simData.agent2.observer.S.plot();
% distObs.plotDisturbanceState(simData);

%Run old method for comparison
% simDataRef = arf.simulate(w, Sr);
% 
% for indAgent = 1:mas.N
%     S_rObs.("agent" + indAgent) = simDataRef.("agent" + (indAgent-1)).S;
%     wObs.("agent" + indAgent) = simDataRef.("agent" + (indAgent-1)).w;
% end
% outputTil = output.backstepping(mas.backsteppingKernel, "inverse", true);
% invLambda = 1/Lambda;
% phi = int(subs(invLambda, "z", "zeta"), "zeta", 0, "z");
% eigVals = [eig([0 1; 0 0]); eig([0 pi; -pi 0])];
% observable = true;
% Ns = [];
% for indEv = 1:length(eigVals)
%     psi = expm(eigVals(indEv)*(phi - phi.subs("z", "zeta")));
%     M = psi.subs("zeta", 0)*(mas.E1+mas.E2*Q0)...
%         - int(psi*invLambda.subs("z", "zeta")*mas.backsteppingKernel.getA0Target().subs("z","zeta"), "zeta", 0, "z");
%     N = outputTil.out(M);
%     Ns = [Ns, det(N)];
%     if det(N) == 0
%        observable = false; 
%     end
% end
%% 

% Export control output
% header = {'t', 'w'};
% M = [tDomain.grid, simData.agent1.observer.w(1).valueDiscrete];
% for indAgent = 1:mas.N
% 	for indDim = 1:2
% 		header = cat(2, header, ("y"+indAgent)+indDim);
% 		M = cat(2, M, simData.("agent"+indAgent).controlOutput(indDim).valueDiscrete);
% 	end
% end
% data.controlOutput = export.dd(...
%     'M', M, ...
%     'header', header, ...
%     'filename', 'controlOutput', ...
%     'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
%     );

% %Export norm ob the observer errors of S and P matrices
% header = {'t'};
% M = [tDomain.grid];
% for indAgent = 1:mas.N
% 	eS = maxNorm(S - simData.("agent"+indAgent).observer.S);
% 	eP = maxNorm(P - simData.("agent"+indAgent).observer.P);
% 	eF = maxNorm(Sr - simData.("agent"+indAgent).observer.Sr);
% 	eR = maxNorm(Pr - simData.("agent"+indAgent).observer.Pr);
% 	header = cat(2, header, "eS"+indAgent);
% 	header = cat(2, header, "eP"+indAgent);
% 	header = cat(2, header, "eF"+indAgent);
% 	header = cat(2, header, "eR"+indAgent);
% 	M = cat(2, M, eS.valueDiscrete);
% 	M = cat(2, M, eP.valueDiscrete);
% 	M = cat(2, M, eF.valueDiscrete);
% 	M = cat(2, M, eR.valueDiscrete);
% end
% data.obsErrMatr = export.dd(...
%     'M', M, ...
%     'header', header, ...
%     'filename', 'obsErrorSignalMatrices', ...
%     'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
%     );
% 
% %Export observed v and w
% header = {'t', 'w1', 'w2', 'v2', 'v3'};
% M = [tDomain.grid, simData.agent1.observer.w(1).valueDiscrete, simData.agent1.observer.w(3).valueDiscrete, v2.valueDiscrete, v3.valueDiscrete];
% for indAgent = 1:mas.N
% 	for indDim = 1:size(Sr, 1)
% 		header = cat(2, header, ("w"+indAgent)+indDim);
% 		M = cat(2, M, simData.("agent"+indAgent).observer.w(indDim).valueDiscrete);
% 	end
% 	for indDim = 1:size(S, 1)
% 		header = cat(2, header, ("v"+indAgent)+indDim);
% 		M = cat(2, M, simData.("agent"+indAgent).observer.v(indDim).valueDiscrete);
% 	end
% end
% data.observedSignals = export.dd(...
%     'M', M, ...
%     'header', header, ...
%     'filename', 'observedSignals', ...
%     'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
%     );
% 
% %Export observer errors for w and v
% header = {'t', 'ev1', 'ev2', 'ev3', 'ew1', 'ew2', 'ew3'};
% 
% ev1 = v1-simData.agent2.observer.v;
% ev2 = v2-simData.agent3.observer.v;
% ev3 = v3-simData.agent4.observer.v;
% 
% ew1 = maxNorm(w-simData.agent2.observer.w);
% ew2 = maxNorm(w-simData.agent3.observer.w);
% ew3 = maxNorm(w-simData.agent4.observer.w);
% 
% M = [tDomain.grid, ev1.valueDiscrete, ev2.valueDiscrete, ev3.valueDiscrete, ew1.valueDiscrete, ew2.valueDiscrete, ew3.valueDiscrete];
% data.obsErrors = export.dd(...
%     'M', M, ...
%     'header', header, ...
%     'filename', 'observerErrors', ...
%     'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
%     );
% 
% 
% %Export some exponential functions to show exponential decay
% eigH = eig(mas.H);
% minEigH = min(eigH(eigH~=0));
% maxEigH = max(eigH)^2;
% e1 = 8*exp(-2*minEigH*linspace(0, 5, 501));
% e2 = 8*exp(-2*minEigH*linspace(0, 5, 501));
% e3 = 80*exp(-2*minEigH*linspace(0, 5, 501));
% e4 = 8*exp(-2*minEigH*linspace(0, 5, 501));
% header = {'t', 'e1', 'e2', 'e3', 'e4'};
% M = [linspace(0, 5, 501).', e1.', e2.', e3.', e4.'];
% data.exonentials = export.dd(...
%     'M', M, ...
%     'header', header, ...
%     'filename', 'exponentials', ...
%     'basepath', "C:\Users\xzb84\Documents\tarikspapers\Cooperative output regulation for networks of hyperbolic systems using an adaptive distributed observer\autosam-2\data" ...
%     );
% 
% export.Data.exportAll(data);