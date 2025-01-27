function simData = simulateInternalModel(simulationSetting, multiAgent, internalModelController, NameValue)

%SIMULATESIMPLE Runs a simulation in a for-loop
%   Detailed explanation goes here
arguments
    simulationSetting;
    multiAgent;
	internalModelController InternalModelController;
    NameValue.exogenousInput = {};
    NameValue.initialCondition = {};
    NameValue.controller (1,1) logical = false;
    NameValue.adaptationDelay = 0;
	NameValue.Q = 1;
	NameValue.R = 1;
	NameValue.c = 1;
end % arguments
simulationSettingStruct = struct(simulationSetting{:});
if isfield(simulationSettingStruct, "t")
    t = simulationSettingStruct.t;
    tStep = t(2)-t(1);
else
    t = [];
end

%Check if c is large enough
if(NameValue.c < 1/(2*min(real(eig(multiAgent.H)))))
	warning("c is too low setting it to" + 1/(2*min(real(eig(multiAgent.H)))));
	NameValue.c = 1/(2*min(real(eig(multiAgent.H))));
end
% Handle exogenous signals
exSig = struct(NameValue.exogenousInput{:});

%Handle initial conditions for v and S
if (isa(multiAgent, "MultiAgent"))
	lengthOutput = multiAgent.output.lengthOutput;
else
	lengthOutput = multiAgent.singleSys{1}.output.lengthOutput;
end
v = zeros(multiAgent.N*lengthOutput*internalModelController.n_w, 1);
S = zeros(internalModelController.n_w, internalModelController.n_w,multiAgent.N-1);
for idxAgent = 1:multiAgent.N
	for idx = 1 : 2 : length(NameValue.initialCondition)
		if NameValue.initialCondition{idx}.contains("agent"+idxAgent)
			if NameValue.initialCondition{idx}.contains(".v")
				v((idxAgent-1)*lengthOutput*internalModelController.n_w+1 : idxAgent*lengthOutput*internalModelController.n_w) = NameValue.initialCondition{idx+1};
			end
			if NameValue.initialCondition{idx}.contains(".S") && idxAgent ~= 1
				S(:,:,idxAgent-1) = NameValue.initialCondition{idx+1};
			end
		end
	end
end

Svec = exSig.S(:);
dSvec = diff(Svec);
ctrbVec = [1; find(max(abs(dSvec.on(t)),[], 2))];
for idx = 1:2:length(ctrbVec)
	for idxb = 1:length(internalModelController.b)
		C0 = ctrb(exSig.S.atIndex(ctrbVec(idx)), internalModelController.b{idxb});
		if length(exSig.S) ~= rank(C0)
			warning("The pair (S,b_"+idxb+") is not controllable for at least one time step");
		end
	end
end

%Get discrete state space representation
plantStateSpace = multiAgent.network.getSimulationModel(simulationSetting);
plantStateSpace = stateSpace.addInput2Output(plantStateSpace);

% stitch initial condition together
x = multiAgent.network.combineInitialCondition(plantStateSpace, NameValue.initialCondition);

% stitch exogenous signals together
u = stateSpace.combineInputSignals(plantStateSpace, t, NameValue.exogenousInput{:});
ud = u.on(t);
r = exSig.reference.on(t);

%Compute control input parameters
uStab = internalModelController.getStabilizingStateFeedback();
% uStabDisc = kron(eye(multiAgent.N), uStab.Cop);
uStabDisc = [];
for idx = 1:multiAgent.N
	uStab.setApproximation(multiAgent.network.agent(idx).grid);
	uStabDisc = blkdiag(uStabDisc, uStab.Cop);
end
xIcOut = plantStateSpace.C*x;
% xIcOut = kron(eye(multiAgent.N), multiAgent.network.agent(1).approximationConvectionShift.') * x;
SObs = internalModelController.reconstructInternalModel(exSig.S, S);
%Take the estimates of S and construct an aggreagted p-copy of them
S = kron(eye(lengthOutput), SObs.agent0);
for idx = 1:multiAgent.N-1
	S = blkdiag(S, kron(eye(lengthOutput), SObs.("agent"+idx)));
end
%Compute control input for first time step
K = [];
Eps = [];
kernel = internalModelController.backsteppingKernel.value;
grid = {};
for idx = 1:multiAgent.N
	[K_agent, Eps_agent] = internalModelController.getControlGains(SObs.("agent" + (idx-1)).atIndex(1), NameValue.Q, NameValue.R, NameValue.c);
	K_agent_out{idx} = K_agent;
	K = blkdiag(K, K_agent);
	if isempty(Eps)
		Eps = Eps_agent;
	else
		Eps = blkdiag(Eps, Eps_agent);
	end
	K_cx{idx} = subs(Eps_agent - int(Eps_agent.subs("z", "zeta")*kernel.subs(["z", "zeta"], "zeta", "z"), "zeta", "z", 1), "z", "zeta");
	grid = cat(1, grid, multiAgent.network.agent(idx).grid);
end
uState = model.Output("control", "C", kron(blkdiag(1,multiAgent.H), eye(size(K_cx{1}, 1)))*blkdiag(K_cx{:}));
uState.setApproximation(grid);
uStateDisc = uState.Cop;
u = ud(1);
stateNameVector = plantStateSpace.OutputName(contains(plantStateSpace.OutputName, "x"));
IdxOffset = 0;
uStateDisc1 = [];
uStabDisc1 = [];
for idx = 1:multiAgent.N
	uStateAgent = uStateDisc(:, contains(stateNameVector, "agent"+idx+".x"));
	uStabAgent = uStabDisc(:, contains(stateNameVector, "agent"+idx+".x"));
% 	B = multiAgent.network.agent(idx).approximationConvectionShift * plantStateSpace.D(contains(plantStateSpace.OutputName, "agent"+idx+".x"), :);
% 	BcPosi = find(B(:,contains(plantStateSpace.InputName, "agent"+idx+".control")))+(idx-1)*2;
% 	BcPosi = find(plantStateSpace.B(:,find(contains(plantStateSpace.InputName, "agent"+idx+".control"))))+(idx-1)*2;
% 	uStateDisc1 = cat(2, uStateDisc1, [uStateAgent(:,1:BcPosi-IdxOffset-1),uStateAgent(:,BcPosi-IdxOffset+2:end)]);
% 	uStabDisc1 = cat(2, uStabDisc1, [uStabAgent(:,1:BcPosi-IdxOffset-1),uStabAgent(:,BcPosi-IdxOffset+2:end)]);
	uStateDisc1 =  cat(2, uStateDisc1, uStateAgent * multiAgent.network.agent(idx).approximationConvectionShift.');
	uStabDisc1 = cat(2, uStabDisc1, uStabAgent * multiAgent.network.agent(idx).approximationConvectionShift.');
	IdxOffset = IdxOffset + size(uStateAgent, 2);
end
controlInput = uStabDisc1 * x - K*(v - uStateDisc1 * x);
ud(1, contains(plantStateSpace.InputName, "control")) = controlInput;
y = plantStateSpace.C*x + plantStateSpace.D*ud(1,:).';
controlInput = uStabDisc * y(contains(plantStateSpace.OutputName, "x")) - K*(v(:, end) - uStateDisc * y(contains(plantStateSpace.OutputName, "x")));
ud(1, contains(plantStateSpace.InputName, "control")) = controlInput;
simDataNy = plantStateSpace.C*x + plantStateSpace.D*ud(1,:).';
% simDataNy = zeros(length(plantStateSpace), 1);
% simDataNy(contains(plantStateSpace.OutputName, "x")) = kron(eye(multiAgent.N), multiAgent.network.agent(1).approximationConvectionShift.')*x;
adaptationDelay = 0;
averageIterationTime = 0;
runningTime = 0;
h = waitbar(0, "Simulation in progress");
for i = 1:length(t)-1
	waitbar(i/(length(t)-1), h, "Simulation in progress, ETA: " + ceil(((length(t)-1-i)*averageIterationTime)) + "s");
	startingTime = tic;
    if adaptationDelay == 0
        adaptationDelay = NameValue.adaptationDelay;
		for indAgent = 1 : multiAgent.N
			if i == 1 ||(i > 1 && any(any(SObs.("agent" + (indAgent-1)).atIndex(i) ~= SObs.("agent" + (indAgent-1)).atIndex(i-1))))
				SObsAgent = SObs.("agent" + (indAgent-1)).atIndex(i);
				[K_agent, Eps_agent] = internalModelController.getControlGains(SObsAgent, NameValue.Q, NameValue.R, NameValue.c);
				K_agent_out{indAgent} = cat(3, K_agent_out{indAgent}, K_agent);
				K((indAgent-1)*size(K_agent, 1)+1:indAgent*size(K_agent, 1), (indAgent-1)*size(K_agent, 2)+1:indAgent*size(K_agent, 2)) = K_agent;
				K_cx{indAgent} = subs(Eps_agent - int(Eps_agent.subs("z", "zeta")*kernel.subs("z", "zeta"), "zeta", "z", 1), "z", "zeta");
			else
				K_agent_out_i = K_agent_out{indAgent};
				K_agent_out{indAgent} = cat(3, K_agent_out{indAgent}, K_agent_out_i(:, :, end));
			end
		end
	else
        adaptationDelay = adaptationDelay - 1;
		for idx = 1:multiAgent.N
			K_agent_out_i = K_agent_out{idx};
			K_agent_out{idx} = cat(3, K_agent_out{idx}, K_agent_out_i(:, :, end));
		end
    end
    %Adapt and update internal model and feedback gains
	v = cat(2, v, v(:, end) + tStep*(SObs.aggregated.atIndex(i)*v(:, end) + ...
		kron(blkdiag(1,internalModelController.H), internalModelController.B)*simDataNy(contains(plantStateSpace.OutputName, "controlOutput"), end)...
		- kron([1;internalModelController.adjacencyMatrix(2:end,1)], internalModelController.B) * r(i, :).'));
	uState = model.Output("control", "C", kron(blkdiag(1,multiAgent.H), eye(size(K_cx{1}, 1)))*blkdiag(K_cx{:}));
	uState.setApproximation(grid);
	uStateDisc = uState.Cop;
    x = plantStateSpace.A*x + plantStateSpace.B*ud(i,:).';
	xOut = plantStateSpace.C*x + plantStateSpace.D*ud(i,:).';
    %Compute control input for this time step
	u = ud(i+1,:);
	u(contains(plantStateSpace.InputName, "control")) = uStabDisc * xOut(contains(plantStateSpace.OutputName, "x")) - K*(v(:, end)...
		- uStateDisc * simDataNy(contains(plantStateSpace.OutputName, "x"), end));
	y = plantStateSpace.C*x + plantStateSpace.D*u.';
	controlInput = uStabDisc * y(contains(plantStateSpace.OutputName, "x")) - K*(v(:, end) - uStateDisc * y(contains(plantStateSpace.OutputName, "x")));
    ud(i+1, contains(plantStateSpace.InputName, "control")) = controlInput;
    %Update the state variable and compute the output using the discretized system
    simDataNy = cat(2, simDataNy, plantStateSpace.C*x + plantStateSpace.D*ud(i+1,:).');
	closingTime = toc(startingTime);
	runningTime = runningTime + closingTime;
	averageIterationTime = runningTime / i;
end
try
	close(h);
catch
	%Do nothing, the waitbar was closed manually
end
%Format raw data into function output
simOutput = simDataNy.';
outputNames = plantStateSpace.OutputName;
grid = multiAgent.network.getSpatialGridsOfSimulation();
simData = stateSpace.simulationOutput2Quantity(...
	simOutput, t, outputNames, "z", grid);
simData.predefinedInput = u;
vObs = quantity.Discrete(permute(v(:,1:end), [2 1]), quantity.Domain("t", t));
for iAgent = 1:multiAgent.N
    simData.("agent"+iAgent).S = SObs.("agent" + (iAgent-1));
    simData.("agent"+iAgent).v = vObs((iAgent-1)*size(SObs.agent0,2)+1:iAgent*size(SObs.agent0,2));
	simData.("agent"+iAgent).K = K_agent_out{iAgent};
end
end

