function simData = simulateInternalModel(simulationSetting, multiAgent, internalModelController, referenceObserver, NameValue)

%SIMULATESIMPLE Runs a simulation in a for-loop
%   Detailed explanation goes here
arguments
    simulationSetting;
    multiAgent;
	internalModelController cell;
    referenceObserver ReferenceObserver;
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
tDomain = quantity.Domain("t", t);

%Check if c is large enough
if(NameValue.c < 1/(2*min(real(eig(referenceObserver.H)))))
	warning("c is too low setting it to" + 1/(2*min(real(eig(multiAgent.H)))));
	NameValue.c = 1/(2*min(real(eig(multiAgent.H))));
end
% Handle exogenous signals
exSig = struct(NameValue.exogenousInput{:});

%Handle initial conditions for v and S
v = zeros(multiAgent.N*referenceObserver.rDimensions*internalModelController{1}.n_w, 1);
S = zeros(internalModelController{1}.n_w, internalModelController{1}.n_w,multiAgent.N-1);
for idxAgent = 1:multiAgent.N
	for idx = 1 : 2 : length(NameValue.initialCondition)
		if NameValue.initialCondition{idx}.contains("agent"+idxAgent)
			if NameValue.initialCondition{idx}.contains(".v")
				v((idxAgent-1)*referenceObserver.rDimensions*internalModelController{1}.n_w+1 : idxAgent*referenceObserver.rDimensions*internalModelController{1}.n_w) = NameValue.initialCondition{idx+1};
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
    for idxAgent = 1:length(internalModelController)
        for idxb = 1:length(internalModelController{idxAgent}.b)
            C0 = ctrb(exSig.S.atIndex(ctrbVec(idx)), internalModelController{idxAgent}.b{idxb});
            if length(exSig.S) ~= rank(C0)
                warning("The pair (S,b_"+idxb+") is not controllable for at least one time step");
            end
        end
    end
end

%Run the reference observer
if ~isfield(exSig, "w")
   exSig.w = quantity.Discrete.zeros([referenceObserver.nDimensions, 1], tDomain); 
end
% F = exSig.S(end-referenceObserver.nDimensions+1:end, end-referenceObserver.nDimensions+1:end);
if ~isfield(exSig, "R")
   exSig.Pr = eye(referenceObserver.nDimensions);
end
simDataRef = referenceObserver.simulate(exSig.w, exSig.F, exSig.R);
rObs = [];
for indAgent = 0 : multiAgent.N-1
    rObs = cat(1, rObs, simDataRef.("agent" + indAgent).P*simDataRef.("agent" + indAgent).w);
end
rObs = rObs.on(tDomain.grid).';

%Build the blkdiag of the B of all controllers
BAgg = {};
for indAgent = 1 : multiAgent.N
    BAgg{end+1} = blkdiag(internalModelController{indAgent}.b{:});
end
BAgg = blkdiag(BAgg{:});

%Get discrete state space representation
plantStateSpace = multiAgent.network.getSimulationModel(simulationSetting);
plantStateSpace = stateSpace.addInput2Output(plantStateSpace);

% stitch initial condition together
x = multiAgent.network.combineInitialCondition(plantStateSpace, NameValue.initialCondition);

% stitch exogenous signals together
u = stateSpace.combineInputSignals(plantStateSpace, t, NameValue.exogenousInput{:});
ud = u.on(t);

%Compute control input parameters
for idx = 1:multiAgent.N
    multiAgent.singleSys{idx}.setBacksteppingKernel();
    uStab{idx} = multiAgent.singleSys{idx}.getStabilizingStateFeedback();
    uStab{idx}.setApproximation(multiAgent.network.agent(idx).grid);
    uStabCop{idx} = uStab{idx}.Cop;
end
uStabDisc = blkdiag(uStabCop{:});
xIcOut = plantStateSpace.C*x;
SObs = internalModelController{1}.reconstructInternalModel(exSig.S, S);
%Take the estimates of S and construct an aggreagted p-copy of them
S = kron(eye(referenceObserver.rDimensions), SObs.agent0);
for idx = 1:multiAgent.N-1
	S = blkdiag(S, kron(eye(referenceObserver.rDimensions), SObs.("agent"+idx)));
end
%Compute control input for first time step
K = [];
Eps = [];
grid = {};
for idx = 1:multiAgent.N
	kernel = multiAgent.singleSys{idx}.backsteppingKernel.value;
	[K_agent, Eps_agent] = internalModelController{idx}.getControlGains(SObs.("agent" + (idx-1)).atIndex(1), NameValue.Q, NameValue.R, NameValue.c);
	K = blkdiag(K, K_agent);
	if isempty(Eps)
		Eps = Eps_agent;
	else
		Eps = blkdiag(Eps, Eps_agent);
	end
	K_cx{idx} = subs(Eps_agent - int(Eps_agent.subs("z", "zeta")*kernel.subs(["z", "zeta"], "zeta", "z"), "zeta", "z", 1), "z", "zeta");
	grid = cat(1, grid, multiAgent.network.agent(idx).grid);
end
uState = model.Output("control", "C", blkdiag(K_cx{:}));
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
				[K_agent, Eps_agent] = internalModelController{indAgent}.getControlGains(SObsAgent, NameValue.Q, NameValue.R, NameValue.c);
				K((indAgent-1)*size(K_agent, 1)+1:indAgent*size(K_agent, 1), (indAgent-1)*size(K_agent, 2)+1:indAgent*size(K_agent, 2)) = K_agent;
				K_cx{indAgent} = subs(Eps_agent - int(Eps_agent.subs("z", "zeta")*kernel.subs("z", "zeta"), "zeta", "z", 1), "z", "zeta");
			end
		end
	else
        adaptationDelay = adaptationDelay - 1;
    end
    %Adapt and update internal model and feedback gains
	v = cat(2, v, v(:, end) + tStep*(SObs.aggregated.atIndex(i)*v(:, end)...
        + BAgg*(simDataNy(contains(plantStateSpace.OutputName, "controlOutput"), end) - rObs(:,i))));
	uState.setApproximation(grid);
	uStateDisc = uState.Cop;
    x = plantStateSpace.A*x + plantStateSpace.B*ud(i,:).';
	xOut = plantStateSpace.C*x + plantStateSpace.D*ud(i,:).';
    %Compute control input for this time step
	u = ud(i+1);
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
	simData.("agent"+iAgent).r = simDataRef.("agent" + (indAgent-1)).P*simDataRef.("agent" + (indAgent-1)).w;
end
end

