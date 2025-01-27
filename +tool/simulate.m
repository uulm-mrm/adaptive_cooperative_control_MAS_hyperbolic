function simData = simulate(simulationSetting, multiAgent, NameValue)

%SIMULATESIMPLE Runs a simulation in a for-loop
%   This function allows to run a simulation for hyperbolic model with an
%   adaptive disturbance observer. 
%	
%		simulate(simulationSetting, multiAgent, "exogenousInput", ..., "initialCondition", ...)
%		runs the simulation with the sprecified exogenous signals and
%		inital conditions. If these optional arguments are not specified
%		they are set to zero.
%	
%		simulate(..., "disturbanceObserver", ..., "adaptationDelayObserver", ...)
%		runs the simulation with the sprecified adaptive disturbance observer of
%		type AdaptiveDisturbanceObserver. The optional argument
%		adaptationDelayObserver determines how frequently the adaption
%		process is executed. If not specified the observer is adapted in
%		every time step.
%
%		simulate(..., "referenceObserver", ..., "Pr", ...)
%		runs the simulation with the sprecified adaptive reference observer of
%		type ReferenceObserver. If Pr is not specified it is set to eye.
%
%		simulate(..., "referenceObserver", ...)
%		runs the simulation with the sprecified adaptive reference observer of
%		type ReferenceObserver.
%
%		simulate(..., "controller", true, "adaptationDelayController")
%		runs the simulation using an output feedback controller. The
%		controller uses the estimates of the observers if they are
%		specified, otherwise uses the actual values. If observers are
%		specified the parameter adaptationDelayController determines how
%		frequently the control gains are adapted.

arguments
    simulationSetting;
    multiAgent;
    NameValue.disturbanceObserver = struct();
    NameValue.referenceObserver = struct();
    NameValue.exogenousInput = {};
    NameValue.initialCondition = {};
    NameValue.controller (1,1) logical = false;
    NameValue.Pr double;
    NameValue.adaptationDelayObserver = 0;
    NameValue.adaptationDelayController = 0;
end % arguments
simulationSettingStruct = struct(simulationSetting{:});
if isfield(simulationSettingStruct, "t")
    t = simulationSettingStruct.t;
    tStep = t(2)-t(1);
else
    t = [];
end
tDomain = quantity.Domain("t", t);

% Handle exogenous signals
exSig = struct(NameValue.exogenousInput{:});

%Prepare for adaptive disturbance observer
if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver")
    measurementOutput = model.Output("measurement", "C1", multiAgent.E2.');
    multiAgent.network.addCommunication2agents(measurementOutput);
    SObs = kron(ones(multiAgent.N,1),zeros(NameValue.disturbanceObserver.n_d));
	PObs = kron(ones(multiAgent.N,1),zeros(NameValue.disturbanceObserver.n_r, NameValue.disturbanceObserver.n_d));
    vObs = zeros(size(SObs, 1), 1);
    icObs = {};
    for idx = 1:2:length(NameValue.initialCondition)
       if contains(NameValue.initialCondition{idx}, "observer.x")
          icObs{end+1} = strrep(NameValue.initialCondition{idx}, "observer.x", "x");
          icObs{end+1} = NameValue.initialCondition{idx+1};
       end
    end
    for indAgent = 1 : multiAgent.N
        stateSpaceCell{indAgent} = NameValue.disturbanceObserver.getStateSpaceApproximation(quantity.Discrete.zeros([multiAgent.n, multiAgent.m], multiAgent.domain),...
            "t", t, "agent", indAgent);
    end
    observerStateSpace = stateSpace.parallel(stateSpaceCell{:});
    observerStateSpaceB = observerStateSpace.B;
    adaptationDelayObserver = 0;
	if isfield(exSig, "P") && isa(exSig.P, 'double')
       exSig.P = exSig.P*quantity.Discrete.ones(1, tDomain);
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

%Compute control input parameters if controller is set to true
if NameValue.controller
    multiAgent.setBacksteppingKernel();
    if ~isfield(exSig, "Sr")
       exSig.Sr = 0;
    end  
    multiAgent.setRegulatorEquationRefCoefficients(size(exSig.Sr,1));
    if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver")
        multiAgent.setRegulatorEquationDistCoefficients(size(exSig.S,1));
    else
        assert(isfield(exSig, "S") && isfield(exSig, "P"), "In order to consider the disturbance state v without an adaptive disturbance observer S and P must be specified");
        multiAgent.setRegulatorEquationDistCoefficients(size(exSig.S,1));
        [PI_v, gain_v] = multiAgent.solveRegulatorEquationsDisturbance(exSig.S, exSig.P);
        Kv = kron(eye(multiAgent.N), gain_v);
    end
    uStab = multiAgent.getStabilizingStateFeedback();
    uStab.setApproximation(multiAgent.network.agent(1).grid);
    uStabDisc = kron(eye(multiAgent.N), uStab.Cop);
    xIcOut = plantStateSpace.C*x;
	if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver")
		controlInput = 0;
	else
		controlInput = uStabDisc * xIcOut(contains(plantStateSpace.OutputName, "x"));
	end
    adaptationDelayController = 0;
end

%Run simulation of reference observer of one is specified
if isa(NameValue.referenceObserver, "ReferenceObserver")
    if ~isfield(exSig, "w")
       exSig.w = quantity.Discrete.zeros([NameValue.referenceObserver.nDimensions, 1], tDomain); 
    end
    if ~isfield(exSig, "Sr")
       exSig.Sr = zeros(NameValue.referenceObserver.nDimensions);
    end
    if ~isfield(exSig, "Pr")
       exSig.Pr = eye(NameValue.referenceObserver.nDimensions);
    end
    simDataRef = NameValue.referenceObserver.simulate(exSig.w, exSig.Sr, exSig.Pr); 
    if NameValue.controller
       wAgg = [];
       for indAgent = 1 - multiAgent.agent1IsLeader : multiAgent.N - multiAgent.agent1IsLeader 
           [PI_w.("agent"+indAgent), gain_w{indAgent + multiAgent.agent1IsLeader}] = multiAgent.solveRegulatorEquationsReference(simDataRef.("agent" + indAgent).S.atIndex(1), simDataRef.("agent" + indAgent).P.atIndex(1));
           wAgg = cat(1, wAgg, simDataRef.("agent"+indAgent).w.atIndex(1));
       end
       Kw = blkdiag(gain_w{:});
       controlInput = controlInput + Kw*wAgg;
    end
else
    if ~isfield(exSig, "w")
       exSig.w = quantity.Discrete.zeros([multiAgent.output.lengthOutput, 1], tDomain); 
    end
    if ~isfield(exSig, "Sr")
       exSig.Sr = 0;
    end  
    if ~isfield(exSig, "Pr")
       exSig.Pr = eye(length(exSig.w));
    end  
    if NameValue.controller
        [PI_w, gain_w] = multiAgent.solveRegulatorEquationsReference(exSig.Sr, exSig.Pr);
        wAgg = kron(ones(multiAgent.N, 1), exSig.w.at(0)); 
        Kw = kron(eye(multiAgent.N), gain_w);
        controlInput = controlInput + Kw*wAgg;
    end
end

if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver")
   %Initialize simulation values
    for i = 1:2:length(NameValue.initialCondition)
       if NameValue.initialCondition{i} == "S"
           SObs = NameValue.initialCondition{i+1};
       end
       if NameValue.initialCondition{i} == "P"
           PObs = NameValue.initialCondition{i+1};
       end
       if NameValue.initialCondition{i} == "v"
           vObs = NameValue.initialCondition{i+1};
       end
    end
   xObs = multiAgent.network.combineInitialCondition(plantStateSpace, icObs); 
   xIcOut = plantStateSpace.C*x;
   xObsIcOut = observerStateSpace.C*xObs;
   measurement = xIcOut(contains(plantStateSpace.OutputName, "agent" + digitsPattern(1) + ".measurement"));
   measurementPrediction = xObsIcOut(contains(observerStateSpace.OutputName, "agent" + digitsPattern(1) + ".observer.measurement"));
   measurementError = measurement - measurementPrediction;
   obsIn = zeros(length(observerStateSpace.InputName), 1);
   obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".error.measurement")) = measurementError;
   obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".measurement")) = measurement;
   obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".disturbance")) = kron(ones(1,multiAgent.N),PObs).*kron(eye(multiAgent.N), ones(size(PObs,1)/multiAgent.N, size(PObs,2)))*vObs(:,end);
   if NameValue.controller
       controlInput = controlInput + uStabDisc * xObsIcOut(contains(plantStateSpace.OutputName, "x"));
       if exist("vObs", 'var') && exist("SObs", 'var')
           for indAgent = 1 : multiAgent.N
                SObsAgent = SObs((indAgent-1)*size(SObs,2)+1:indAgent*size(SObs,2), :);
                PObsAgent = PObs((indAgent-1)*(size(PObs,1)/multiAgent.N)+1:indAgent*(size(PObs,1)/multiAgent.N), :);
                if NameValue.controller
                    [PI_v.("agent"+indAgent), gain_v{indAgent}] = multiAgent.solveRegulatorEquationsDisturbance(SObsAgent, PObsAgent);
                end
           end
            Kv = blkdiag(gain_v{:});
            controlInput = controlInput + Kv*vObs;
       end
       obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".control")) = controlInput;
   end
   simDataObs = observerStateSpace.C*xObs + observerStateSpace.D*obsIn;
   yObs = simDataObs;
elseif isfield(exSig, "v") && NameValue.controller
    controlInput = controlInput + Kv*exSig.v.atIndex(1);
end

if NameValue.controller
    ud(1, contains(plantStateSpace.InputName, "control")) = controlInput;
end

simDataNy = plantStateSpace.C*x + plantStateSpace.D*ud(1,:).';
if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver") && multiAgent.agent1IsLeader
   SLead = exSig.S.on(t);
   SObs(1:size(SObs, 2),:) = SLead(1, :, :); 
   PLead = exSig.P.on(t);
   PObs(1:size(PObs, 1)/multiAgent.N,:) = PLead(1, :, :); 
end

for i = 1:length(t)-1
    if NameValue.controller && adaptationDelayController == 0
        solveRegEq = true;
        adaptationDelayController = NameValue.adaptationDelayController;
    elseif NameValue.controller
        adaptationDelayController = adaptationDelayController - 1;
        solveRegEq = false;
    end
    %Adapt observer parameters if one is specified
    if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver")
        if adaptationDelayObserver == 0
            L_v = {};
            for indAgent = 1 : multiAgent.N
                SObsAgent = SObs((indAgent-1)*size(SObs,2)+1:indAgent*size(SObs,2), :, i);
                PObsAgent = PObs((indAgent-1)*(size(PObs,1)/multiAgent.N)+1:indAgent*(size(PObs,1)/multiAgent.N), :, i);
                [L_v{end+1}, L] = NameValue.disturbanceObserver.getObserverGains(SObsAgent, PObsAgent);
                gainLocal = model.Input("agent" + indAgent + ".error.measurement", "B", L);
                B = tStep*multiAgent.network.agent(1).domain2indomain * gainLocal.getBDiscrete(multiAgent.network.agent(1).grid);
                if B == 0
                   B = zeros(size(multiAgent.network.agent(1).domain2indomain, 1), size(L, 2)); 
                end
                observerStateSpaceB(:, contains(observerStateSpace.InputName, "agent" + indAgent + ".error.measurement")) = kron([zeros(indAgent-1, 1);1; zeros(multiAgent.N-indAgent, 1)], B);
                if NameValue.controller
                    if solveRegEq
                        [PI_v.("agent"+indAgent), gain_v{indAgent}] = multiAgent.solveRegulatorEquationsDisturbance(SObsAgent, PObsAgent);
                    end
                end
            end
            L_vAgg = blkdiag(L_v{:});
            adaptationDelayObserver = NameValue.adaptationDelayObserver;
        else
            adaptationDelayObserver = adaptationDelayObserver - 1;
        end
        
        measurement = simDataNy(contains(plantStateSpace.OutputName, "agent" + digitsPattern(1) + ".measurement"), end);
        measurementPrediction = simDataObs(contains(observerStateSpace.OutputName, "agent" + digitsPattern(1) + ".observer.measurement"), end);
        measurementError = measurement - measurementPrediction;
        obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".measurement")) = measurement;
        %Get observer input
        obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".error.measurement")) = measurementError;
        %Compute updated state variables and output
        SObs = cat(3, SObs, SObs(:,:,i)+tStep*NameValue.disturbanceObserver.mu*kron(multiAgent.H_ext, eye(size(SObs, 2)))*(kron(ones(multiAgent.N, 1),squeeze(SLead(i,:,:)))-SObs(:,:,i)));
		if multiAgent.agent1IsLeader
			SObs(1:NameValue.disturbanceObserver.n_d,:,end) = SLead(i,:,:);
		end
		PObs = cat(3, PObs, PObs(:,:,i)+tStep*NameValue.disturbanceObserver.mu_P*kron(multiAgent.H_ext, eye(size(PObs, 1)/multiAgent.N))*(kron(ones(multiAgent.N, 1),reshape(PLead(i,:,:),[],size(PLead,3)))-PObs(:,:,i)));
		if multiAgent.agent1IsLeader
			PObs(1:NameValue.disturbanceObserver.n_r,:,end) = PLead(i,:,:);
		end
        vObs = cat(2, vObs, vObs(:,i)+tStep*(SObs(:,:,end)*kron(ones(1, multiAgent.N), eye(size(SObs, 2))).*kron(eye(multiAgent.N),...
            ones(size(SObs, 2)))*vObs(:,end) + L_vAgg*measurementError));
        if NameValue.controller
            obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".control")) = controlInput;
        end
        obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".disturbance")) = kron(ones(1,multiAgent.N),PObs(:,:,end)).*kron(eye(multiAgent.N), ones(size(PObs,1)/multiAgent.N, size(PObs,2)))*vObs(:,end);
        xObs = observerStateSpace.A*xObs + observerStateSpaceB*obsIn; 
        
        if NameValue.controller
            controlInput = 0;
            if isa(NameValue.referenceObserver, "ReferenceObserver")
               controlInput = controlInput + Kw*wAgg;
            end
            controlInput = controlInput + uStabDisc * yObs(contains(observerStateSpace.OutputName, "x"));
%             Kv = HkronIr * blkdiag(gain_v{:});
            Kv = blkdiag(gain_v{:});
            controlInput = controlInput + Kv*vObs(:, end);
            obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".control")) = controlInput;
        end
%         simDataObs = cat(2, simDataObs, observerStateSpace.C*xObs + observerStateSpace.D*obsIn);
        yObs = observerStateSpace.C*xObs + observerStateSpace.D*obsIn;
    end
    
    %Compute control input for this time step if NameValue.controller is
    %set to true
    x = plantStateSpace.A*x + plantStateSpace.B*ud(i,:).';
    if NameValue.controller
        controlInput = 0;
        if isa(NameValue.referenceObserver, "ReferenceObserver")
           wAgg = [];
           for indAgent = 1 - multiAgent.agent1IsLeader : multiAgent.N - multiAgent.agent1IsLeader 
               if solveRegEq
                   [PI_w.("agent"+indAgent), gain_w{indAgent + multiAgent.agent1IsLeader}] = multiAgent.solveRegulatorEquationsReference(simDataRef.("agent" + indAgent).S.atIndex(i), simDataRef.("agent" + indAgent).P.atIndex(i));
               end
               wAgg = cat(1, wAgg, simDataRef.("agent"+indAgent).w.atIndex(i));
           end
           Kw = blkdiag(gain_w{:});
           controlInput = controlInput + Kw*wAgg;
        else
            wAgg = kron(ones(multiAgent.N, 1), exSig.w.atIndex(i+1));
            controlInput = controlInput + Kw*wAgg;
        end
        if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver")
           controlInput = controlInput + uStabDisc * yObs(contains(observerStateSpace.OutputName, "x"));
           Kv = blkdiag(gain_v{:});
           controlInput = controlInput + Kv*vObs(:, end);
           obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".control")) = controlInput;
		else
		   u = ud(i+1);
		   u(contains(plantStateSpace.InputName, "control")) = controlInput + uStabDisc * simDataNy(contains(plantStateSpace.OutputName, "x"), end);
		   y = plantStateSpace.C*x + plantStateSpace.D*u.';
		   controlInput = controlInput + uStabDisc * y(contains(plantStateSpace.OutputName, "x"));
		   if isfield(exSig, "v")
			   controlInput = controlInput + Kv * exSig.v.atIndex(i + 1);
		   end
        end
        ud(i+1, contains(plantStateSpace.InputName, "control")) = controlInput;
	end
    %Update the state variable and compute the output using the discretized system
    simDataNy = cat(2, simDataNy, plantStateSpace.C*x + plantStateSpace.D*ud(i+1,:).');
    if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver")
        measurement = simDataNy(contains(plantStateSpace.OutputName, "agent" + digitsPattern(1) + ".measurement"), end);
        obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".measurement")) = measurement;
        simDataObs = cat(2, simDataObs, observerStateSpace.C*xObs + observerStateSpace.D*obsIn);
    end
    
end
%Format raw data into function output
simOutput = simDataNy.';
outputNames = plantStateSpace.OutputName;
grid = multiAgent.network.getSpatialGridsOfSimulation();
if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver")
    simOutput = [simOutput, simDataObs.'];
    outputNames = [outputNames; observerStateSpace.OutputName];
    grid{1} = cat(1, grid{1}, replace(grid{1}, "x", "observer.x"));
    grid{2} = cat(1, grid{2}, grid{2});
end
simData = stateSpace.simulationOutput2Quantity(...
	simOutput, t, outputNames, "z", grid);
simData.predefinedInput = u;
if isa(NameValue.disturbanceObserver, "AdaptiveDisturbanceObserver")
    SObs = quantity.Discrete(permute(SObs(:,:,1:end), [3 1 2]), quantity.Domain("t", t));
	PObs = quantity.Discrete(permute(PObs(:,:,1:end), [3 1 2]), quantity.Domain("t", t));
    vObs = quantity.Discrete(permute(vObs(:,1:end), [2 1]), quantity.Domain("t", t));
    for iAgent = 1:multiAgent.N
        simData.("agent"+iAgent).observer.S = SObs((iAgent-1)*size(SObs,2)+1:iAgent*size(SObs,2), :);
		simData.("agent"+iAgent).observer.P = PObs((iAgent-1)*size(PObs,1)/multiAgent.N+1:iAgent*size(PObs,1)/multiAgent.N, :);
        simData.("agent"+iAgent).observer.v = vObs((iAgent-1)*size(SObs,2)+1:iAgent*size(SObs,2));
    end
end
if isa(NameValue.referenceObserver, "ReferenceObserver")
    for iAgent = 1 - multiAgent.agent1IsLeader : multiAgent.N - multiAgent.agent1IsLeader 
        simData.("agent"+(iAgent + multiAgent.agent1IsLeader)).observer.Sr = simDataRef.("agent"+iAgent).S;
		simData.("agent"+(iAgent + multiAgent.agent1IsLeader)).observer.Pr = simDataRef.("agent"+iAgent).P;
        simData.("agent"+(iAgent + multiAgent.agent1IsLeader)).observer.w = simDataRef.("agent"+iAgent).w;
    end
end
end

