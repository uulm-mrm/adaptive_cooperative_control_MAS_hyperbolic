classdef MultiAgentHeterogenous
	%MULTIAGENTHETEROGENOUS Summary of this class goes here
	%   Detailed explanation goes here
	
	properties (SetAccess = protected)
		adjacencyMatrix double;
		network (1,1) model.Network;
		singleSys cell;
		informedAgents (:,1);
	end
	properties (Dependent = true)
        H double;   %Leader-follower matrix
        L double;   %Laplacian matrix
        N (1,1) {mustBeInteger, mustBeNonnegative};	% number of agents
		H_ext double;
		agent1IsLeader
    end
	methods
		function obj = MultiAgentHeterogenous(hyperbolicSystems, adjacencyMatrix, NameValue)
			%MULTIAGENTHETEROGENOUS Construct an instance of the class
			%MultiAgentHeterogenous which allows to simulate heterogenous
			%multi agent systems of hyperbolic agents
			arguments
				hyperbolicSystems {cell} = {};
				adjacencyMatrix double = 0;
				NameValue.nCoefRef (1,1) {mustBeInteger, mustBePositive} = 20;
                NameValue.nCoefDist (1,1) {mustBeInteger, mustBePositive} = 20;
                NameValue.informedAgents (:,1) = adjacencyMatrix(:,1);
			end
			obj.adjacencyMatrix = adjacencyMatrix;
			obj.informedAgents = NameValue.informedAgents;
			
			obj.singleSys = {obj.N};
			agents = [];
			for idx = 1:length(hyperbolicSystems)
				unitVector = [zeros(1, idx-1), 1, zeros(1, obj.N - idx)];
				agent = hyperbolicSystems{idx}.copy;

				obj.singleSys{idx} = MultiAgent(0, agent.Lambda, "A", agent.A, "G1", agent.input.input(2).B, "G2", agent.input.input(2).B0, ...
					"G3", agent.input.input(2).B1, "G4", agent.input.input(2).D.value, "diffusiveDisturbance", false, "Q0", agent.Q0, "Q1", agent.Q1, ...
					"output", agent.output.output, "nCoefRef", NameValue.nCoefRef, "nCoefDist", NameValue.nCoefDist);
				obj.singleSys{idx}.network.agent.output.add(agent.output.output);
				obj.singleSys{idx}.network.agent.output.remove(obj.singleSys{idx}.network.agent.output.output(1).type);
				distInput = model.Input("disturbance", "B", kron(unitVector, agent.input.input(2).B), "B0", kron(unitVector, agent.input.input(2).B0),...
					"B1", kron(unitVector, agent.input.input(2).B1), ...
					"D", misc.Gain("disturbance", kron(unitVector, agent.input.input(2).D.value), "outputType", "agent" + idx + ".controlOutput"));
				controlInput = agent.input.input(1).strrepType("control", "agent" + idx + ".control");
				measurementOutput = model.Output("agent" + idx + ".measurement", "C1", agent.e2.');
				agent = hyperbolicSystems{idx}.copyAndReplace("input", controlInput+distInput, "output", agent.output + measurementOutput);
				agents = cat(1, agents, agent);
			end
			obj.network = model.Network(agents, obj.adjacencyMatrix, obj.adjacencyMatrix(:, 1));
		end
		
		function simData = adaptiveControlSimulation(obj, simulationSetting, distModelDim, distOutDim, referenceObserver, NameValue)
			%ADAPTIVECONTROLSIMULATION Runs a simulation in a for-loop
			%   A controller and observer is designed for each agent
			arguments
				obj;
    			simulationSetting;
				distModelDim (1,1) double {mustBeInteger, mustBeNonnegative};	%Dimension of v(t)
				distOutDim (1,1) double {mustBeInteger, mustBeNonnegative};		%Dimension of d(t)
				referenceObserver (1,1) ReferenceObserver;						%The reference observer to be used
    			NameValue.exogenousInput = {};
    			NameValue.initialCondition = {};
				%Parameters of the adaptive disturbance observer
				NameValue.mu_S (1,1) double = 1;				%Observer gain for S
				NameValue.mu_P (1,1) double = 1;				%Observer gain for P
    			NameValue.adaptationDelayObserver = 0;			%Determines how frequently the decoupling eqautions are solved
				NameValue.distCoef (1,1) double {mustBeInteger, mustBeNonnegative} = 5;	%Determines the number of coefficients used to solve the decoupling equations
				NameValue.targetEvOde (1,:) double = -1*[1, 2];		%The target eigenvalues of the ODE of e_v(t)
				%Paramters of the controller
    			NameValue.adaptationDelayController = 0;		%Determines how frequently the regulator eqautions are solved
				NameValue.checkControllability = false;
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
			if isa(exSig.S, "double")
				exSig.S = quantity.Discrete.ones(1, tDomain)*exSig.S;
			end

			%Prepare an adaptive disturbance observer for each agent
			disturbanceObserver = {obj.N};
			for idx = 1:obj.N
				disturbanceObserver{idx} = AdaptiveDisturbanceObserver(obj.singleSys{idx}, distModelDim, ...
					distOutDim, NameValue.mu_S, "mu_P", NameValue.mu_P, "nCoef", NameValue.distCoef, "targetEvOde", NameValue.targetEvOde);
			end
			SObs = kron(ones(obj.N,1),zeros(distModelDim));
			PObs = kron(ones(obj.N,1),zeros(distOutDim, distModelDim));
			vObs = zeros(size(SObs, 1), 1);
			icObs = {};
			for idx = 1:2:length(NameValue.initialCondition)
   				if contains(NameValue.initialCondition{idx}, "observer.x")
      				icObs{end+1} = strrep(NameValue.initialCondition{idx}, "observer.x", "x");
      				icObs{end+1} = NameValue.initialCondition{idx+1};
   				end
			end
			for indAgent = 1 : obj.N
				%Add measurement output to MulitAgent objects
    			measurementOutput = model.Output("measurement", "C1", obj.network.agent(indAgent).e2.');
    			obj.singleSys{indAgent}.network.agent.output.add(measurementOutput);
    			stateSpaceCell{indAgent} = disturbanceObserver{indAgent}.getStateSpaceApproximation...
					(quantity.Discrete.zeros([obj.singleSys{indAgent}.n, obj.singleSys{indAgent}.m], obj.singleSys{indAgent}.domain),...
        			"t", t, "agent", indAgent);
			end
			observerStateSpace = stateSpace.parallel(stateSpaceCell{:});
			observerStateSpaceB = observerStateSpace.B;
			adaptationDelayObserver = 0;
			if isfield(exSig, "P") && isa(exSig.P, 'double')
				exSig.P = exSig.P*quantity.Discrete.ones(1, tDomain);
			end
			if ~isfield(exSig, "Sr")
				exSig.Sr = 0;
			end  
			
			%Get discrete state space representation
			plantStateSpace = obj.network.getSimulationModel(simulationSetting);
			plantStateSpace = stateSpace.addInput2Output(plantStateSpace);
			
			% stitch initial condition together
			x = obj.network.combineInitialCondition(plantStateSpace, NameValue.initialCondition);
			
			% stitch exogenous signals together
			u = stateSpace.combineInputSignals(plantStateSpace, t, NameValue.exogenousInput{:});
			ud = u.on(t);
			
			%Compute control input parameters
			for idx = 1:obj.N
    			obj.singleSys{idx}.setBacksteppingKernel();
    			obj.singleSys{idx}.setRegulatorEquationRefCoefficients(size(exSig.Sr,1));
				obj.singleSys{idx}.setRegulatorEquationDistCoefficients(size(exSig.S,1));
    			uStab{idx} = obj.singleSys{idx}.getStabilizingStateFeedback();
    			uStab{idx}.setApproximation(obj.network.agent(idx).grid);
				uStabCop{idx} = uStab{idx}.Cop;
			end
			uStabDisc = blkdiag(uStabCop{:});
			controlInput = 0;
			adaptationDelayController = 0;
			
			%Run simulation of reference observer
			if ~isfield(exSig, "w")
				exSig.w = quantity.Discrete.zeros([referenceObserver.nDimensions, 1], tDomain); 
			end
			if ~isfield(exSig, "Sr")
				exSig.Sr = zeros(referenceObserver.nDimensions);
			end
			if ~isfield(exSig, "Pr")
				exSig.Pr = eye(referenceObserver.nDimensions);
			end
			simDataRef = referenceObserver.simulate(exSig.w, exSig.Sr, exSig.Pr); 
   			wAgg = [];
   			for indAgent = 0 : obj.N-1
       			[PI_w.("agent"+indAgent), gain_w{indAgent + 1}] = obj.singleSys{indAgent+1}.solveRegulatorEquationsReference(simDataRef.("agent" + indAgent).S.atIndex(1), simDataRef.("agent" + indAgent).P.atIndex(1));
       			wAgg = cat(1, wAgg, simDataRef.("agent"+indAgent).w.atIndex(1));
   			end
   			Kw = blkdiag(gain_w{:});
   			controlInput = controlInput + Kw*wAgg;
			
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
   			xObs = obj.network.combineInitialCondition(plantStateSpace, icObs); 
   			xIcOut = plantStateSpace.C*x;
   			xObsIcOut = observerStateSpace.C*xObs;
   			measurement = xIcOut(contains(plantStateSpace.OutputName, "agent" + digitsPattern(1) + ".measurement"));
   			measurementPrediction = xObsIcOut(contains(observerStateSpace.OutputName, "agent" + digitsPattern(1) + ".observer.measurement"));
   			measurementError = measurement - measurementPrediction;
   			obsIn = zeros(length(observerStateSpace.InputName), 1);
   			obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".error.measurement")) = measurementError;
   			obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".measurement")) = measurement;
   			obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".disturbance")) = kron(ones(1,obj.N),PObs).*kron(eye(obj.N), ones(size(PObs,1)/obj.N, size(PObs,2)))*vObs(:,end);
   			controlInput = controlInput + uStabDisc * xObsIcOut(contains(plantStateSpace.OutputName, "x"));
   			if exist("vObs", 'var') && exist("SObs", 'var')
       			for indAgent = 1 : obj.N
            			SObsAgent = SObs((indAgent-1)*size(SObs,2)+1:indAgent*size(SObs,2), :);
            			PObsAgent = PObs((indAgent-1)*(size(PObs,1)/obj.N)+1:indAgent*(size(PObs,1)/obj.N), :);

            			[PI_v.("agent"+indAgent), gain_v{indAgent}] = obj.singleSys{indAgent}.solveRegulatorEquationsDisturbance(SObsAgent, PObsAgent);
       			end
        			Kv = blkdiag(gain_v{:});
        			controlInput = controlInput + Kv*vObs;
   			end
   			obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".control")) = controlInput;
   			simDataObs = observerStateSpace.C*xObs + observerStateSpace.D*obsIn;
   			yObs = simDataObs;
			
			ud(1, contains(plantStateSpace.InputName, "control")) = controlInput;
			simDataNy = plantStateSpace.C*x + plantStateSpace.D*ud(1,:).';
   			SLead = exSig.S.on(t);
   			SObs(1:size(SObs, 2),:) = SLead(1, :, :);
   			PLead = exSig.P.on(t);
   			PObs(1:size(PObs, 1)/obj.N,:) = PLead(1, :, :);
        	L_v = {};
			L_x = {};
			averageIterationTime = 0;
			runningTime = 0;
			h = waitbar(0, "Simulation in progress");
			for i = 1:length(t)-1
				waitbar(i/(length(t)-1), h, "Simulation in progress, ETA: " + ceil(((length(t)-1-i)*averageIterationTime)) + "s");
				startingTime = tic;
    			if adaptationDelayController == 0
        			solveRegEq = true;
        			adaptationDelayController = NameValue.adaptationDelayController;
				else
        			adaptationDelayController = adaptationDelayController - 1;
        			solveRegEq = false;
    			end
    			%Adapt observer parameters if one is specified
    			if adaptationDelayObserver == 0
					endidxprev = 0;
        			for indAgent = 1 : obj.N
            			SObsAgent = SObs((indAgent-1)*size(SObs,2)+1:indAgent*size(SObs,2), :, i);
            			PObsAgent = PObs((indAgent-1)*(size(PObs,1)/obj.N)+1:indAgent*(size(PObs,1)/obj.N), :, i);
						if NameValue.checkControllability
							%Check controllability condition
							outputTil = obj.singleSys{indAgent}.output.backstepping(obj.singleSys{indAgent}.backsteppingKernel, "inverse", true);
							invLambda = 1/obj.singleSys{indAgent}.Lambda;
							phi = int(subs(invLambda, "z", "zeta"), "zeta", 0, "z");
							eigVals = eig(SObsAgent);
							controllable = true;
							for indEv = 1:length(eigVals)
						    	psi = expm(eigVals(indEv)*(phi - phi.subs("z", "zeta")));
						    	M = psi.subs("zeta", 0)*(obj.singleSys{indAgent}.E1+obj.singleSys{indAgent}.E2*obj.singleSys{indAgent}.Q0)...
						        	- int(psi*invLambda.subs("z", "zeta")*obj.singleSys{indAgent}.backsteppingKernel.getA0Target().subs("z","zeta"), "zeta", 0, "z");
						    	N = outputTil.out(M);
						    	if det(N) == 0
									controllable = false; 
						    	end
							end
							if controllable
								[L_v{indAgent}, L_x{indAgent}] = disturbanceObserver{indAgent}.getObserverGains(SObsAgent, PObsAgent);
							else
								warning("Controllability issue on time step " + i);
							end
						else
							[L_v{indAgent}, L_x{indAgent}] = disturbanceObserver{indAgent}.getObserverGains(SObsAgent, PObsAgent);
						end

            			gainLocal = model.Input("agent" + indAgent + ".error.measurement", "B", L_x{indAgent});
            			B = tStep*obj.network.agent(indAgent).domain2indomain * gainLocal.getBDiscrete(obj.network.agent(indAgent).grid);
            			if B == 0
							B = zeros(size(obj.network.agent(indAgent).domain2indomain, 1), size(L_x{indAgent}, 2)); 
            			end
						observerStateSpaceB(:, contains(observerStateSpace.InputName, "agent" + indAgent + ".error.measurement")) = [zeros(endidxprev, size(B, 2)); B; zeros(size(observerStateSpaceB, 1)-length(B) - endidxprev, size(B, 2))];
						endidxprev = endidxprev + length(B);
            			if solveRegEq
                			[PI_v.("agent"+indAgent), gain_v{indAgent}] = obj.singleSys{indAgent}.solveRegulatorEquationsDisturbance(SObsAgent, PObsAgent);
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
    			SObs = cat(3, SObs, SObs(:,:,i)+tStep*NameValue.mu_S*kron(obj.H_ext, eye(size(SObs, 2)))*(kron(ones(obj.N, 1),squeeze(SLead(i,:,:)))-SObs(:,:,i)));
				SObs(1:distModelDim,:,end) = SLead(i,:,:);
				
				PObs = cat(3, PObs, PObs(:,:,i)+tStep*NameValue.mu_P*kron(obj.H_ext, eye(size(PObs, 1)/obj.N))*(kron(ones(obj.N, 1),reshape(PLead(i,:,:),[],size(PLead,3)))-PObs(:,:,i)));
				PObs(1:distOutDim,:,end) = PLead(i,:,:);
    			vObs = cat(2, vObs, vObs(:,i)+tStep*(SObs(:,:,end)*kron(ones(1, obj.N), eye(size(SObs, 2))).*kron(eye(obj.N),...
        			ones(size(SObs, 2)))*vObs(:,end) + L_vAgg*measurementError));
    			obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".control")) = controlInput;
    			obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".disturbance")) = kron(ones(1,obj.N),PObs(:,:,end)).*kron(eye(obj.N), ones(size(PObs,1)/obj.N, size(PObs,2)))*vObs(:,end);
    			xObs = observerStateSpace.A*xObs + observerStateSpaceB*obsIn; 
    			
    			controlInput = 0;
       			controlInput = controlInput + Kw*wAgg;
    			controlInput = controlInput + uStabDisc * yObs(contains(observerStateSpace.OutputName, "x"));
	%             Kv = HkronIr * blkdiag(gain_v{:});
    			Kv = blkdiag(gain_v{:});
    			controlInput = controlInput + Kv*vObs(:, end);
    			obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".control")) = controlInput;
		%         simDataObs = cat(2, simDataObs, observerStateSpace.C*xObs + observerStateSpace.D*obsIn);
    			yObs = observerStateSpace.C*xObs + observerStateSpace.D*obsIn;
			
				%Compute control input for this time step if NameValue.controller is
				%set to true
				controlInput = 0;
   				wAgg = [];
   				for indAgent = 0 : obj.N - 1 
       				if solveRegEq
						if NameValue.checkControllability
							%Check controllability condition
							outputTil = obj.singleSys{indAgent+1}.output.backstepping(obj.singleSys{indAgent+1}.backsteppingKernel, "inverse", true);
							invLambda = 1/obj.singleSys{indAgent+1}.Lambda;
							phi = int(subs(invLambda, "z", "zeta"), "zeta", 0, "z");
							eigVals = eig(simDataRef.("agent" + indAgent).S.atIndex(i));
							controllable = true;
							for indEv = 1:length(eigVals)
						    	psi = expm(eigVals(indEv)*(phi - phi.subs("z", "zeta")));
						    	M = psi.subs("zeta", 0)*(obj.singleSys{indAgent+1}.E1+obj.singleSys{indAgent+1}.E2*obj.singleSys{indAgent+1}.Q0)...
						        	- int(psi*invLambda.subs("z", "zeta")*obj.singleSys{indAgent+1}.backsteppingKernel.getA0Target().subs("z","zeta"), "zeta", 0, "z");
						    	N = outputTil.out(M);
						    	if det(N) == 0
									controllable = false; 
						    	end
							end
							if controllable
								[PI_w.("agent"+indAgent), gain_w{indAgent + 1}] = ...
									obj.singleSys{indAgent + 1}.solveRegulatorEquationsReference(simDataRef.("agent" + indAgent).S.atIndex(i), simDataRef.("agent" + indAgent).P.atIndex(i));
							else
								warning("Controllability issue on time step " + i);
							end
						else
							[PI_w.("agent"+indAgent), gain_w{indAgent + 1}] = ...
									obj.singleSys{indAgent + 1}.solveRegulatorEquationsReference(simDataRef.("agent" + indAgent).S.atIndex(i), simDataRef.("agent" + indAgent).P.atIndex(i));
						end
					end
       				wAgg = cat(1, wAgg, simDataRef.("agent"+indAgent).w.atIndex(i));
   				end
   				Kw = blkdiag(gain_w{:});
   				controlInput = controlInput + Kw*wAgg;
   				controlInput = controlInput + uStabDisc * yObs(contains(observerStateSpace.OutputName, "x"));
   				Kv = blkdiag(gain_v{:});
   				controlInput = controlInput + Kv*vObs(:, end);
   				obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".control")) = controlInput;
				
				ud(i+1, contains(plantStateSpace.InputName, "control")) = controlInput;
				%Update the state variable and compute the output using the discretized system
				x = plantStateSpace.A*x + plantStateSpace.B*ud(i,:).';
				simDataNy = cat(2, simDataNy, plantStateSpace.C*x + plantStateSpace.D*ud(i+1,:).');
				measurement = simDataNy(contains(plantStateSpace.OutputName, "agent" + digitsPattern(1) + ".measurement"), end);
				obsIn(contains(observerStateSpace.InputName, "agent" + digitsPattern(1) + ".measurement")) = measurement;
				simDataObs = cat(2, simDataObs, observerStateSpace.C*xObs + observerStateSpace.D*obsIn);
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
			grid = obj.network.getSpatialGridsOfSimulation();
			simOutput = [simOutput, simDataObs.'];
			outputNames = [outputNames; observerStateSpace.OutputName];
			grid{1} = cat(1, grid{1}, replace(grid{1}, "x", "observer.x"));
			grid{2} = cat(1, grid{2}, grid{2});
			simData = stateSpace.simulationOutput2Quantity(...
				simOutput, t, outputNames, "z", grid);
			simData.predefinedInput = u;
			SObs = quantity.Discrete(permute(SObs(:,:,1:end), [3 1 2]), quantity.Domain("t", t));
			PObs = quantity.Discrete(permute(PObs(:,:,1:end), [3 1 2]), quantity.Domain("t", t));
			vObs = quantity.Discrete(permute(vObs(:,1:end), [2 1]), quantity.Domain("t", t));
			for iAgent = 1:obj.N
    			simData.("agent"+iAgent).observer.S = SObs((iAgent-1)*size(SObs,2)+1:iAgent*size(SObs,2), :);
				simData.("agent"+iAgent).observer.P = PObs((iAgent-1)*size(PObs,1)/obj.N+1:iAgent*size(PObs,1)/obj.N, :);
    			simData.("agent"+iAgent).observer.v = vObs((iAgent-1)*size(SObs,2)+1:iAgent*size(SObs,2));
			end
			for iAgent = 0 : obj.N - 1 
    			simData.("agent"+(iAgent + 1)).observer.Sr = simDataRef.("agent"+iAgent).S;
				simData.("agent"+(iAgent + 1)).observer.Pr = simDataRef.("agent"+iAgent).P;
    			simData.("agent"+(iAgent + 1)).observer.w = simDataRef.("agent"+iAgent).w;
			end
		end

		%% get dependent properties
        function L = get.L(obj)
           L =  diag(sum(obj.adjacencyMatrix, 2)) - obj.adjacencyMatrix;
        end % get.L()

        function agent1IsLeader = get.agent1IsLeader(obj)
            agent1IsLeader = isempty(find(obj.adjacencyMatrix(1,:)));
        end % agent1IsLeader
        
        function H = get.H(obj)
            if obj.agent1IsLeader
               H = obj.H_ext(2:end, 2:end);
            else
               H = obj.H_ext; 
            end
        end
        
		function H_ext = get.H_ext(obj)
			H_ext = obj.L + diag(obj.informedAgents);
		end % get.H()

        function N = get.N(obj)
           N = size(obj.adjacencyMatrix, 1);
        end % get.N()
	end
end

