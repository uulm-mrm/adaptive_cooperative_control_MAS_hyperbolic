classdef MultiAgent < handle & matlab.mixin.Copyable
    % MultiAgent represents a System composed by N agents that communicate
    % according to a communication graph G.
	% The agents, denoted by i, are defined by eqations of the form
    %   
	%   xi_t(z,t) = Lambda(z) xi_z(z,t) + A(z) xi(z,t) + G1(z)(sum_j=1^N a_ij(s_ij d_i - d_j) + a_i0 d_i)
    %   xi^+(0) = Q0 xi^-(0,t) + G2 (sum_j=1^N a_ij(s_ij d_i - d_j) + a_i0 d_i)
    %   xi^-(1) = Q1 xi^+(1,t) + u_i + G3 (sum_j=1^N a_ij(s_ij d_i - d_j) + a_i0 d_i)
    %   y_i = output[xi(z,t)] + G4  (sum_j=1^N a_ij(s_ij d_i - d_j) + a_i0 d_i)
    %   ni = xi^-(0)
    %
    %       Where s_ij = sgn(a_ij)
    %
	%
	% and the initial condition xi(z,0) = xi0(z).
	% The distributed state is x(z, t) \in R^n, with the spatial variable z in [0, 1] and the 
	% temporal variable t. The state vector can be splitted into x = col(x_1, x_2), with
	%
	%	xi^-(z, t) = E1.' x(z, t) \in R^p describing the states related to transport from z=1 to z=0,
	%	xi^+(z, t) = E2.' x(z, t) \in R^m describing the states related to transport from z=0 to z=1.
	%
    %The aggregation x = col(x1, ..., xN) leads to the aggregated System
    %
    %   x_t(z,t) = (I kron Lambda(z)) x_z(z,t) + (I kron A(z)) x(z,t) + (H kron G1(z)) d
    %   x^+(0) = (I_N kron Q0) x^-(0,t) + (H kron G2) d
    %   x^-(1) = (I_N kron Q1) x^+(1,t) + u + (H kron G3) d
    %   y = (I kron output)[x(z,t)] + (H kron G4) d
    %   n = x^-(0)
	%
	% See also model.Transport.Transport (constructor), model.Dps, model.TransportOde, model.Wave,
	%	model.WaveOde.
    
    properties (SetAccess = protected)
        Lambda quantity.Discrete {mustBe.quadraticMatrix};
        A quantity.Discrete {mustBe.quadraticMatrix};
        Q0 double;
        Q1 double;

        adjacencyMatrix double; %adjacency matrix A of the communication Graph G
        informedAgents (:,1);
        
        E1 double;
        E2 double;
        
        G1 quantity.Discrete;
        G2 double;
        G3 double;
        G4 double;

        output (1,1) model.Output;
        
        name;			% name of object
        prefix string;	% is added to the type of output, the state name and simulation signals
        
        network;
        aggregatedSystem (1,1) model.Transport;
        
        backsteppingKernel;
        G1Backstepped;
        A0Backstepped;
        outputBackstepped;
        
        nCoefRef;       %Number of coefficients used for the numeric calculation of the regulator equations reference
        coefRegEq1 cell;
        coefRegEq2 cell;
        polynomialSFactor double;
        
        nCoefDist;      %Number of coefficients used for the numeric calculation of the regulator equations disturbance
        coefRegEqD1 cell;
        coefRegEqD2 cell;
        coefRegEqD3 cell;
        coefRegEqD4 cell;
		
		diffusiveDisturbance (1,1) logical;
    end
    
    properties (Dependent = true)
        H double;   %Leader-follower matrix
        H_ext double;   %Leader-follower matrix extended by zero rows and columns for the leader in case agent 1 is the leader
        agent1IsLeader (1, 1) logical;  %Is true if agent1 is treated as the leader
        L double;   %Laplacian matrix
                
        n (1,1) {mustBeInteger, mustBeNonnegative};	% number of PDEs, i.e. size(obj.Lambda, 1)
		p (1,1) {mustBeInteger, mustBeNonnegative};	% number of PDEs with positive lambda_i
		m (1,1) {mustBeInteger, mustBeNonnegative};	% number of PDEs with negative lambda_i
        N (1,1) {mustBeInteger, mustBeNonnegative};	% number of agents
        
        domain (1,1) quantity.Domain;

    end
    
    methods
        function obj = MultiAgent(adjacencyMatrix, Lambda, NameValue)
            %MultiAgent constructs an object of type MultiAgent, which models a cooperative 
            %system composed by N agents that communicate according to a
            %fixed graph G with the corresponding adjacency matrix A
            %
            %   obj = MultiAgent(adjacencyMatrix, Lambda) obj is a system of hyperbolic agents
            %       with transport velocities defined by Lambda. Usually, Lambda is diagonal, and positive diagonal-elements 
			%		describe from z=1 to z=0, while negative diagonal elements describe transport
			%		from z=0 to z=1. The communication between the agents
			%		is indicated by adjacencyMatrix
            %
            %   obj = MultiAgent([...], "A", A) sets the distributed matrix A. A must be a
			%		quantity.Discrete with size and domain fitting to Lambda.
            %
            %   obj = MultiAgent([...], "Q0", Q0) sets the boundary operator for z=0
            %
            %   obj = MultiAgent([...], "Q1", Q1) sets the boundary operator for z=1
            %
            %   obj = MultiAgent([...], "Gx", Gx) sets input gains 
            %       of the disturbance d
            %
            %   obj = MultiAgent([...], "output", output) sets output operator
            %
			%	model.Transport([...], "name", name) sets obj.name.
			%
			%	model.Transport([...], "prefix", prefix) sets obj.prefix, which is a string, which
			%		is added to the type of output, the state name and simulation signals to be able
			%		to differ between different systems.
            arguments
				adjacencyMatrix {mustBe.quadraticMatrix(adjacencyMatrix)} = [];
                Lambda {mustBe.quadraticMatrix(Lambda)} = [];
                NameValue.A {mustBe.quadraticMatrix(NameValue.A)} = quantity.Discrete.zeros(size(Lambda), Lambda.getDomain());
                NameValue.Q0 double = zeros(length(find(Lambda.at(0)<0)), length(find(Lambda.at(0)>0)));
                NameValue.Q1 double = zeros(length(find(Lambda.at(0)>0)), length(find(Lambda.at(0)<0)));
                NameValue.G1 = quantity.Discrete.zeros(size(Lambda), Lambda(1,1).domain);
                NameValue.G2 double = zeros(size(find(diag(Lambda.at(0))<=0), 1), size(Lambda, 2));
                NameValue.G3 double = zeros(size(find(diag(Lambda.at(0))>0), 1), size(Lambda, 2));
                NameValue.G4 double = zeros(size(find(diag(Lambda.at(0))>0), 1), size(Lambda, 2));
                NameValue.output (1,1) model.Output = model.Output('measurement', "C", quantity.Discrete.zeros(length(Lambda), Lambda.getDomain()));
				NameValue.name (1,1) string = "global.reference";
				NameValue.prefix (1,1) string = "multiAgent";
                NameValue.informedAgents (:,1) = adjacencyMatrix(:,1);
                NameValue.diffusiveDisturbance (1,1) logical = false;
                NameValue.nCoefRef (1,1) {mustBeInteger, mustBePositive} = 10;
                NameValue.nCoefDist (1,1) {mustBeInteger, mustBePositive} = 10;
			end % arguments
            assert(nargin >= 2, "adjacencyMatrix and Lambda are obligatory arguments");
            obj.adjacencyMatrix = adjacencyMatrix;
            obj.Lambda = Lambda;
            obj.A = NameValue.A;
            obj.Q0 = NameValue.Q0;
            obj.Q1 = NameValue.Q1;
            obj.G1 = NameValue.G1;
            obj.G2 = NameValue.G2;
            obj.G3 = NameValue.G3;
            obj.G4 = NameValue.G4;
            obj.output = NameValue.output;
            obj.name = NameValue.name;
            obj.prefix = NameValue.prefix;
            obj.informedAgents = NameValue.informedAgents;
            obj.nCoefRef = NameValue.nCoefRef;
            obj.nCoefDist = NameValue.nCoefDist;
			obj.diffusiveDisturbance = NameValue.diffusiveDisturbance;
            assert(det(obj.H) ~= 0 || obj.adjacencyMatrix == 0, "Graph must be connected");
            % Matrices to select states with positive or negative lambda_i:
            I = eye(size(Lambda));
            if isdiag(Lambda.at(0))
                positiveLambda = diag(Lambda.at(0)) > 0;
                obj.E1 = I(:, positiveLambda);
                obj.E2 = I(:, ~positiveLambda);
            elseif (obj.p == obj.m) && misc.iseye(Lambda(1:obj.p, (obj.p+1):end).at(0))
                % energy coordinates
                obj.E1 = eye(obj.n, obj.p);
                obj.E2 = [zeros(obj.m, obj.m); eye(obj.m, obj.m)];
            else
                warning("unrecognized structure of Lambda");
                obj.E1 = [];
                obj.E2 = [];
            end
            
            %Define the model.Network object
			if obj.N > 1
            	distCoupling = obj.H_ext^NameValue.diffusiveDisturbance;
            	aggG1 = kron(distCoupling, obj.G1);
            	aggG2 = kron(distCoupling, obj.G2);
            	aggG3 = kron(distCoupling, obj.G3);
            	aggG4 = kron(distCoupling, obj.G4);
            	agents = [];
            	outputType = obj.output.type;
            	for idx = 1:obj.N
                	controlInput = model.Input("agent" + idx + ".control", "B1", eye(obj.p));
                	disturbanceInput = model.Input("disturbance", "B", aggG1((idx - 1)*obj.n + 1:idx*obj.n, :),...
                    	"B0", aggG2((idx - 1)*size(obj.G2, 1) + 1:idx*size(obj.G2, 1), :), "B1", aggG3((idx - 1)*size(obj.G3, 1) + 1:idx*size(obj.G3, 1), :),...
                    	"D", misc.Gain("disturbance", aggG4((idx - 1)*size(obj.G4, 1) + 1:idx*size(obj.G4, 1), :), "outputType", "agent" + idx + "." + outputType));
                	input = controlInput + disturbanceInput;                
                	output = obj.output.copy();
                	output.strrepType(obj.output.type, "agent" + idx + "." + outputType);
                	agent = model.Transport(obj.Lambda, "A", obj.A, "input", input, "output", output, "Q0", obj.Q0, "Q1", obj.Q1, "prefix", "agent" + idx);
                	agents = cat(1, agents, agent);
            	end
            	obj.output.strrepType(obj.output.type, outputType);
            	obj.network = model.Network(agents, adjacencyMatrix, adjacencyMatrix(:, 1));
            	
            	%Define the aggregated System as a model.Transport object
            	I = eye(obj.N);
            	controlInput = model.Input('control', "B1", eye(obj.N*obj.p));
            	aggregatedDisturbance = model.Input('disturbance', "B", aggG1, "B0", aggG2, ...
                	"B1", aggG3, "D", misc.Gain('disturbance', aggG4));
            	obj.aggregatedSystem = model.Transport(kron(I, obj.Lambda), "A", kron(I, obj.A), "input", controlInput+aggregatedDisturbance, ...
                	"output", kron(I, obj.output), "Q0", kron(I, obj.Q0), "Q1", kron(I, obj.Q1));
			else
            	outputType = obj.output.type;
            	controlInput = model.Input("agent" + 1 + ".control", "B1", eye(obj.p));
            	disturbanceInput = model.Input("disturbance", "B", obj.G1,...
                	"B0", obj.G2, "B1", obj.G3,...
                	"D", misc.Gain("disturbance", obj.G4, "outputType", "agent" + 1 + "." + outputType));
            	input = controlInput + disturbanceInput;                
            	output = obj.output.copy();
            	output.strrepType(obj.output.type, "agent" + 1 + "." + outputType);
            	agent = model.Transport(obj.Lambda, "A", obj.A, "input", input, "output", output, "Q0", obj.Q0, "Q1", obj.Q1, "prefix", "agent" + 1);
            	obj.output.strrepType(obj.output.type, outputType);
				obj.network.agent = agent;
			end
            
        end
        
        function [simData, myStateSpace] = simulate(obj, ...
				simulationSetting, initialCondition, controller, exogenousSignal, observer, normalDisturbance)
            %calls model.Network.simulate
            arguments
               obj;
               simulationSetting cell = {'t', linspace(0, 1, obj.domain.n)};
               initialCondition cell = obj.getNeutralInitialConditions();
               controller = struct();
               exogenousSignal cell = cell(0);
               observer = struct([]);
               normalDisturbance (1,1) logical = false;
            end % arguments
			assert(obj.adjacencyMatrix ~= 0, "This function is not available for single-agent systems")
            for ind = 1:length(simulationSetting)
               if  (isstring(simulationSetting{ind})||ischar(simulationSetting{ind})) && contains(simulationSetting{ind}, "t")
                  tDomain = quantity.Domain("t", simulationSetting{ind+1}); 
               end
            end
            exogenousSignalTransformed = exogenousSignal;
            %Transform the reference signals using the solution of the regulator equations for reference tracking
            if isfield(controller, "PI_w")
                for ind = 1:length(exogenousSignalTransformed)
                   if isstring(exogenousSignalTransformed{ind}) && contains(exogenousSignalTransformed{ind}, "reference")
                      splitted = split(exogenousSignalTransformed{ind}, ".");
                      origRef.(splitted(1)) = exogenousSignalTransformed{ind + 1};
                      exogenousSignalTransformed{ind + 1} = obj.E1.' * subs(controller.PI_w.(splitted(1)), "z", 1) * exogenousSignalTransformed{ind + 1};
                   end
                end
            end
            
            %If disturbanceState is not specified for each controller
            for ind = 1:length(exogenousSignalTransformed)
               if strcmp(exogenousSignalTransformed{ind}, "disturbanceState")
                   dim = length(exogenousSignalTransformed{ind + 1})/obj.N;
                   for indAgent = 1:obj.N
                       exogenousSignalTransformed = [exogenousSignalTransformed, ...
                           {"agent" + indAgent + ".disturbanceState", exogenousSignalTransformed{ind + 1}((indAgent - 1)*dim + 1:indAgent*dim)}];
                       disturbanceStateAggregated = kron(obj.H_ext, eye(size(obj.G1))) * exogenousSignalTransformed{ind + 1};
                       exogenousSignalTransformed = [exogenousSignalTransformed, ...
                           {"agent" + indAgent + ".network.disturbanceState", disturbanceStateAggregated((indAgent - 1)*dim + 1:indAgent*dim)}];
                   end
               elseif strcmp(exogenousSignalTransformed{ind}, "disturbance") && normalDisturbance
                   %Extract individual disturbances from aggregated
                   %disturbance signal
                   for indAgent = 1:obj.N
                       exogenousSignalTransformed = [exogenousSignalTransformed, ...
                           {"agent" + indAgent + ".disturbance", exogenousSignalTransformed{ind+1}((indAgent - 1)*size(obj.G2, 2) + 1:indAgent*size(obj.G2, 2))}];
                   end
               end
            end
            
            if isfield(controller, "input")
                controlInput = controller.input;
            else
                controlInput = {};
            end
            if ~isempty(observer)
                if ~isempty(controlInput)
                    for indAgent = 1:obj.N
                        observerFeedback = controlInput.("agent" + indAgent){2};
                        observerFeedback = observerFeedback.copyAndReplace("input", ...
                            observerFeedback.input.strrepInputType("agent" + indAgent + ".network.disturbanceState", "agent" + indAgent + ".network.observer.disturbanceState"));
                        controlInput.("agent" + indAgent){2} = observerFeedback;
                        controlInput.("agent" + indAgent) = [controlInput.("agent" + indAgent), observer.observer.("agent" + indAgent)];
                    end
                else
                   controlInput = observer.observer; 
                end
                network = obj.network.copyAndReplace("sharedSignals", ["observer.disturbanceState", "error.measurement"], "agent", copy(obj.network.agent));
                simulationSetting = [simulationSetting, {"errorSignalPrefix", true}];
                measurementOutput = model.Output("measurement", "C0", obj.E1.');
                network.addCommunication2agents(measurementOutput);
            else
                network = obj.network.copy();
            end
            if normalDisturbance
               for indAgent = 1:obj.N
                  distInput = model.Input("agent" + indAgent + ".disturbance", "B", obj.G1, "B0", obj.G2, "B1", obj.G3, ...
                      "D", misc.Gain("agent" + indAgent + ".disturbance",...
                      obj.G4, "outputType", "agent" + indAgent + "." + obj.output.type));
                  network.agent(indAgent).input.remove("disturbance");
                  network.agent(indAgent).input.add(distInput);
               end
            end
            [simData, myStateSpace] = network.simulate(simulationSetting, initialCondition, controlInput, exogenousSignalTransformed);
            
            %Change the reference signal in the simData if it was altered
            %before
            if exist("origRef", 'var') == 1
                for indAgent = 1:obj.N
                   simData.("agent" + indAgent).reference = origRef.("agent" + indAgent); 
                end
            end
            
        end
        
        function setBacksteppingKernel(obj)
            %Computes the backstepping kernel and saves it as well as G1,
            %A0 and the output operator of the target system
            obj.backsteppingKernel = kernel.Feedback.tryToLoadKernel("kernel.Feedback", obj.Lambda, "A", obj.A, "Q0", obj.Q0);
            obj.G1Backstepped = obj.backsteppingKernel.backsteppingTransformation(obj.G1)+ obj.backsteppingKernel.getValue.subs("zeta", 0) * obj.Lambda.at(0)...
							* obj.E2 * obj.G2;
            obj.A0Backstepped = obj.backsteppingKernel.getA0Target();
            obj.outputBackstepped = obj.output.backstepping(obj.backsteppingKernel, "inverse", true);
        end
        
        function [systemTarget, feedforwardControl] = feedforwardControl(obj, S_r, P_r, S, P)
            %feedforwardControl computes a feedback using Backstepping
            %transformation as well as a controller that ensures
            %disturbance rejection and reference tracking for the
            %Parameters S_r, P_r, S and P. Using the control input causes every
            %agent to show the behaviour if systemTarget.
            %
            %   systemTarget = feedforwardControl(obj, S_r, P_r, S, P)
            %       returns the target System as a model.Transport object.
            %       S_r can be either a numberic matrix or a struct
            %       containing individual, time variant S_r for each agent
            %       as a quantity.Discrete object like 
            %           S_r.agent1 = quantity.Discrete(...)
            %
            %   [systemTarget, feedforwardControl] = feedforwardControl(obj, S_r, P_r, S, P)
            %       returns the target System and a struct containing the
            %       feedback input and several interim results. The struct
            %       feedforwardControl can be given to the simulate
            %       function directly.
            
            %Compute Backstepping Kernel
            obj.setBacksteppingKernel();
            feedforwardControl.backsteppingKernel = obj.backsteppingKernel;
            
            %Save variables that appear in the regulator equations
            A0tilde = obj.A0Backstepped;      
            
            %Solve regulator equations (disturbance)
            obj.setRegulatorEquationDistCoefficients(size(S,1));
            [PI_v, disturbanceGain] = obj.solveRegulatorEquationsDisturbance(S, P);
            feedforwardControl.PI_v = PI_v;
            
            %Solve regulator equations (reference)
            if isstruct(S_r)              
                %Compute uDash for each agent
                uDash = {};
                for indAgent = 1:obj.N
                   [PI_w, PI_wPolynomCoefficients] = obj.solveRegulatorEquationReferenceExplicit(S_r.("agent"+indAgent), P_r);
                   feedforwardControl.disturbanceGain.("agent" + indAgent) = disturbanceGain;
                   uDash{indAgent, 1} = misc.Gain("agent" + indAgent + ".reference", eye(size(obj.E1, 2)), "outputType", "agent" + indAgent + "." + 'control') + ...
                       misc.Gain("agent" + indAgent + "." + obj.network.prefix + ".disturbanceState",...
                       disturbanceGain, "outputType", "agent" + indAgent + "." + 'control');
                   feedforwardControl.PI_w.("agent" + indAgent) = PI_w;    
                end
                
                
            else
                PI_w = tool.solveGenericBvpFirstOrderCheb(S_r, obj.Lambda, quantity.Discrete.zeros(size(obj.Lambda), obj.domain),...
                    quantity.Discrete.zeros([size(obj.Lambda, 2), size(S_r, 1)], obj.domain), A0tilde, ...
                    obj.Q0, zeros(size(obj.E1, 2), size(S_r, 1)), obj.output.backstepping(obj.backsteppingKernel, "inverse", true), P_r);
                for indAgent = 1:obj.N
                   feedforwardControl.PI_w.("agent" + indAgent) = PI_w;    
                end
                               
                %Compute uDash for each agent
                uDash = {};
                for indAgent = 1:obj.N
%                    disturbanceGain = disturbanceMatAggregated((indAgent - 1) * size(G3dash, 1) + 1 : indAgent *  size(G3dash, 1), :);
                   feedforwardControl.disturbanceGain.("agent" + indAgent) = disturbanceGain;
                   uDash{indAgent, 1} = misc.Gain("agent" + indAgent + ".reference", eye(size(obj.E1, 2)), "outputType", "agent" + indAgent + "." + 'control') + ...
                       misc.Gain("agent" + indAgent + "." + obj.network.prefix + ".disturbanceState", disturbanceGain, "outputType", "agent" + indAgent + "." + 'control');
                end
            end
            
            %Compute u (controlInput) and save it in feedforwardControl
            for indAgent = 1:obj.N
                stable = model.Output("agent" + indAgent + "." + 'control', "input", uDash{indAgent});
                feedforwardControl.stable.("agent" + indAgent) = {'feedback', stable};
                feedforwardControl.input.("agent" + indAgent) = {'feedback', obj.backsteppingKernel.getFeedback(obj.Q1).series("input", stable)};
            end
            
            systemTarget = model.Transport(obj.Lambda, "A0", A0tilde, "Q0", obj.Q0, ...
                "output", obj.output.backstepping(obj.backsteppingKernel, "inverse", true), ...
                "input", model.Input("empty", "B", quantity.Discrete.zeros([obj.n, size(P, 1)], obj.domain)));
            
        end % feedforwardControl
        
        function u = getStabilizingStateFeedback(obj)
            u = obj.backsteppingKernel.getFeedback(obj.Q1);
        end
        
        function [PI_v, gain] = solveRegulatorEquationsDisturbance(obj, S, P)
            %Computes the solution of the regulator equations
            %   Lambda(z)dzPI_v(z) + A0Tilde(z)E1.'PI_v(0) - PI_v(z)S = -G1Tilde(z)*P
            %   E2.'PI_v(0) = Q0 E1-'PI_v(0) + G2*P
            %   outputBackstepped[PI_v] = -G4*P
            %using numeric methods. Also returns the resulting gain value
            %   gain = -1 * (obj.G3 * P - obj.E1.' * PI_v(1));
            
            S_1 = kron(S.', eye(obj.n));
            S_k = eye(size(S_1));
            
            S1_1 = kron(S.', eye(obj.output.lengthOutput));
            S1_k = eye(size(S1_1));
            
            poly1 = S_k*obj.coefRegEqD1{1};
            poly2 = S1_k*obj.coefRegEqD2{1};
            poly3 = S1_k*obj.coefRegEqD3{1};
            poly4 = S_k*obj.coefRegEqD4{1};
            G4vec = obj.G4 * P;
            G4vec = G4vec(:);
            for idx = 1:length(obj.coefRegEqD1)-1
                S_k = S_k * S_1;
                S1_k = S1_k * S1_1;
                poly1 = poly1 + S_k*obj.coefRegEqD1{idx+1};
                poly2 = poly2 + S1_k*obj.coefRegEqD2{idx+1};
                poly3 = poly3 + S1_k*obj.coefRegEqD3{idx+1};
                poly4 = poly4 + S_k*obj.coefRegEqD4{idx+1};
            end
            pi_v = -poly1/poly2*(G4vec + poly3*P(:)) + poly4*P(:);
            PI_v = reshape(pi_v, size(obj.E1, 1), []);
            
            gain = -1 * (obj.G3 * P - obj.E1.' * PI_v.at(1));
        end
        
        function [PI_w, gain] = solveRegulatorEquationsReference(obj, S_r, P_r)
            %Computes the solution of the regulator equations
            %   Lambda(z)dzPI_w(z) + A0Tilde(z)E1.'PI_w(0) - PI_w(z)S_r = 0
            %   E2.'PI_v(0) = Q0 E1.'PI_v(0)
            %   outputBackstepped[PI_w] = P_r
            %by computing a polynomial. Also returns the resulting gain
            %value 
            %   gain = E1.'*PI_w(1)
            
            S_1 = kron(S_r.', eye(obj.n));
            S_k = eye(size(S_1));
            poly1 = obj.coefRegEq1{1}*S_k*obj.polynomialSFactor;
            poly2 = obj.coefRegEq2{1}*S_k*obj.polynomialSFactor;
            for idx = 1:length(obj.coefRegEq1)-1
                S_k = S_k * S_1;
                poly1 = poly1 + obj.coefRegEq1{idx+1}*S_k*obj.polynomialSFactor;
                poly2 = poly2 + obj.coefRegEq2{idx+1}*S_k*obj.polynomialSFactor;
			end
			%Dilemma: pinv() makes more sense than inv(), but doesn`t give
			%a warning if poly2 is not invertible :/
            pi_w = poly1/poly2*P_r(:);
            PI_w = reshape(pi_w, size(obj.E1, 1), []);
            gain = obj.E1.'*PI_w.subs("z", 1);
        end
        
        function setRegulatorEquationRefCoefficients(obj, n_r)
            %Determins the coefficients of the polynoms to solve the
            %regulator equations for PI_w
            obj.coefRegEq1 = {};
            obj.coefRegEq2 = {};
            obj.polynomialSFactor = kron(eye(n_r), obj.E1 + obj.E2*obj.Q0);
            
            diagPhi = int(subs(inv(obj.Lambda), "z", "z_tilde"), "z_tilde", "zeta", "z");
            for idx = 0:obj.nCoefRef-1
                obj.coefRegEq1{idx+1} = kron(eye(n_r), (diagPhi.subs("zeta", 0)^idx ...
                    - int(diagPhi^idx/obj.Lambda.subs("z", "zeta")*obj.A0Backstepped.subs("z", "zeta")*obj.E1.', "zeta", 0, "z"))/factorial(idx));
                out_tilde = obj.outputBackstepped;
                Omega = int(subs(kron(eye(n_r), out_tilde.C) * obj.coefRegEq1{idx+1}, "z", "zeta"), "zeta", 0, 1) ...
                    + kron(eye(n_r), out_tilde.C0) * subs(obj.coefRegEq1{idx+1}, "z", 0) ...
                    + kron(eye(n_r), out_tilde.C1) * subs(obj.coefRegEq1{idx+1}, "z", 1);
                for indOut = 1:out_tilde.k
                    Omega = Omega + kron(eye(n_r), out_tilde.Ck(:, :, indOut)) * subs(obj.coefRegEq1{idx+1}, "z", out_tilde.zk(indOut));
                end
                obj.coefRegEq2{idx+1} = Omega;
            end
        end
        
        function setRegulatorEquationDistCoefficients(obj, n_d)
            %Determins the coefficients of the polynoms to solve the
            %regulator equations for PI_w
            obj.coefRegEqD1 = {};
            obj.coefRegEqD2 = {};
            obj.coefRegEqD3 = {};
            obj.coefRegEqD4 = {};
            
            out_tilde = obj.outputBackstepped;
            
            diagPhi = int(subs(inv(obj.Lambda), "z", "z_tilde"), "z_tilde", "zeta", "z");
            for idx = 0:obj.nCoefDist-1
                diagPhi_n = diagPhi^idx/factorial(idx);
                M = kron(eye(n_d), diagPhi_n.subs("zeta", 0)-int(diagPhi_n/obj.Lambda.subs("z", "zeta")*obj.A0Backstepped.subs("z", "zeta")*obj.E1.', "zeta", 0, "z"));
                M2 = M*kron(eye(n_d), obj.E2*obj.G2) - int(kron(eye(n_d), diagPhi_n/obj.Lambda.subs("z", "zeta")*obj.G1Backstepped.subs("z", "zeta")), "zeta", 0, "z");
                obj.coefRegEqD1{idx+1} = M*kron(eye(n_d), obj.E1 + obj.E2*obj.Q0);
                obj.coefRegEqD4{idx+1} = M2;
                Omega = int(subs(kron(eye(n_d), out_tilde.C) * M, "z", "zeta"), "zeta", 0, 1) ...
                    + kron(eye(n_d), out_tilde.C0) * M.at(0) ...
                    + kron(eye(n_d), out_tilde.C1) * M.at(1);
                M3 = int(subs(kron(eye(n_d), out_tilde.C) * M2, "z", "zeta"), "zeta", 0, 1) ...
                    + kron(eye(n_d), out_tilde.C0) * M2.at(0) ...
                    + kron(eye(n_d), out_tilde.C1) * M2.at(1);
                for indOut = 1:out_tilde.k
                    Omega = Omega + kron(eye(n_d), out_tilde.Ck(:, :, indOut)) * subs(M, "z", out_tilde.zk(indOut));
                    M3 = M3  + kron(eye(n_d), out_tilde.Ck(:, :, indOut)) * subs(M2, "z", out_tilde.zk(indOut));
                end
                Omega = Omega * kron(eye(n_d), obj.E1 + obj.E2*obj.Q0);
                obj.coefRegEqD2{idx+1} = Omega;
                obj.coefRegEqD3{idx+1} = M3;
            end
        end
        
        function plotState(obj, simData)
            %Creates a figure that displays the course of the system state
            %of each agent by plotting the data stored in
            %simData.(obj.agentPrefix).x for all elements of obj.agentPrefix
            figure('Position', [10 10 obj.N*200 300 + obj.output.lengthOutput*250]);
            ax = gobjects(fliplr([obj.N, obj.n]));
            for indAgent = 1:obj.N
               for indDim = 1:obj.n
                  ax((indAgent - 1)*obj.n + indDim) = subplot(obj.N, obj.n, (indAgent - 1)*obj.n + indDim);
                  simDomain = simData.("agent"+indAgent).x.domain;
                  surf(simDomain(2).grid, simDomain(1).grid, simData.("agent"+indAgent).x(indDim).valueDiscrete,...
                      'FaceAlpha', 0.5, 'edgeColor', 'interp');
                  title("x_" + indDim);
                  xlabel(simData.("agent"+indAgent).x(indDim).domain(2).name);
                  zlabel("x_" + indDim + "(" + simData.("agent"+indAgent).x(indDim).domain(1).name + ...
                      "," + simData.("agent"+indAgent).x(indDim).domain(2).name + ")");
                  ylabel(simData.("agent"+indAgent).x(indDim).domain(1).name);
               end
            end
            ax = ax';
            
            % Get upper position of each row of axes, normalized coordinates
            axPos = cell2mat(get(ax(:,1), 'Position'));
            axUpperPos = sum(axPos(:,[2,4]),2);  %upper pos.
            % Get center position for 1st row (assumes all rows have same center)
            axPos = cell2mat(get(ax(1,[1,end]),'Position'));
            axCenterPos = mean([axPos(1,1), sum(axPos(2,[1,3]))]);
            axPos = cell2mat(get(ax(:,1), 'Position'));
            titleHandles = gobjects(obj.N); 
            for indAgent = 1:obj.N
               %Increases space between subplots
               for indDim = 1:obj.n
                   pos = get(ax(indAgent, indDim), 'position');
                   set(ax(indAgent, indDim), 'Position', pos + [0 0.02 0 -0.04]);
               end
               titleHandles = annotation('textbox','String',"State of Agent " + indAgent, ...
                   'Position', [axCenterPos, axUpperPos(indAgent), 0, 0], ... 
                   'HorizontalAlignment', 'center','VerticalAlignment','bottom',...
                   'LineStyle','none','FitBoxToText','on', ...
                   'FontWeight',ax(1).Title.FontWeight, ... % matches title property
                   'FontSize', ax(1).Title.FontSize, ...    % matches title property
                   'FontName', ax(1).Title.FontName, ...    % matches title property
                   'Color', ax(1).Title.Color);             % matches title property
           end
        end
        
        function plotOutput(obj, simData)
           %Creates a figure that displays the output
           %of each agent by plotting the data stored in
           %simData.(obj.agentPrefix).(obj.output.type) for all elements of obj.agentPrefix
           figure('Position', [10 10 obj.N*200 300 + obj.output.lengthOutput*250]);
           ax = gobjects(fliplr([obj.N, obj.output.lengthOutput]));
           for indAgent = 1:obj.N
               for indDim = 1:obj.output.lengthOutput
                   ax((indAgent - 1)*obj.output.lengthOutput + indDim) = subplot(obj.N, obj.output.lengthOutput, (indAgent - 1)*obj.output.lengthOutput + indDim);
                   simDomain = simData.("agent"+indAgent).(obj.output.type).domain;
                   plot(simDomain.grid, simData.("agent"+indAgent).(obj.output.type)(indDim).valueDiscrete);
                   if isfield(simData.("agent"+indAgent), "reference")
                      hold on; 
                      plot(simDomain.grid, simData.("agent"+indAgent).reference(indDim).valueDiscrete, '--g');
                      hold off;
                   end
                   title("Output " + indDim);
                   xlabel(simData.("agent"+indAgent).(obj.output.type)(indDim).domain.name);
                   ylabel("y_" + indDim + "(" + simData.("agent"+indAgent).(obj.output.type)(indDim).domain.name + ")");
               end
           end
           ax = ax';
           
           % Get upper position of each row of axes, normalized coordinates
           axPos = cell2mat(get(ax(:,1), 'Position'));
           axUpperPos = sum(axPos(:,[2,4]),2);  %upper pos.
           % Get center position for 1st row (assumes all rows have same center)
           axPos = cell2mat(get(ax(1,[1,end]),'Position'));
           axCenterPos = mean([axPos(1,1), sum(axPos(2,[1,3]))]);
           axPos = cell2mat(get(ax(:,1), 'Position'));
           titleHandles = gobjects(obj.N); 
           for indAgent = 1:obj.N
               %Increases space between subplots
               for indDim = 1:obj.output.lengthOutput
                   pos = get(ax(indAgent, indDim), 'position');
                   set(ax(indAgent, indDim), 'Position', pos + [0 0.02 0 -0.04]);
               end
               titleHandles = annotation('textbox','String',"Output of Agent " + indAgent, ...
                   'Position', [axCenterPos, axUpperPos(indAgent), 0, 0], ... 
                   'HorizontalAlignment', 'center','VerticalAlignment','bottom',...
                   'LineStyle','none','FitBoxToText','on', ...
                   'FontWeight',ax(1).Title.FontWeight, ... % matches title property
                   'FontSize', ax(1).Title.FontSize, ...    % matches title property
                   'FontName', ax(1).Title.FontName, ...    % matches title property
                   'Color', ax(1).Title.Color);             % matches title property
           end
        end
        
        
        function [PI_w, PI_wPolynomCoefficients] = solveRegulatorEquationReferenceExplicit(obj, S_r, P_r)
            %solveRegulatorEquationReferenceExplicit solves the regulator
            %equations for a given S_r, P_r and nacksteppingKernel
            f_M = subs(inv(obj.Lambda)*obj.A0Backstepped*obj.E1.', "z", "zeta");
            
            S_1 = kron(S_r.', eye(size(obj.Lambda)));
            S_n = S_1^0;
            
            diagPhi_1 = int(subs(inv(obj.Lambda), "z", "z_tilde"), "z_tilde", "zeta", "z");
            diagPhi0_1 = subs(diagPhi_1, "zeta", 0);
            diagPhi0_n = diagPhi0_1^0;
            diagPhi_n = diagPhi_1^0;
            n_fak = 1;

            M = kron(S_r^0, diagPhi0_n - int(diagPhi_n * f_M, "zeta", 0, "z"));
            
            ceilPsi = ceil(max(max(diagPhi_1), [], 'all') * max(max(S_r), [], 'all'));
            ind = 1;
            PI_wPolynomCoefficients = {};
            while ceilPsi^ind/n_fak > 10^(-5)
                diagPhi0_n = diagPhi0_n * diagPhi0_1;
                diagPhi_n = diagPhi_n * diagPhi_1;
                n_fak = n_fak * ind;
                Q = kron(eye(size(S_r)), (diagPhi0_n - int(diagPhi_n * f_M, "zeta", 0, "z")) / n_fak);
                PI_wPolynomCoefficients{end+1} = Q;
                S_n = S_n * S_1;
                K = Q * S_n;
                M = M + K;
                ind = ind + 1;
            end
            M_tilde = M * kron(eye(size(S_r)), obj.E1 + obj.E2*obj.Q0);
            out_tilde = obj.outputBackstepped;
            Omega = int(subs(kron(eye(size(S_r)), out_tilde.C) * M_tilde, "z", "zeta"), "zeta", 0, "z") ...
                + kron(eye(size(S_r)), out_tilde.C0) * subs(M_tilde, "z", 0) ...
                + kron(eye(size(S_r)), out_tilde.C1) * subs(M_tilde, "z", 1);
            for indOut = 1:out_tilde.k
                Omega = Omega + out_tilde.Ck(:, :, indOut) * subs(M_tilde, "z", out_tilde.zk(indOut));
            end
            p_r = P_r(:);
            PI_w = M_tilde * inv(Omega) * p_r;
            PI_w = reshape(PI_w, size(obj.E1, 1), []);
        end
        
        %% get dependent properties
        function L = get.L(obj)
           L =  diag(sum(obj.adjacencyMatrix, 2)) - obj.adjacencyMatrix;
        end % get.L()
        
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
        
        function agent1IsLeader = get.agent1IsLeader(obj)
            agent1IsLeader = isempty(find(obj.adjacencyMatrix(1,:)));
        end % agent1IsLeader
        
        function n = get.n(obj)
           n = size(obj.Lambda, 1);
        end % get.n()
        
        function p = get.p(obj)
           p = sum(eig(obj.Lambda.at(0)) > 0);
        end % get.p()
        
        function m = get.m(obj)
           m = obj.n - obj.p;
        end % get.m()
        
        function N = get.N(obj)
           N = size(obj.adjacencyMatrix, 1);
        end % get.N()
        
        function domain = get.domain(obj)
           domain = obj.Lambda.domain;
        end % get.domain()

    end
    
    methods (Access = private)
       
        function ic = getNeutralInitialConditions(obj)
            ic = {};
            for agentPrefix = obj.network.agentPrefix.'
                ic{end + 1} = agentPrefix + ".x";
                ic{end + 1} = quantity.Symbolic(zeros(obj.n, 1), obj.domain);
            end          
        end
        
    end
end

