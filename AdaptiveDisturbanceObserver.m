classdef AdaptiveDisturbanceObserver
    %ADAPTIVEDISTURBANVEOBSERVER Summary of this class goes here
    %systems with distributed parameters
    %
    %   Implements a disturbance observer for the disturbance model
    %
    %       vi_t(t) = S vi(t)
    %       d_i(t) = P vi(t)
    %   
    %   using the cooperative observer model
    %  
    %       ^Si_t = mu (sum_j=1^N a_ij(s_ij ^Si - ^Sj) + a_i0 ^Si)
    %       ^vi_t = ^S ^vi + L_d ei^+(1)
    %       ^xi_t = Lambda(z) ^xi_z(z,t) + A(z) ^xi(z,t) +
    %               G1dash(z)(sum_j=1^N a_ij(s_ij ^vi - ^vj) + a_i0 ^vi) +
    %               L(z)ei^-(0) + M(z)(sum_j=1^N a_ij(s_ij ei^-(0) - ej^-(0)) + a_i0 ei^-(0))
    %       ^xi^+(0) = Q0 ni + G2dash(sum_j=1^N a_ij(s_ij ^vi - ^vj) + a_i0 ^vi)
    %       ^xi^-(1) = Q1 ^xi^+(1) + u_i + G3dash(sum_j=1^N a_ij(s_ij ^vi - ^vj) + a_i0 ^vi)
    %
    %   with the aggregated system
    %
    %       ^v_t = (I kron S)^v + (I kron L_d) e^-(0)
    %       ^x_t = (I kron Lambda(z))^x_z + (I kron A(z))^x + (H kron G1dash)^v
    %                                     + (I kron L(z))e^-(0) + (H kron M(z))e^-(0)
    %       ^x^+(0) = (I kron Q0)n + (H kron G2dash)^v
    %       ^x^-(1) = (I kron Q1)^x^+(1) + (H kron G3dash)^v
    %
    %    where
    %       ei^-(0) = ni - ^xi^-(0)
    %       ni = xi^+(1)
    %       Gidash = Gi * P
    
    properties
        multiAgentSystem MultiAgent;
        n_d (1,1) double {mustBeInteger, mustBeNonnegative}; %Number of dimensions of v
		n_r (1,1) double {mustBeInteger, mustBeNonnegative}; %Number of dimensions of r
        
        F;
        G1Tilde;
        P_I;
        backsteppingKernel;
        
        mu (1,1) double;
		mu_P (1,1) double;
        c (1,1) double;
        observer struct;
        
        %Properties for the solution of the decoupling equations
        EpsPolynomCoefficients1 cell;
        EpsPolynomCoefficients2 cell;
        EpsPolynomCoefficients3 cell;
        phi_z1 cell;
        
        targetEvOde (1, :) double; % Eigenvalues of the Ode of the target system
    end
    
    methods
        function obj = AdaptiveDisturbanceObserver(multiAgentSystem, n_d, n_r, mu, NameValue)
            %ADAPTIVEDISTURBANCEOBSERVER creates a disturbance observer for the
            %given multiAgentSystem and signal ode
            %
            %   Computes the parameters L(z) and M(z) of a
            %   disturbance observer for the multiAgentSystem
            arguments
                multiAgentSystem;
                n_d (1,1) double {mustBeInteger, mustBeNonnegative};
                n_r double;
                mu double = 1;
                NameValue.nCoef (1,1) double {mustBeInteger, mustBeNonnegative} = 10;
                NameValue.targetEvOde (1, :) double = -1*(1:n_d);
				NameValue.mu_P (1,1) double = 1;
            end
            obj.n_d = n_d;
            obj.multiAgentSystem = multiAgentSystem;
            obj.n_r = n_r;
            obj.mu = mu;
			obj.mu_P = NameValue.mu_P;
            obj.targetEvOde = NameValue.targetEvOde;
            if mu <= 0
               warning("To stabilize the observer ODE mu must be positive") 
			end
			if NameValue.mu_P <= 0
               warning("To stabilize the observer ODE mu_P must be positive") 
            end
            
            %Solve Backstepping Kernel equations for inverse transformation
            %and compute backstepping Kernel
            P_I = AdaptiveDisturbanceObserver.solveObserverKernelEquations(multiAgentSystem.Lambda, multiAgentSystem.A, multiAgentSystem.Q0);
            obj.P_I = P_I;
            backsteppingKernel = misc.Backstepping("kernelInverse", P_I, ...
                "signOfIntegralTermInverse", -1, ...
                "integralBounds", {"z", 1});
            F = (multiAgentSystem.Q0*multiAgentSystem.E1.' -  multiAgentSystem.E2.') * P_I.subs("z", 0);
            obj.F = F;
            backsteppingKernel = backsteppingKernel.invert(+1);
            obj.backsteppingKernel = backsteppingKernel;
            
            %Compute transformed G1
            G1Tilde = backsteppingKernel.transform(P_I.subs("zeta", 1)*multiAgentSystem.Lambda.subs("z", 1)*multiAgentSystem.E1*multiAgentSystem.G3 + multiAgentSystem.G1);
            obj.G1Tilde = G1Tilde;

            %Determine polynom coefficients
            I = eye(obj.n_d);
            diagPhi = int(subs(inv(obj.multiAgentSystem.Lambda), "z", "z_tilde"), "z_tilde", "zeta", "z");
            b = kron(I, multiAgentSystem.Lambda^(-1)*G1Tilde);
            obj.EpsPolynomCoefficients1 = {};
            obj.EpsPolynomCoefficients2 = {};
            obj.EpsPolynomCoefficients3 = {};
            obj.phi_z1 = {};
            for idx = 0:NameValue.nCoef-1
                obj.EpsPolynomCoefficients1{idx+1} = int(kron(I, diagPhi^idx/factorial(idx))*b.subs("z", "zeta"), "zeta", "z", 1);
                obj.EpsPolynomCoefficients2{idx+1} = kron(I, (multiAgentSystem.Q0*multiAgentSystem.E1.'-multiAgentSystem.E2.')*diagPhi.subs(["z", "zeta"], 0, 1)^idx/factorial(idx)*multiAgentSystem.E2...
                    - int(obj.F*diagPhi.subs(["z", "zeta"], "zeta", 1)^idx/factorial(idx)*multiAgentSystem.E2, "zeta", 0, 1));
                obj.EpsPolynomCoefficients3{idx+1} = kron(I, (multiAgentSystem.Q0*multiAgentSystem.E1.'-multiAgentSystem.E2.')*diagPhi.subs(["z", "zeta"], 0, 1)^idx/factorial(idx)*multiAgentSystem.E1*obj.multiAgentSystem.G3)...
                    + kron(I, (multiAgentSystem.Q0*multiAgentSystem.E1.'-multiAgentSystem.E2.'))*int(kron(I, diagPhi.subs("z", 0)^idx/factorial(idx))*b.subs("z", "zeta"), "zeta", 0, 1)...
                    - int(kron(I, obj.F*diagPhi.subs(["z", "zeta"], "zeta", 1)^idx/factorial(idx)*multiAgentSystem.E1*multiAgentSystem.G3) ...
                    + int(kron(I, obj.F*diagPhi.subs(["z", "zeta"], "zeta", "z_Tilde")^idx/factorial(idx))*b.subs("z", "z_Tilde"), "z_Tilde", "zeta", 1), "zeta", 0, 1);
                obj.phi_z1{idx+1} = kron(I, diagPhi.subs("zeta", 1)^idx/factorial(idx));
            end
            
        end
        
        function stateSpace = getStateSpaceApproximation(obj, L, NameValue)
            arguments
                obj;
                L;
                NameValue.t = linspace(0, 1, 201);
                NameValue.prefix = "observer";
                NameValue.agent = 1;
            end
            %Assuming a homogenous system it would be possible to only get
            %the closed loop simulation model once and just changing output
            %names of the state space
            dynamics = copyAndReplace(obj.multiAgentSystem.network.agent(1), "bc1", model.Output("plant.bc1", "C1", - obj.multiAgentSystem.E1.'), ...
                "prefix", "agent" + NameValue.agent + "." + NameValue.prefix);
            gainLocal = model.Input("agent" + NameValue.agent + ".error.measurement", "B", L);
            measurementInput = model.Input("agent" + NameValue.agent +".measurement", "B1", obj.multiAgentSystem.Q1);
            controlInput = model.Input("agent" + NameValue.agent + ".control", "B1", eye(obj.multiAgentSystem.p));
            disturbanceInput = model.Input("agent" + NameValue.agent + ".disturbance", "B", obj.multiAgentSystem.G1,...
                    "B0", obj.multiAgentSystem.G2, "B1", obj.multiAgentSystem.G3);
            observerInput = gainLocal + measurementInput + controlInput + disturbanceInput;
            dynamics = copyAndReplace(dynamics, "input", observerInput, "bc1", model.Output("bc1", "C1",...
                - (obj.multiAgentSystem.E1).'));
            stateSpace = dynamics.getClosedLoopSimulationModel('t', NameValue.t);
		end
        
        function Eps = solveDecouplingEquations(obj, S, P)
            %SOLVEDECOUPLINGEQUATIONS calculates the solutions and
            %Epsilon(z) of the observer decoupling equations
            %   Lambda(z) dzEpsilon(z) - Epsilon(z)S = -G1Tilde(z)
            %   (Q0 E1.'-E2.') Epsilon(0) = int_0^1 F(zeta)Epsilon(zeta)dzeta - G2*P
            %   E1.'Epsilon(1) = G3*P
            
%             Omega = tool.solveGenericBvpFirstOrderCheb(S, obj.multiAgentSystem.Lambda, quantity.Discrete.zeros(size(obj.multiAgentSystem.Lambda), obj.multiAgentSystem.domain), ...
%                 -obj.G1Tilde, quantity.Discrete.zeros([size(obj.multiAgentSystem.Lambda, 1), size(obj.multiAgentSystem.E1, 2)], obj.multiAgentSystem.domain), 0, obj.multiAgentSystem.("G"+(2+obj.measurement))*obj.P, ...
%                 model.Output("", "C1", obj.multiAgentSystem.("E" + (2-obj.measurement)).' - obj.multiAgentSystem.Q1*obj.multiAgentSystem.("E" + (1+obj.measurement)).', "C", obj.F), obj.multiAgentSystem.("G"+(3-obj.measurement))*obj.P);
%             Eps = tool.solveGenericBvpFirstOrderCheb(S, -obj.multiAgentSystem.Lambda.flipDomain("z"), quantity.Discrete.zeros(size(obj.multiAgentSystem.Lambda), obj.multiAgentSystem.domain), ...
%                 -obj.G1Tilde.flipDomain("z"), quantity.Discrete.zeros([size(obj.multiAgentSystem.Lambda, 1), size(obj.multiAgentSystem.E1, 2)], obj.multiAgentSystem.domain),...
%                 0, zeros(size(obj.multiAgentSystem.("G"+(2+obj.measurement))*obj.P)), ...
%                 model.Output("", "C1", obj.multiAgentSystem.("E" + (2-obj.measurement)).' - obj.multiAgentSystem.Q1*obj.multiAgentSystem.("E" + (1+obj.measurement)).', "C", obj.F),...
%                 zeros(size(obj.multiAgentSystem.("G"+(3-obj.measurement))*obj.P)));
            
            %Set initial values
            S1_1 = kron(S.', eye(size(obj.multiAgentSystem.E2.', 1))); %ep
            SAgg1 = eye(size(S1_1));
            S_1 = kron(S.', eye(obj.multiAgentSystem.n));
            SAgg = eye(size(S_1));
            IE2 = kron(eye(obj.n_d), obj.multiAgentSystem.E2);
            M0 = SAgg*obj.EpsPolynomCoefficients1{1};
            M1 = SAgg1*obj.EpsPolynomCoefficients2{1};
            M2 = -SAgg1*obj.EpsPolynomCoefficients3{1};
            psi_z1 = SAgg*obj.phi_z1{1};
            %Calculate the polynomial
            for idx = 1:length(obj.EpsPolynomCoefficients1)-1
                SAgg = SAgg*S_1;
                SAgg1 = SAgg1 * S1_1;
                M0 = M0 + SAgg*obj.EpsPolynomCoefficients1{idx+1};
                M1 = M1 + SAgg1*obj.EpsPolynomCoefficients2{idx+1};
                M2 = M2 -SAgg1*obj.EpsPolynomCoefficients3{idx+1};
                psi_z1 = psi_z1 + SAgg*obj.phi_z1{idx+1};
            end
            G3vec = obj.multiAgentSystem.G3*P;
            G3vec = G3vec(:);
			G2vec = obj.multiAgentSystem.G2*P;
            G2vec = G2vec(:);
            eps = psi_z1*(kron(eye(obj.n_d), obj.multiAgentSystem.E1)*G3vec + IE2/M1*(M2*P(:) - G2vec)) + M0*P(:);
            Eps = reshape(eps, obj.multiAgentSystem.n, size(S, 1));
        end
        
        function [L_v, L] = getObserverGains(obj, S, P, method, NameValue)
            arguments
                obj;
                S;
				P;
				method = "riccati";
				NameValue.Q = 1;
				NameValue.R = 1;
				NameValue.c = 1;
            end
           %Solve decoupling equations for S_init
           [Eps] = obj.solveDecouplingEquations(S, P);
           %Check if system is controllable and return zeros in case it is not
           if rank(ctrb(S.', (obj.multiAgentSystem.E2.'*Eps.at(1)).')) ~= length(S)
              L = zeros(obj.multiAgentSystem.n, obj.multiAgentSystem.m);
              L_v = zeros(obj.n_d, obj.multiAgentSystem.m);
           else
               %Calculate L_v
			   if strcmp(method, "place")
					L_v = (place(S.', (obj.multiAgentSystem.E2.'*Eps.at(1)).', obj.targetEvOde)).';  
			   elseif strcmp(method, "riccati")
					B_ = (obj.multiAgentSystem.E2.'*Eps.at(1)).';
					P_ = icare(S.', B_, NameValue.Q, NameValue.R, zeros(size(B_)), eye(size(S)), zeros(size(S)));
					if isempty(P_)
						warning("Nearly unobservable, setting L_v to zero")
						L_v = zeros(obj.n_d, obj.multiAgentSystem.m);
					else
						L_v = NameValue.c*P_*B_/NameValue.R;
					end
			   else
				   error("Supported methods to determine the observer gain L_v are 'place' and 'riccati'.")
			   end         

               %Calculate L
               LTilde = Eps*L_v;
               L = obj.backsteppingKernel.transformInverse(LTilde) + obj.P_I.subs("zeta", 1)*obj.multiAgentSystem.Lambda.subs("z", 1)*obj.multiAgentSystem.E2;
           end
        end

		function observable = checkObservability(obj, S, P)
			invLambda = 1/obj.multiAgentSystem.Lambda;
            phi = int(subs(invLambda, "z", "zeta_"), "zeta_", "z", "zeta");
            [eigVec, eigVal] = eig(S);
            
            observable = true;
            for indEv = 1:size(eigVec, 2)
                psi = expm(-1*eigVal(indEv, indEv)*(phi));
                M1 = psi.subs("zeta", 1)*obj.multiAgentSystem.E1*obj.multiAgentSystem.G3*P...
                    + int(psi*invLambda.subs("z", "zeta")*obj.G1Tilde.subs("z", "zeta")*P, "zeta", "z", 1);
                N = (obj.multiAgentSystem.E2.' - obj.multiAgentSystem.Q0 * obj.multiAgentSystem.E1.')*M1.subs("z", 0) ...
                    -int(obj.F*M1.subs("z", "zeta"), "zeta", 0, 1) + obj.multiAgentSystem.G2*P;
                if N*(P*eigVec(:, indEv)) == 0
                   observable = false; 
                end
            end
		end
        
        function plotDisturbanceState(obj, simData)
           %Creates a figure that displays the output
           %of each agent by plotting the data stored in
           %simData.(obj.agentPrefix).(obj.output.type) for all elements of obj.agentPrefix
           figure('Position', [10 10 obj.multiAgentSystem.N*200 300 + obj.n_d*250]);
           ax = gobjects(fliplr([obj.multiAgentSystem.N, obj.n_d]));
           for indAgent = 1:obj.multiAgentSystem.N
               for indDim = 1:obj.n_d
                   ax((indAgent - 1)*obj.n_d + indDim) = subplot(obj.multiAgentSystem.N, obj.n_d, (indAgent - 1)*obj.n_d + indDim);
                   simDomain = simData.("agent" + indAgent).observer.x(1).domain;
%                    plot(simDomain.grid, simData.("agent" + indAgent).observer.x(indDim).valueDiscrete);
                   surf(simDomain(2).grid, simDomain(1).grid, simData.("agent"+indAgent).observer.x(indDim).valueDiscrete,...
                      'FaceAlpha', 0.5, 'edgeColor', 'none');
                   title("x" + indDim);
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
           titleHandles = gobjects(obj.multiAgentSystem.N); 
           for indAgent = 1:obj.multiAgentSystem.N
               %Increases space between subplots
               for indDim = 1:obj.n_d
                   pos = get(ax(indAgent, indDim), 'position');
                   set(ax(indAgent, indDim), 'Position', pos + [0 0.02 0 -0.04]);
               end
               titleHandles = annotation('textbox','String',"Observed Disturbance State of Agent " + indAgent, ...
                   'Position', [axCenterPos, axUpperPos(indAgent), 0, 0], ... 
                   'HorizontalAlignment', 'center','VerticalAlignment','bottom',...
                   'LineStyle','none','FitBoxToText','on', ...
                   'FontWeight',ax(1).Title.FontWeight, ... % matches title property
                   'FontSize', ax(1).Title.FontSize, ...    % matches title property
                   'FontName', ax(1).Title.FontName, ...    % matches title property
                   'Color', ax(1).Title.Color);             % matches title property
           end
        end
    end
    
    methods(Static)
        function P_I = solveObserverKernelEquations(Lambda, A, Q0)
             %Solves the observer kernel equations by transforming them
             %into the same form as the regulator kernel equations and
             %using the regulator kernel equations solver. See
             %hyperbolic/network.observer.backstepping.getDualKernel
			 
			backsteppingKernel = kernel.Observer.tryToLoadKernel("kernel.Observer", Lambda, "A", A, "Q0", Q0);
			P_I = backsteppingKernel.getValue();
		end
        
    end
end

