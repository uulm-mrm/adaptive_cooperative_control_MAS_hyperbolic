classdef DisturbanceObserver
    %DISTURBANCEOBSERVER cooperative disturbance observer for multi agent
    %systems with distributed parameters
    %
    %   Implements a disturbance observer for the disturbance model
    %
    %       vi_t(t) = S vi(t)
    %       d_i(t) = P vi(t)
    %   
    %   using the cooperative observer model
    %  
    %       ^vi_t = S ^vi + L_d ei^-(0)
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
    %       ni = xi^-(0)
    %       Gidash = Gi * P
    
    properties (SetAccess = protected)
        multiAgentSystem MultiAgent;
        L_d double;
        S;
        P;
        
        F;
        Omega;
        G1Tilde;
        MTilde;
        backsteppingKernel;
        
        L quantity.Discrete;
        M quantity.Discrete;
        observer struct;
    end
    
    methods
        function obj = DisturbanceObserver(mas, S, P, NameValue)
            %DISTURBANCEOBSERVER creates a disturbance observer for the
            %given multiAgentSystem and signal ode
            %
            %   Computes the parameters L(z) and M(z) of a
            %   disturbance observer for the multiAgentSystem
            arguments
                mas;
                S double;
                P double;
                NameValue.c_ = 1;
                NameValue.Q_ = eye(size(S));
                NameValue.R_ = eye(size(mas.E1, 2));
            end
            obj.multiAgentSystem = mas;
            obj.S = S; 
            obj.P = P;
            eigen = real(eig(mas.H));
            minC_ = 1/(2*min(eigen));
            if NameValue.c_ < minC_
                disp("Setting c_ to " + minC_ + ", which equals 1/2*min(real(eig(mas.H)))");
                NameValue.c = minC_;
            end
            
            %Solve Backstepping Kernel equations for inverse transformation
            %and compute backstepping Kernel
            P_I = DisturbanceObserver.solveObserverKernelEquations(mas.Lambda, mas.A, mas.Q1);
            backsteppingKernel = misc.Backstepping("kernelInverse", P_I, ...
                "signOfIntegralTermInverse", -1, ...
                "integralBounds", {0, "z"});
            F = -(mas.E1.' - mas.Q1 * mas.E2.') * P_I.subs("z", 1);
            obj.F = F;
            backsteppingKernel = backsteppingKernel.invert(+1);
            obj.backsteppingKernel = backsteppingKernel;
            
            %Compute transformed G1
            G1Tilde = backsteppingKernel.transform(-P_I.subs("zeta", 0)*mas.Lambda.subs("z", 0)*mas.E2*mas.G2*P + mas.G1*P);
            obj.G1Tilde = G1Tilde;
            %Solve decoupling equations
            Omega = tool.solveGenericBvpFirstOrderCheb(S, mas.Lambda, quantity.Discrete.zeros(size(mas.Lambda), mas.domain), ...
                -G1Tilde, quantity.Discrete.zeros([size(mas.Lambda, 1), size(mas.E1, 2)], mas.domain), 0, mas.G2*P, ...
                model.Output("", "C1", mas.E1.' - mas.Q1*mas.E2.', "C", F), mas.G3*P);
            obj.Omega = Omega;
            %Determine observer gain values
            %   Check if the pair (E1.' Omega(0), S) is observable
            if ~obj.checkObservability()
               warning("The pair (E1.' Omega(0), S) is not observable"); 
               L_d = zeros(size(S, 1), size(mas.E1, 2));
            else
                Qw = mas.E1.'*Omega.at(0);
                P_ = icare(S.', Qw.', NameValue.Q_, NameValue.R_, zeros(size(Qw.')), eye(size(S)), zeros(size(S)));
                L_d = NameValue.c_*P_*Qw.'/NameValue.R_;
            end
            MTilde = Omega*L_d;
            obj.MTilde = MTilde;
            obj.M = backsteppingKernel.transformInverse(MTilde);
            obj.L = -subs(P_I, "zeta", 0)*subs(mas.Lambda, "z", 0)*mas.E1;
            obj.L_d = L_d;
            
            %Define observers as model.Observer for simulation
            network = mas.network.copy();
            measurementOutput = model.Output("measurement", "C0", mas.E1.');
            network.addCommunication2agents(measurementOutput);
            
            observers = struct();
            for indAgent = 1:mas.N
               gainOde = misc.Gain("agent" + indAgent + ".error.measurement", L_d, "outputType", "agent" + indAgent + ".observer" + ".disturbanceModelInput");
               gainLocal = model.Input("agent" + indAgent + ".error.measurement", "B", obj.L);
               gainGlobal = model.Input("agent" + indAgent + ".network.error.measurement", "B", obj.M);
               dynamics = ss(S, eye(size(S)), P, [], "OutputName", "agent" + indAgent + ".observer" + ".disturbanceState", ...
                   "InputName", "agent" + indAgent + ".observer" + ".disturbanceModelInput");
               ode = stateSpace.Odes(dynamics, "agent" + indAgent + ".observer" + ".disturbanceState");
               
%                plant = network.agent(indAgent).copy;
               plant = copyAndReplace(network.agent(indAgent), "bc0", model.Output("plant.bc0", "C0", -mas.E2.'));
               plant.input.remove("disturbance");
               observerInput = model.Input("agent" + indAgent + ".network.observer.disturbanceState", ...
                   "B", mas.G1 * P, "B0", mas.G2 * P, "B1", mas.G3 * P);
               
               observer = model.Observer(plant, "gainOde", gainOde, "gain", gainLocal,...
                   "ode", ode, "prefixSignalModel", "agent" + indAgent, "prefix", "agent" + indAgent + ".observer", "keepBc1", true);
               observer.dynamics.input.add(gainGlobal);
               observer.dynamics.input.add(observerInput);
               measurementInput = model.Input("agent" + indAgent + ".measurement",...
                   "B0", mas.Q0);
               observer.dynamics.input.add(measurementInput);
               observers.("agent" + indAgent) = {'observer', observer};
            end
            obj.observer = observers;
            
        end
        function plotDisturbanceState(obj, simData)
           %Creates a figure that displays the output
           %of each agent by plotting the data stored in
           %simData.(obj.agentPrefix).(obj.output.type) for all elements of obj.agentPrefix
           figure('Position', [10 10 obj.multiAgentSystem.N*200 300 + size(obj.S, 1)*250]);
           ax = gobjects(fliplr([obj.multiAgentSystem.N, size(obj.S, 1)]));
           for indAgent = 1:obj.multiAgentSystem.N
               for indDim = 1:size(obj.S, 1)
                   ax((indAgent - 1)*size(obj.S, 1) + indDim) = subplot(obj.multiAgentSystem.N, size(obj.S, 1), (indAgent - 1)*size(obj.S, 1) + indDim);
                   simDomain = simData.("agent" + indAgent).observer.disturbanceState(1).domain;
                   plot(simDomain.grid, simData.("agent" + indAgent).observer.disturbanceState(indDim).valueDiscrete);
                   title("DisturbanceState. " + indDim);
                   xlabel(simData.("agent" + indAgent).observer.disturbanceState(indDim).domain.name);
                   ylabel("y_" + indDim + "(" + simData.("agent" + indAgent).observer.disturbanceState(indDim).domain.name + ")");
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
               for indDim = 1:size(obj.S, 1)
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
        
        function observable = checkObservability(obj)
            %checkObservability checks if the pair (E1.' Omega(0), S) is observable
            arguments
               obj;
            end
            
            invLambda = 1/obj.multiAgentSystem.Lambda;
            phi = int(subs(invLambda, "z", "zeta"), "zeta", 0, "z");
            [eigVec, eigVal] = eig(obj.S);
            
            observable = true;
            for indEv = 1:size(eigVec, 2)
                psi = expm(eigVal(indEv, indEv)*(phi - phi.subs("z", "zeta")));
                M1 = psi.subs("zeta", 0)*obj.multiAgentSystem.E2*obj.multiAgentSystem.G2*obj.P...
                    - int(psi*invLambda.subs("z", "zeta")*obj.G1Tilde.subs("z", "zeta"), "zeta", 0, "z");
                N = -(obj.multiAgentSystem.E1.' - obj.multiAgentSystem.Q1 * obj.multiAgentSystem.E2.')*M1.subs("z", 1) ...
                    -int(obj.F*M1, "z", 0, 1) + obj.multiAgentSystem.G3*obj.P;
                if N*eigVec(:, indEv) == 0
                   observable = false; 
                end
            end
        end
    end
    
    
    
    methods(Static)
        function P_I = solveObserverKernelEquations(Lambda, A, Q1)
             %Solves the observer kernel equations by transforming them
             %into the same form as the regulator kernel equations and
             %using the regulator kernel equations solver
             L_ = Lambda.flipDomain("z");
             A_ = L_*A.flipDomain("z").'/L_;
             dualBacksteppingKernel = kernel.Feedback.tryToLoadKernel(L_, "A", A_, "Q0", -Q1.');
             N = dualBacksteppingKernel.getValue();
             P_I_TildeT = L_^(-1) * subs(N, ["z", "zeta"], "zeta", "z") * subs(L_, "z", "zeta"); %
            
             P_I = flipDomain(P_I_TildeT.', ["z", "zeta"]);
        end
        
    end
end

