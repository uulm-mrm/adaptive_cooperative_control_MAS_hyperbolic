classdef ReferenceObserver
    %REFERENCEOBSERVER adaptive reference observer for multi
    %agent systems with distributed parameters
    %
    %   Implements an adaptive reference observer for the
    %   global reference model
    %
    %       w_t(t) = S_r w(t)
    %       r(t) = P_r w(t)
    
    properties
        nDimensions (1,1) double {mustBeInteger, mustBeNonnegative}; %Number of dimensions of w
		rDimensions (1,1) {mustBeInteger, mustBeNonnegative}; %Number of dimensions of r
        S_init (:,:,:) double; %Initial value of predicted S
        w_init (:,:) double; %Initial value of predicted state w
		P_init (:,:,:) double; %Initial value of predicted P
        L_local (:, :,:) double; %local observer gain of S for informed agents
        l_local (:, :,:) double; %local observer gain of w for informed agents
        L_global (1,1) double; %global observer gain of S for uninformed agents
        l_global (1,1) double; %global observer gain of w for uninformed agents
		Lp_local (:, :,:) double; %local observer gain of P for informed agents
        Lp_global (1,:) double; %local observer gain of P for informed agents
        
        adjacencyMatrix double; %adjacency matrix A of the communication Graph G
        
        agent1IsLeader (1,1) logical; %True if agent1 is the leader
    end
    
    properties (Dependent = true)
       H double;   %Leader-follower matrix
       Laplacian double;   %Laplacian matrix 
       N double {mustBeInteger, mustBeNonnegative}; %Number of agents
       nInformed double {mustBeInteger, mustBeNonnegative}; %Number of informed agents
       
    end
    
    methods
        function obj = ReferenceObserver(adjacencyMatrix, nDimensions, L_local, l_local, L_global, l_global, rDimensions, Lp_local, Lp_global, S_init, w_init, P_init)
            %REFERENCEOBSERVER Construct an instance of this class
            %   Detailed explanation goes here
            arguments
               adjacencyMatrix {mustBe.quadraticMatrix(adjacencyMatrix)} = [];
               nDimensions (1,1) {mustBeInteger, mustBeNonnegative} = 0;
               L_local (:, :,:) double = 1;
               l_local (:, :,:) double = 1;
               L_global (1,1) double = 1;
               l_global (1,1) double = 1;
			   rDimensions (1,1) {mustBeInteger, mustBeNonnegative} = nDimensions;
			   Lp_local (:, :,:) double = 1;
			   Lp_global (1,1) double = 1;
               S_init (:,:,:) double = zeros(nDimensions, nDimensions, size(adjacencyMatrix, 1));
               w_init (:,:) double = zeros(nDimensions, size(adjacencyMatrix, 1));
			   P_init (:,:,:) double = zeros(rDimensions, nDimensions, size(adjacencyMatrix, 1));
            end
            assert(size(S_init, 1) == size(S_init, 2), "Initial value of S must be quadratic");
            if isempty(find(adjacencyMatrix(1,:)))
                obj.agent1IsLeader = true;
                obj.adjacencyMatrix = adjacencyMatrix;
            else
                obj.agent1IsLeader = false;
                warning("Informed agents are not specified. Assuming agent 1 is informed.");
                obj.adjacencyMatrix = blkdiag(0, adjacencyMatrix);
                obj.adjacencyMatrix(2, 1) = 1;
                S_init = cat(3, zeros(nDimensions, nDimensions, 1), S_init);
				P_init = cat(3, zeros(rDimensions, nDimensions, 1), P_init);
                w_init = cat(2, zeros(nDimensions, 1), w_init);
            end
            
            obj.nDimensions = nDimensions;
			obj.rDimensions = rDimensions;
            assert(obj.nInformed == size(L_local, 3) || length(L_local) == 1, ...
                "Number of local observer gains must equal the number of informed agents or 1");
            if length(L_local) == 1
                obj.L_local = [];
				for idx = 1:sum(obj.adjacencyMatrix(:, 1)~=0)
					obj.L_local = cat(3, obj.L_local, eye(nDimensions)*L_local);
				end
            else
                obj.L_local = L_local;
			end
			assert(obj.nInformed == size(Lp_local, 3) || length(Lp_local) == 1, ...
                "Number of local observer gains must equal the number of informed agents or 1");
            if length(Lp_local) == 1
                obj.Lp_local = [];
				for idx = 1:sum(obj.adjacencyMatrix(:, 1)~=0)
					obj.Lp_local = cat(3, obj.Lp_local, eye(rDimensions)*Lp_local);
				end
            else
                obj.Lp_local = Lp_local;
            end
            assert(obj.nInformed == size(l_local, 3) || length(l_local) == 1, ...
                "Number of local observer gains must equal the number of informed agents or 1");
            if length(l_local) == 1
                obj.l_local = [];
				for idx = 1:sum(obj.adjacencyMatrix(:, 1)~=0)
					obj.l_local = cat(3, obj.l_local, eye(nDimensions)*l_local);
				end
            else
                obj.l_local = l_local;
            end
            obj.L_global = L_global;
            obj.l_global = l_global;
			obj.Lp_global = Lp_global;
            assert(size(S_init, 3) == obj.N || size(S_init, 3) == 1, ...
                "Number of inital values of S must equal the number of agents or 1");
            if length(l_local) == 1
                obj.S_init = ones(1, 1, obj.N).*S_init;
            else
                obj.S_init = S_init;
			end
			
			assert(size(P_init, 3) == obj.N || size(P_init, 3) == 1, ...
                "Number of inital values of P must equal the number of agents or 1");
            if length(l_local) == 1
                obj.P_init = ones(1, 1, obj.N).*P_init;
            else
                obj.P_init = P_init;
            end
            
            assert(size(w_init, 2) == obj.N || size(w_init, 2) == 1, ...
                "Number of inital values of w must equal the number of agents or 1");
            if length(l_local) == 1
                obj.w_init = ones(1, obj.N).*w_init;
            else
                obj.w_init = w_init;
            end
            
            assert(det(obj.H(2:end, 2:end)) ~= 0, "Graph must be connected");
            
        end
        
        function simData = simulate(obj, w_r, S_r, P_r)
            %simulate computes the observed S and w of each agent for a
            %S_r and w communicated by a leader.
            %   Simulates the course of the obvserved reference state w and
            %   reference system matrix S of each agent. S_r must be a 
            %   quadratic matrix of the expected dimension.
            %   w is an object of type quantity.Discrete. 
            %   The function returns a struct simData containing the
            %   predicted values of S_r and w of each agent as objects of 
            %   type quantity.Discrete with a domain matching that of w.
            arguments
               obj;
               w_r quantity.Discrete;
               S_r {mustBe.quadraticMatrix} = obj.S_init(:,:,1); 
			   P_r = obj.P_init(:,:,1);
            end
            assert(size(S_r, 1) == obj.nDimensions, "Dimension of S_r must match the oberserver`s expected dimension");
            assert(size(w_r, 1) == obj.nDimensions, "Dimension of w must match the oberserver`s expected dimension");
			assert(size(P_r, 1) == obj.rDimensions, "Dimension of P_r must match the oberserver`s expected dimension");
			assert(size(P_r, 2) == obj.nDimensions, "Dimension of P_r must match the oberserver`s expected dimension");
            simData.w_r = w_r;
            if isnumeric(S_r)
                S_r = S_r * quantity.Discrete.ones(1, w_r(1).domain);
			end
			if isnumeric(P_r)
                P_r = P_r * quantity.Discrete.ones(1, w_r(1).domain);
            end
            %Vectorizse initial conditions for the S and P observer. The ic of
            %agent1 is S_r if it is the leader
            ic = S_r(:).at(0);
            for indAgent = 2:obj.N
               S_i = obj.S_init(:,:,indAgent);
               ic = cat(1, ic, S_i(:)); 
			end
            %Compute observed S
            [t_S, S_pred] = ode45(@(t,S) SOde(t,S,obj,S_r), [w_r(1).domain.lower, w_r(1).domain.upper], ic); 
            %Set prediction of agent0 to S
            diagS_pred = S_r;
            simData.("agent0").S = S_r;
			
            %Vectorizse initial conditions for the S and P observer. The ic of
            %agent1 is S_r if it is the leader
			ic = P_r(:).at(0);
            for indAgent = 2:obj.N
               P_i = obj.P_init(:,:,indAgent);
               ic = cat(1, ic, P_i(:)); 
            end
			%Compute observed P
            [t_P, P_pred] = ode45(@(t,P) POde(t,P,obj,P_r), [w_r(1).domain.lower, w_r(1).domain.upper], ic); 
            %Set prediction of agent0 to P
            simData.("agent0").P = P_r;
            
            for indAgent = 2:obj.N
               %Transform the dimensions of predicted S to what is expected
               S_i = S_pred(:, (indAgent - 1) * obj.nDimensions^2 + 1 : indAgent * obj.nDimensions^2);
               S_i = reshape(S_i, [length(t_S), obj.nDimensions, obj.nDimensions]);
               S_i_disc = quantity.Discrete(interp1(t_S, S_i, w_r(1).domain.grid), w_r(1).domain);
               simData.("agent" + (indAgent - 1)).S = S_i_disc;
               %Build diag(S_i) for the w observer
               if isempty(diagS_pred)
                   diagS_pred = S_i_disc;
               else
                   diagS_pred = blkdiag(diagS_pred, S_i_disc);
			   end
			   %Transform the dimensions of predicted P to what is expected
			   P_i = P_pred(:, (indAgent - 1) * obj.nDimensions*obj.rDimensions + 1 : indAgent * obj.nDimensions*obj.rDimensions);
               P_i = reshape(P_i, [length(t_P), obj.rDimensions, obj.nDimensions]);
               P_i_disc = quantity.Discrete(interp1(t_P, P_i, w_r(1).domain.grid), w_r(1).domain);
               simData.("agent" + (indAgent - 1)).P = P_i_disc;
            end
            %Concatenate initial conditions for the w observer. The ic of
            %agent0 is the value of w at t=0 if it is the leader
            ic = w_r.at(0);
            for indAgent = 2:obj.N
               w_i = obj.w_init(:,indAgent);
               ic = cat(1, ic, w_i); 
            end
            %Compute observed w
            [t_w, w_pred] = ode45(@(t,w) wOde(t,w,obj,diagS_pred,w_r), [w_r(1).domain.lower, w_r(1).domain.upper], ic);
            %Set prediction of agent0 to w 
            simData.("agent0").w = w_r;
            for indAgent = 2:obj.N
               %Extract observed w for each agent
               w_i = w_pred(:, (indAgent - 1) * obj.nDimensions + 1 : indAgent * obj.nDimensions);
               simData.("agent" + (indAgent - 1)).w = quantity.Discrete(interp1(t_w, w_i, w_r(1).domain.grid), w_r(1).domain);
            end
        end
        
        function plotObservation(obj, simData)
           %Plots S and w stored in simData for each agent
            
           %Plot S
           figure('Position', [10 10 obj.N*200 300 + obj.nDimensions*250]);
           ax = gobjects(fliplr([obj.N, obj.nDimensions]));
           for indAgent = 0:obj.N - 1
               for indDim = 1:obj.nDimensions
                   ax((indAgent)*obj.nDimensions + indDim) = subplot(obj.N, obj.nDimensions, (indAgent)*obj.nDimensions + indDim);
                   simDomain = simData.("agent"+indAgent).S.domain;
                   plot(simDomain.grid, simData.("agent"+indAgent).w(indDim).valueDiscrete);
                   hold on;
                   plot(simDomain.grid, simData.w_r(indDim).valueDiscrete, '--');
                   title("w_" + indDim);
                   xlabel(simData.("agent"+indAgent).w(indDim).domain.name);
                   ylabel("w_" + indDim + "(" + simData.("agent"+indAgent).w(indDim).domain.name + ")");
                   hold off;
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
               for indDim = 1:obj.nDimensions
                   pos = get(ax(indAgent, indDim), 'position');
                   set(ax(indAgent, indDim), 'Position', pos + [0 0.02 0 -0.04]);
               end
               titleHandles = annotation('textbox','String',"Observed reference signal state w of Agent " + (indAgent - 1), ...
                   'Position', [axCenterPos, axUpperPos(indAgent), 0, 0], ... 
                   'HorizontalAlignment', 'center','VerticalAlignment','bottom',...
                   'LineStyle','none','FitBoxToText','on', ...
                   'FontWeight',ax(1).Title.FontWeight, ... % matches title property
                   'FontSize', ax(1).Title.FontSize, ...    % matches title property
                   'FontName', ax(1).Title.FontName, ...    % matches title property
                   'Color', ax(1).Title.Color);             % matches title property
           end
            
        end
        
        %% get dependent properties
        function Laplacian = get.Laplacian(obj)
           Laplacian =  diag(sum(obj.adjacencyMatrix, 2)) - obj.adjacencyMatrix;
        end % get.L()

        function H = get.H(obj)
            H = obj.Laplacian + diag(obj.adjacencyMatrix(:, 1));
            if obj.agent1IsLeader
               H = H(2:end, 2:end); 
            end
        end % get.H()
        
        function N = get.N(obj)
           N = size(obj.adjacencyMatrix, 1);
        end % get.N()
        
        function nInformed = get.nInformed(obj)
           nInformed = sum(obj.adjacencyMatrix(:,1)~=0);
        end % get.nInformed()
        
    end
    
    methods (Access = private)
        
        function dSdt = SOde(t, S, obj, S_leader)
            %This helper function transforms the differential equation for
            %the calculation of the observed S into a form that can be
            %solved with the ode45 solver. For this the equation is first
            %vectorized, then aggregated using the Laplacian matrix of the
            %communication graph.
            
            %vec(S_i_t) = L sum a_ij I_n^2 vec(S_j-S_i)
            %vec S = diag(diag(Li), diag(L)) * (-Laplacian kron I_n^2) vec(S_j-S_i) 
			L_loc1 = [];
			L_loc2 = [];
			a0 = obj.adjacencyMatrix(:, 1);
			a0 = a0(a0~=0);
			for idx = 1:obj.nInformed
				L = obj.L_local(:,:,idx);
				L_loc1 = [L_loc1; kron(eye(obj.nDimensions), L*a0(idx))];
				L_loc2 = blkdiag(L_loc2, kron(eye(obj.nDimensions), L*a0(idx)));
			end
			L_loc = [L_loc1, -L_loc2];
            L_glo = eye(obj.N - obj.nInformed-1) * obj.L_global;
			redLapl = obj.Laplacian(obj.Laplacian(:,1)==0, :);
            weightedLaplacian = L_glo * redLapl(2:end, :);
            gainMat = kron(-weightedLaplacian, eye(obj.nDimensions^2));
            L_loc = [L_loc, zeros(size(L_loc, 1), size(gainMat, 2)-size(L_loc, 2))];
            M = [zeros(obj.nDimensions^2, size(L_loc, 2)); L_loc; gainMat];
            dSdt = M*cat(1, S_leader(:).at(t), S(obj.nDimensions^2 + 1:end));
		end % dSdt()
		
		function dPdt = POde(t, P, obj, P_leader)
            %This helper function transforms the differential equation for
            %the calculation of the observed S into a form that can be
            %solved with the ode45 solver. For this the equation is first
            %vectorized, then aggregated using the Laplacian matrix of the
            %communication graph.
            
			L_loc1 = [];
			L_loc2 = [];
			a0 = obj.adjacencyMatrix(:, 1);
			a0 = a0(a0~=0);
			for idx = 1:obj.nInformed
				L = obj.Lp_local(:,:,idx);
				L_loc1 = [L_loc1; kron(eye(obj.nDimensions), L*a0(idx))];
				L_loc2 = blkdiag(L_loc2, kron(eye(obj.nDimensions), L*a0(idx)));
			end
            L_loc = [L_loc1, -L_loc2];
            L_glo = eye(obj.N - obj.nInformed-1) * obj.Lp_global;
            redLapl = obj.Laplacian(obj.Laplacian(:,1)==0, :);
            weightedLaplacian = L_glo * redLapl(2:end, :);
            gainMat = kron(-weightedLaplacian, eye(obj.nDimensions*obj.rDimensions));
            L_loc = [L_loc, zeros(size(L_loc, 1), size(gainMat, 2)-size(L_loc, 2))];
            M = [zeros(obj.nDimensions*obj.rDimensions, size(L_loc, 2)); L_loc; gainMat];
            dPdt = M*cat(1, P_leader(:).at(t), P(obj.nDimensions*obj.rDimensions + 1:end));
        end % dPdt()
        
        function dwdt = wOde(t, w, obj, diagS_pred, w_leader)
            %This helper function transforms the differential equation for
            %the calculation of the observed w into a form that can be
            %solved with the ode45 solver. For this the observed S and 
            
            %w_ti = S_i w_i + l I_n sum a_ij (wj - wi)
            %w_t = diag(S_i) w + diag(diag(li), diag(l)) * (-Laplacian kron I_n) w
			L_loc1 = [];
			L_loc2 = [];
			a0 = obj.adjacencyMatrix(:, 1);
			a0 = a0(a0~=0);
			for idx = 1:obj.nInformed
				L = obj.l_local(:,:,idx);
				L_loc1 = [L_loc1; L*a0(idx)];
				L_loc2 = blkdiag(L_loc2, L*a0(idx));
			end
            L_loc = [L_loc1, -L_loc2];
            L_glo = eye(obj.N - obj.nInformed-1) * obj.l_global;
            redLapl = obj.Laplacian(obj.Laplacian(:,1)==0, :);
            weightedLaplacian = L_glo * redLapl(2:end, :);
            gainMat = kron(-weightedLaplacian, eye(obj.nDimensions));
            L_loc = [L_loc, zeros(size(L_loc, 1), size(gainMat, 2)-size(L_loc, 2))];
            M = diagS_pred.at(t) + [zeros(obj.nDimensions, size(L_loc, 2)); L_loc; gainMat];
            dwdt = M * cat(1, w_leader.at(t), w(obj.nDimensions + 1:end));
        end % dwdt()
        
    end
end

