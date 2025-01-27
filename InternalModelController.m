classdef InternalModelController
	%INTERNALMODELCONTROLLER Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
        Lambda quantity.Discrete {mustBe.quadraticMatrix};
        A quantity.Discrete {mustBe.quadraticMatrix};
        Q0 double;
        Q1 double;
        output (1,1) model.Output;
		E1 double;
		E2 double;

		b cell;
		B (:,:) double;
		mu (1,1) double;
		L(:,:,:) double;

        adjacencyMatrix double; %adjacency matrix A of the communication Graph G
        
        backsteppingKernel;
        A0Backstepped;
        outputBackstepped;

		%Properties for the solution of the decoupling equations for
		%constant Lambda
        DecPolynomCoefficients1 cell;
        DecPolynomCoefficients2 cell;
        DecPolynomCoefficients3 cell;
        DecPolynomCoefficients4 cell;
        phi_z1 cell;

		%Properties for the solution of the decoupling equations for
		%non-constant Lambda
        phi1 quantity.Discrete;
		phi2 quantity.Discrete;
		nCoef;
		ih;
	end

	properties (Dependent = true)
       H double;   %Leader-follower matrix
       Laplacian double;   %Laplacian matrix 
       N double {mustBeInteger, mustBeNonnegative}; %Number of agents
            
        n (1,1) {mustBeInteger, mustBeNonnegative};	% number of PDEs, i.e. size(obj.Lambda, 1)
		p (1,1) {mustBeInteger, mustBeNonnegative};	% number of PDEs with positive lambda_i
		m (1,1) {mustBeInteger, mustBeNonnegative};	% number of PDEs with negative lambda_i

        n_w (1,1) {mustBeInteger, mustBeNonnegative};	% number of dimensions of the internal model
    end
	
	methods
		function obj = InternalModelController(adjacencyMatrix, Lambda, b, mu, L, NameValue)
			%INTERNALMODELCONTROLLER Construct an instance of
			%InternalModelController. This class implements an output
			%feedback controller based on the internal model principle.
			%   Detailed explanation goes here
			arguments
				adjacencyMatrix {mustBe.quadraticMatrix(adjacencyMatrix)} = [];
                Lambda {mustBe.quadraticMatrix(Lambda)} = [];
				b cell = {};
				mu (1,1) double = 1;
				L (:,:,:) double = 1;
                NameValue.A {mustBe.quadraticMatrix(NameValue.A)} = quantity.Discrete.zeros(size(Lambda), Lambda.getDomain());
                NameValue.Q0 double = zeros(length(find(Lambda.at(0)<0)), length(find(Lambda.at(0)>0)));
                NameValue.Q1 double = zeros(length(find(Lambda.at(0)>0)), length(find(Lambda.at(0)<0)));
                NameValue.output (1,1) model.Output = model.Output('measurement', "C", quantity.Discrete.zeros(length(Lambda), Lambda.getDomain()));
				NameValue.nCoef = 5;
			end % arguments
			obj.adjacencyMatrix = adjacencyMatrix;
            obj.Lambda = Lambda;
            obj.A = NameValue.A;
            obj.Q0 = NameValue.Q0;
            obj.Q1 = NameValue.Q1;
			obj.b = b;
			obj.B = blkdiag(b{:});
			obj.adjacencyMatrix = adjacencyMatrix;
			obj.mu = mu;
            obj.output = NameValue.output;
            if size(L,3) == 1
                obj.L = [];
				for idx = 1:sum(obj.adjacencyMatrix(:, 1)~=0)
					obj.L = cat(3, obj.L, L);
				end
            else
                obj.L = L;
			end

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
			
			obj.backsteppingKernel = kernel.Feedback.tryToLoadKernel("kernel.Feedback", obj.Lambda, "A", obj.A, "Q0", obj.Q0);
			obj.A0Backstepped = obj.backsteppingKernel.getA0Target();
            obj.outputBackstepped = obj.output.backstepping(obj.backsteppingKernel, "inverse", true);

        	ih = -obj.B*obj.outputBackstepped.C;
			ih = ih(:);
			ih = ih.subs("z", "zeta");
			%Determine polynom coefficients
			n_out = length(find(Lambda.at(0) > 0));
			n_w = size(obj.L, 1);
        	diagPhi = int(subs(kron(inv(Lambda), eye(n_out*n_w)), "z", "z_tilde"), "z_tilde", "z", "zeta");
        	obj.DecPolynomCoefficients1 = {};
        	obj.DecPolynomCoefficients2 = {};
        	obj.DecPolynomCoefficients3 = {};
        	obj.DecPolynomCoefficients4 = {};
        	obj.phi_z1 = {};
			C1 = kron((obj.E1.'+obj.Q0.'*obj.E2.'), eye(n_out*n_w));
			D = kron(obj.A0Backstepped.subs("z", "zeta").'/obj.Lambda.subs("z", "zeta"), eye(n_out*n_w));
        	for idx = 0:NameValue.nCoef-1
				diagPhi_k = diagPhi^idx/factorial(idx);
				obj.DecPolynomCoefficients1{idx+1} = diagPhi_k.subs("zeta", 1);
            	obj.DecPolynomCoefficients2{idx+1} = C1*diagPhi_k.subs(["z", "zeta"], 0, 1)...
                	- int(D*diagPhi_k.subs(["z", "zeta"], "zeta", 1), "zeta", 0, 1);
            	obj.DecPolynomCoefficients3{idx+1} = int(C1*diagPhi_k.subs("z", 0) * ih - ...
					D*int(diagPhi_k.subs(["z", "zeta"], "zeta", "zeta2")*ih.subs("zeta", "zeta2"), "zeta2", "zeta", 1), "zeta", 0, 1);
            	obj.DecPolynomCoefficients4{idx+1} = int(diagPhi^idx/factorial(idx)*ih, "zeta", "z", 1);
			end
		end

		function u = getStabilizingStateFeedback(obj)
            u = obj.backsteppingKernel.getFeedback(obj.Q1);
        end
		
		function [K, Eps] = getControlGains(obj, S, Q, R, c, limit)
            arguments
                obj;
                S;
				Q = 1;
				R = 1;
				c = 1;
				limit = max(abs(S), [], 'all')*100;
            end
           %Solve decoupling equations for S_init
           [Eps] = obj.solveDecouplingEquations(S);
		   S_ = kron(eye(length(find(obj.Lambda.at(0) > 0))), S);
		   B_ = Eps.at(1)*obj.Lambda.at(1)*obj.E1;
		   P_ = icare(S_, B_, Q, R, zeros(size(B_)), eye(size(S_)), zeros(size(S_)));
		   if isempty(P_) || any(max(P_)>limit, 'all')
			   P_ = zeros(size(B_, 1));
		   end
		   K = -c/R*B_.'*P_;
			
		end

		function Eps = solveDecouplingEquations(obj, S)
            %SOLVEDECOUPLINGEQUATIONS calculates the solutions and
            %Epsilon(z) of the observer decoupling equations
            %   dz(Lambda(z) Epsilon(z)) + Epsilon(z)S = -By Ctilde(z)
            %   Epsilon(0)Lambda(0)(E1-E2 Q0)  = int_0^1 Epsilon(zeta)A_0(zeta)dzeta
            %   Epsilon(1)E2.'Lambda(1) = 0
			S_ = kron(eye(obj.n*obj.p), S);
			S__ = kron(eye(obj.p^2), S);
			b1 = obj.B*obj.outputBackstepped.C1*obj.E2;
			b1 = b1(:);
			b0 = obj.B*(obj.outputBackstepped.C0*obj.E1 + obj.outputBackstepped.C0*obj.E2*obj.Q0);
			b0 = b0(:);
        	S_1 = eye(size(S_));
        	S__1 = eye(size(S__));
        	M1 = S_1*obj.DecPolynomCoefficients1{1};
        	M2 = S__1*obj.DecPolynomCoefficients2{1};
        	M3 = S__1*obj.DecPolynomCoefficients3{1};
        	M4 = S_1*obj.DecPolynomCoefficients4{1};
        	%Calculate the polynomial
        	for idx = 1:length(obj.DecPolynomCoefficients1)-1
            	S_1 = S_1*S_;
            	S__1 = S__1*S__;
            	M1 = M1 + S_1*obj.DecPolynomCoefficients1{idx+1};
            	M2 = M2 + S__1*obj.DecPolynomCoefficients2{idx+1};
            	M3 = M3 + S__1*obj.DecPolynomCoefficients3{idx+1};
            	M4 = M4 + S_1*obj.DecPolynomCoefficients4{idx+1};
			end
			M3 = M3 - b0;
			b1_ = kron(obj.E2, eye(length(find(obj.Lambda.at(0) > 0))*size(obj.L, 1)))*b1;
        	eps = M1*(kron(obj.E1, eye(obj.p*obj.n_w))*(inv(M2*kron(obj.E1, eye(obj.p*obj.n_w)))*(M3 - M2*b1_) )+ b1_ ) - M4;
% 			A = kron(inv(obj.Lambda), S);
% 			z = quantity.Discrete(M1.getDomain().grid, M1.getDomain());
% 			eps = expm(-A*z)*eps.at(0);
        	Eps = reshape(eps, size(S, 1)*obj.p, obj.n);
			Eps = Eps/obj.Lambda;
		end

		function S = reconstructInternalModel(obj, S_Leader, S_init) 
			arguments
				obj;
				S_Leader quantity.Discrete;
				S_init = zeros([size(S_Leader), obj.N-1]);
			end
			ic = S_Leader(:).at(0);
            for indAgent = 2:obj.N
               S_i = S_init(:,:,indAgent-1);
               ic = cat(1, ic, S_i(:)); 
			end
			[t_S, S_pred] = ode45(@(t,S) SOde(t,S,obj,S_Leader), [S_Leader(1).domain.lower, S_Leader(1).domain.upper], ic);
			diagS_pred = kron(eye(length(find(obj.Lambda.at(0) > 0))), S_Leader);
			S.agent0 = S_Leader;
			for indAgent = 2:obj.N
               %Transform the dimensions of predicted S to what is expected
               S_i = S_pred(:, (indAgent - 1) * size(obj.L, 1)^2 + 1 : indAgent * size(obj.L, 1)^2);
               S_i = reshape(S_i, [length(t_S), size(obj.L, 1), size(obj.L, 1)]);
               S_i_disc = quantity.Discrete(interp1(t_S, S_i, S_Leader(1).domain.grid), S_Leader(1).domain);
               S.("agent" + (indAgent - 1)) = S_i_disc;
               %Build diag(S_i) for the w observer
               if isempty(diagS_pred)
                   diagS_pred = kron(eye(length(find(obj.Lambda.at(0) > 0))), S_i_disc);
               else
                   diagS_pred = blkdiag(diagS_pred, kron(eye(length(find(obj.Lambda.at(0) > 0))), S_i_disc));
			   end
			end
			S.aggregated = diagS_pred;
		end

		%% get dependent properties
        function Laplacian = get.Laplacian(obj)
           Laplacian =  diag(sum(obj.adjacencyMatrix(2:end, 2:end), 2)) - obj.adjacencyMatrix(2:end, 2:end);
        end % get.L()

        function H = get.H(obj)
            H = obj.Laplacian + diag(obj.adjacencyMatrix(2:end, 1));
        end % get.H()
        
        function N = get.N(obj)
           N = size(obj.adjacencyMatrix, 1);
        end % get.N()

        function n = get.n(obj)
           n = size(obj.Lambda, 1);
        end % get.n()
        
        function p = get.p(obj)
           p = sum(eig(obj.Lambda.at(0)) > 0);
        end % get.p()
        
        function m = get.m(obj)
           m = obj.n - obj.p;
        end % get.m()

		function n_w = get.n_w(obj)
           n_w = size(obj.L,1);
        end % get.m()
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
			nInformed = length(find(obj.adjacencyMatrix(:, 1)));
			for idx = 1:nInformed
				L = obj.L(:,:,idx);
				L_loc1 = [L_loc1; kron(eye(size(obj.L, 1)), L*a0(idx))];
				L_loc2 = blkdiag(L_loc2, kron(eye(size(obj.L, 1)), L*a0(idx)));
			end
			L_loc = [L_loc1, -L_loc2];
            L_glo = eye(obj.N - nInformed-1) * obj.mu;
			Laplacian =  diag(sum(obj.adjacencyMatrix, 2)) - obj.adjacencyMatrix;
			redLapl = Laplacian(Laplacian(:,1)==0, :);
            weightedLaplacian = L_glo * redLapl(2:end, :);
            gainMat = kron(-weightedLaplacian, eye(size(obj.L, 1)^2));
            L_loc = [L_loc, zeros(size(L_loc, 1), size(gainMat, 2)-size(L_loc, 2))];
            M = [zeros(size(obj.L, 1)^2, size(L_loc, 2)); L_loc; gainMat];
            dSdt = M*cat(1, S_leader(:).at(t), S(size(obj.L, 1)^2 + 1:end));
        end % predictS()
	end
end

