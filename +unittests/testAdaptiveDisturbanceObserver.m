function [tests] = testAdaptiveDisturbanceObserver()
tests = functiontests(localfunctions);
end % testMultiAgent()

function testConstructorNoFailure(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 0;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [1, 0];
G3 = [0, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

failureHappened = false;
try
    distObserver = AdaptiveDisturbanceObserver(mas, 2, 2, 2);
catch
	failureHappened = true;
end

tc.verifyFalse(failureHappened);
end % testConstructorNoFailure()

function testSolveDecouplingEquationsSimple(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 0;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [1, 0];
G3 = [0, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);
distObserver = AdaptiveDisturbanceObserver(mas, 2, 2, 2);
[Omega, Eps] = distObserver.solveDecouplingEquations(S, P);
tc.verifyEqual(Omega.median(), zeros(size(mas.Lambda,2), size(S, 1)), "AbsTol", 10^-3);
shouldBeZero = mas.Lambda*Eps.diff()-Eps*S + distObserver.G1Tilde;
tc.verifyEqual(shouldBeZero.median(), zeros(size(Omega)), "AbsTol", 10^-3);
tc.verifyEqual((mas.Q0*mas.E1.'-mas.E2.')*Eps.subs("z", 0), int(distObserver.F*Eps.subs("z", "zeta"), "zeta", 0, 1) - G2*P, "AbsTol", 10^-3);
tc.verifyEqual(mas.E1.'*Eps.subs("z", 1), G3*P, "AbsTol", 10^-3);
end % testSolveDecouplingEquationsSimple

function testSolveDecouplingEquationsAdvanced(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([2+z, -(0.5+exp(-z))]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 1+z/4; 1, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [1, 0; 1, 0]*quantity.Symbolic(1+z, zDomain, "name", "G1");
G2 = [0, 1];
G3 = [1, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);
distObserver = AdaptiveDisturbanceObserver(mas, 2, 2, 2);
[Omega, Eps] = distObserver.solveDecouplingEquations(S, P);
tc.verifyEqual(Omega.median(), zeros(size(mas.Lambda,2), size(S, 1)), "AbsTol", 10^-3);
shouldBeZero = mas.Lambda*Eps.diff()-Eps*S + distObserver.G1Tilde;
tc.verifyEqual(shouldBeZero.median(), zeros(size(Omega)), "AbsTol", 10^-3);
tc.verifyEqual((mas.Q0*mas.E1.'-mas.E2.')*Eps.subs("z", 0), int(distObserver.F*Eps.subs("z", "zeta"), "zeta", 0, 1) - G2*P, "AbsTol", 10^-3);
tc.verifyEqual(mas.E1.'*Eps.subs("z", 1), G3*P, "AbsTol", 10^-3);
end % testSolveDecouplingEquationsAdvanced

function testSolveDecouplingEquations4x4(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([2, 1, -1, -2]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0, 0, 1; 1, 0, 0, 0; 0, -1, 0, 1; 1, 0, -1, 0], zDomain, "name", "A");
Q0 = [1 0; 0 1];
Q1 = [1 0; 0 1];

output = model.Output("controlOutput", "C0", [1 1 0 0; 0 0 0 1]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
G1 = [0, 1; 0, 2; 0, 0; -1, 0]*quantity.Symbolic(1 + z, zDomain, "name", "G1");
G2 = [0 , 0; 1 0];
G3 = [1, 0; 0 1];
G4 = [0, 0; 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);
distObserver = AdaptiveDisturbanceObserver(mas, 2, 2, 2);
[Omega, Eps] = distObserver.solveDecouplingEquations(S,P);
tc.verifyEqual(Omega.median(), zeros(size(mas.Lambda,2), size(S, 1)), "AbsTol", 10^-3);
shouldBeZero = mas.Lambda*Eps.diff()-Eps*S + distObserver.G1Tilde;
tc.verifyEqual(shouldBeZero.median(), zeros(size(Omega)), "AbsTol", 10^-3);
tc.verifyEqual((mas.Q0*mas.E1.'-mas.E2.')*Eps.subs("z", 0), int(distObserver.F*Eps.subs("z", "zeta"), "zeta", 0, 1) - G2*P, "AbsTol", 10^-3);
tc.verifyEqual(mas.E1.'*Eps.subs("z", 1), G3*P, "AbsTol", 10^-3);
end % testSolveDecouplingEquations4x4

function testSolveDecouplingEquationsNonQuadraticP(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 0;

G1 = [0; 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [1];
G3 = [0];
G4 = [0];

output = model.Output("controlOutput", "C0", [1 1 0 0; 0 0 0 1]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = [1 0];
distObserver = AdaptiveDisturbanceObserver(mas, 2, 2, 2, "targetEvOde", [-4, -6]);
[Omega, Eps] = distObserver.solveDecouplingEquations(S,P);
tc.verifyEqual(Omega.median(), zeros(size(mas.Lambda,2), size(S, 1)), "AbsTol", 10^-3);
shouldBeZero = mas.Lambda*Eps.diff()-Eps*S + distObserver.G1Tilde*P;
tc.verifyEqual(shouldBeZero.median(), zeros(size(Omega)), "AbsTol", 10^-3);
tc.verifyEqual((mas.Q0*mas.E1.'-mas.E2.')*Eps.subs("z", 0), int(distObserver.F*Eps.subs("z", "zeta"), "zeta", 0, 1) - G2*P, "AbsTol", 10^-3);
tc.verifyEqual(mas.E1.'*Eps.subs("z", 1), G3*P, "AbsTol", 10^-3);
end % testSolveDecouplingEquationsNonQuadraticP

function testSolveDecouplingEquationsUnstable(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([2, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 10; 0, 0], zDomain, "name", "A");
Q0 = 0;
Q1 = 0;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [1, 0];
G3 = [0, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);
distObserver = AdaptiveDisturbanceObserver(mas, 2, 2, 2);
[Omega, Eps] = distObserver.solveDecouplingEquations(S,P);
tc.verifyEqual(Omega.median(), zeros(size(mas.Lambda,2), size(S, 1)), "AbsTol", 10^-3);
shouldBeZero = mas.Lambda*Eps.diff()-Eps*S + distObserver.G1Tilde;
tc.verifyEqual(shouldBeZero.median(), zeros(size(Omega)), "AbsTol", 10^-3);
tc.verifyEqual((mas.Q0*mas.E1.'-mas.E2.')*Eps.subs("z", 0), int(distObserver.F*Eps.subs("z", "zeta"), "zeta", 0, 1) - G2*P, "AbsTol", 10^-3);
tc.verifyEqual(mas.E1.'*Eps.subs("z", 1), G3*P, "AbsTol", 10^-3);
end % testSolveDecouplingEquationsUnstable

function testSolveKernelEquations(tc)
zDomain = quantity.Domain("z", linspace(0, 1, 201));
syms z t;
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
Q0 = 2;
Q1 = 1;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 0];
G3 = [0, 0];
G4 = [0, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

P = eye(2);
distObserver = AdaptiveDisturbanceObserver(mas, 2, 2, 2);
shouldBeZero = diff(distObserver.P_I*Lambda.subs("z", "zeta"), "zeta") + Lambda*distObserver.P_I.diff("z") + A*distObserver.P_I;
tc.verifyEqual(shouldBeZero.median(), zeros(size(distObserver.P_I)), "AbsTol", 10^-2);
tc.verifyEqual(median((mas.Q0*mas.E1.'-mas.E2.')*distObserver.P_I.subs("z", 0)), median(distObserver.F), "AbsTol", 10^-3);
tc.verifyEqual(median(distObserver.P_I.subs("zeta", "z")*Lambda - Lambda*distObserver.P_I.subs("zeta", "z")), median(A), "AbsTol", 10^-3);
end % testSolveKernelEquations

% function testSolveKernelEquationsAdvanced(tc)
% zDomain = quantity.Domain("z", linspace(0, 1, 101));
% syms z t;
% Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
% A = quantity.Symbolic([0, 1; 0, 0], zDomain, "name", "A");
% Q0 = 2;
% Q1 = 1;
% 
% output = model.Output("plant.controlOutput", "C0", [1 0]);
% 
% adjacencyMatrix = [0 1; 1 0];
% G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
% G2 = [0, 0];
% G3 = [0, 0];
% G4 = [0, 0];
% 
% mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
%     "G1", G1, "G2", G2, "G3", G3, "G4", G4);
% 
% P = eye(2);
% distObserver = AdaptiveDisturbanceObserver(mas, 2, P, 2);
% shouldBeZero = diff(distObserver.P_I*Lambda.subs("z", "zeta"), "zeta") + Lambda*distObserver.P_I.diff("z") + A*distObserver.P_I;
% tc.verifyEqual(shouldBeZero.median(), zeros(size(distObserver.P_I)), "AbsTol", 10^-2);
% tc.verifyEqual(median((mas.Q0*mas.E1.'-mas.E2.')*distObserver.P_I.subs("z", 0)), median(distObserver.F), "AbsTol", 10^-3);
% tc.verifyEqual(median(distObserver.P_I.subs("zeta", "z")*Lambda - Lambda*distObserver.P_I.subs("zeta", "z")), median(A), "AbsTol", 10^-3);
% end % testSolveKernelEquationsAdvanced
% 
% function testSolveKernelEquations4x4(tc)
% zDomain = quantity.Domain("z", linspace(0, 1, 101));
% syms z t;
% Lambda = quantity.Symbolic(diag([2, -1, -2, -3]), zDomain, "name", "Lambda");
% A = quantity.Symbolic([0, 0, 0, 1; 0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0], zDomain, "name", "A");
% Q0 = [2; 0; 1];
% Q1 = [1 0 0];
% 
% output = model.Output("plant.controlOutput", "C0", [1 0]);
% 
% adjacencyMatrix = [0 1; 1 0];
% G1 = [0, 0; 0, 0; 0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
% G2 = [0, 0; 0, 0; 0, 0];
% G3 = [0, 0];
% G4 = [0, 0; 0, 0];
% 
% mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
%     "G1", G1, "G2", G2, "G3", G3, "G4", G4);
% 
% P = eye(2);
% distObserver = AdaptiveDisturbanceObserver(mas, 2, P, 2);
% shouldBeZero = diff(distObserver.P_I*Lambda.subs("z", "zeta"), "zeta") + Lambda*distObserver.P_I.diff("z") + A*distObserver.P_I;
% tc.verifyEqual(shouldBeZero.median(), zeros(size(distObserver.P_I)), "AbsTol", 10^-2);
% tc.verifyEqual(median((mas.Q0*mas.E1.'-mas.E2.')*distObserver.P_I.subs("z", 0)), median(distObserver.F), "AbsTol", 10^-3);
% tc.verifyEqual(median(distObserver.P_I.subs("zeta", "z")*Lambda - Lambda*distObserver.P_I.subs("zeta", "z")), median(A), "AbsTol", 10^-3);
% end % testSolveKernelEquations4x4