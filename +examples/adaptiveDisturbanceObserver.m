zDomain = quantity.Domain("z", linspace(0, 1, 11));
syms z t;
Lambda = quantity.Symbolic(diag([2, -1, -3]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 1 0; 1, 0 0; 0 0 -0.5], zDomain, "name", "A")*0;
Q0 = [1; 0];
Q1 = [1 0];

G1 = [0, 0; 0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 1; 0 0];
G3 = [1, 0];
G4 = [0, 0];

output = model.Output("controlOutput", "C0", [1 0 0]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "G1", G1, "G2", G2, "G3", G3, "G4", G4, "Q0", Q0, "Q1", Q1);

S = [0 1; -1 0];
P = eye(2);

distObs = AdaptiveDisturbanceObserver(mas, 2, P, 1);

[L_v, L, M] = distObs.getObserverGains(S);
timeGrid = linspace(0, 1, 21);
stateSpace = distObs.getStateSpaceApproximation(L, M, "t", timeGrid);