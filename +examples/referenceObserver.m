tDomain = quantity.Domain("t", linspace(0, 20, 4001));
adjacencyMatrix = [0 1 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
%adjacencyMatrix = [0 0; 1 0];
% S_r = [0 1; 0 0];
% w = quantity.Discrete([(tDomain.grid + 2), ones(tDomain.n, 1)], tDomain);
% S_r = [0 1; -1 0];
splitDomain = tDomain.split(10);
S_r = quantity.Piecewise({[0 1; -1 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 2; -2 0] * quantity.Discrete.ones(1, splitDomain(2))});
S_r = quantity.Discrete(S_r);
P_r = quantity.Piecewise({[0 1] * quantity.Discrete.ones(1, splitDomain(1)), [1 0] * quantity.Discrete.ones(1, splitDomain(2))});
P_r = quantity.Discrete(P_r);
% w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
w = quantity.Piecewise({quantity.Discrete([sin(splitDomain(1).grid), cos(splitDomain(1).grid)], splitDomain(1)), ...
    quantity.Discrete([sin(2*splitDomain(2).grid), cos(2*splitDomain(2).grid)], splitDomain(2))});
w = quantity.Discrete(w);

% arf = ReferenceObserver(adjacencyMatrix, 2, 4*eye(2), 4*eye(2), 4, 4, 1, 2, 2);
arf = ReferenceObserver(adjacencyMatrix, 2, [3 0; 0 2], 4, 4, 4, 1, 2, 2);

simData = arf.simulate(w, S_r, P_r);
arf.plotObservation(simData);