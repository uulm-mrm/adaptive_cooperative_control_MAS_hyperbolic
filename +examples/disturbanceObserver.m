zDomain = quantity.Domain("z", linspace(0, 1, 101));
tDomain = quantity.Domain("t", linspace(0, 10, 2001));
syms z t;
% Lambda = quantity.Symbolic(diag([2+z, -(0.5+exp(-z))]), zDomain, "name", "Lambda");
Lambda = quantity.Symbolic(diag([3, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, -1+z/4; 1.5, 0], zDomain, "name", "A");
Q0 = 0;
Q1 = 1;

output = model.Output("plant.controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 1; 1 0];
G1 = [0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
G2 = [0, 0];
G3 = [1, 0];
G4 = [0, 0];

% Lambda = quantity.Symbolic(diag([1+z/3, 0.9*cos(0.8*z), -1, -3*exp(-z)]), zDomain, "name", "Lambda");
% A = quantity.Symbolic([0, 0.8, -1+z/4, 0.2; 1/(z+1), 0, -0.4, 0; 0, 0, 0, -1; -2, z/2, 0, 0], zDomain, "name", "A");
% Q0 = [1, 0; 0 0];
% Q1 = [1 1; 1 1];
% 
% output = model.Output("controlOutput", "C0", [1 0 0 0; 0 1 0 0]);
% 
% adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
% 
% G1 = [0, 0; 0, 0; 0, 0; 0, 0]*quantity.Symbolic(1, zDomain, "name", "G1");
% G2 = [0, 0; 0 0];
% G3 = [1, 0; 0 0];
% G4 = [0, 0; 0 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

S = [0 1; -1 0];
P = eye(2);

distObserver = DisturbanceObserver(mas, S, P, "c_", 6);
% [systemTarget, controller] = mas.feedforwardControl(S, [1 0], S, P);

simulationSetting = {'t', tDomain.grid};
% ic = {"agent1.x", quantity.Symbolic(zeros(4, 1), zDomain),...
%     "agent2.x", quantity.Symbolic(zeros(4, 1), zDomain),...
%     "agent3.x", quantity.Symbolic(zeros(4, 1), zDomain),...
%     "agent4.x", quantity.Symbolic(zeros(4, 1), zDomain),...
%     "agent1.observer.x", quantity.Symbolic(zeros(4, 1), zDomain)...
%     "agent2.observer.x", quantity.Symbolic(zeros(4, 1), zDomain),...
%     "agent3.observer.x", quantity.Symbolic(zeros(4, 1), zDomain),...
%     "agent4.observer.x", quantity.Symbolic(zeros(4, 1), zDomain)};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain)};
w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
v = quantity.Discrete.zeros([2, 1], tDomain);
% dAgg = [v; w; v; w];
dAgg = [w; v];
exogenousSignals = {"disturbance", dAgg};
simData = mas.simulate(simulationSetting, ic, struct(), exogenousSignals, distObserver);

% distObserver.plotDisturbanceState(simData);
figure;
for indAgent = 1:mas.N
    for indDim = 1:length(w)
        ax = subplot(mas.N, length(v), (indAgent - 1)*length(v) + indDim);
        simDomain = simData.("agent" + indAgent).observer.disturbanceState.domain;
        plot(simDomain.grid, simData.("agent" + indAgent).observer.disturbanceState(indDim).valueDiscrete);
        hold on; 
        d = dAgg((indAgent - 1)*length(v) + 1 : indAgent * length(v));
        plot(simDomain.grid, d(indDim).valueDiscrete, '--g');
        hold off;
    end
end

% In = eye(2);
% ode = ss(kron(In, S) - kron(mas.H, distObserver.L_d*mas.E1.'*distObserver.Omega.at(1)), eye(4), eye(4), []);
% targetSystem = model.TransportOde(kron(In, Lambda), "ode", ode,...
%     "bc0", model.Output("plant.bc0", "C0", -kron(In, mas.E2.')), ...
%     "bc1", model.Output("plant.bc1", "C1", kron(In, Q1*mas.E2.' - mas.E1.'), "C", -kron(In, distObserver.F)),...
%     "pde2ode", model.Output("plant.pde2ode", "C0", -kron(In, distObserver.L_d*mas.E1.')),...
%     "ode2pde", model.Input("plant.ode2pde", "B0", zeros(2, 4)), "input", model.Input("empty", "B0", zeros([2, 1])));
% simDataTarget = targetSystem.simulate("t", tDomain.grid, targetSystem.stateName.pde(1), quantity.Symbolic.ones([4, 1], zDomain));
% simDataTarget.plant.w.plot();

pde2ode = [];
In = eye(size(mas.H, 1));
A = kron(In, S) - kron(mas.H, distObserver.L_d*mas.E1.'*distObserver.Omega.at(0));
B = kron(In, distObserver.L_d);
ev0 = [];
e0 = [];
eTarg = [];
emOrig = [];
evOrig = [];
ind0 = 1;
if mas.agent1IsLeader
   ind0 = 2; 
end
for indAgent = ind0:size(adjacencyMatrix, 1)
   ev0 = [ev0; simData.disturbance((indAgent - 1)*size(P, 1) + 1:indAgent*size(P, 1)).at(0) - simData.("agent"+indAgent).observer.disturbanceState.at(0)];
   e0 = [e0; distObserver.backsteppingKernel.transform(simData.("agent"+indAgent).x.subs("t", 0) - simData.("agent"+indAgent).observer.x.subs("t", 0))];
end
e0Tilde = e0 - kron(mas.H, distObserver.Omega)*ev0;
em0Targ = [];
em0Orig = [];
for indAgent = ind0:size(adjacencyMatrix, 1)
   systemTarget = model.Transport(Lambda, ...
       "bc0", model.Output("plant.bc0", "C0", -mas.E2.'), ...
       "bc1", model.Output("plant.bc1", "C1", Q1*mas.E2.' - mas.E1.', "C", -distObserver.F),...
       "input", model.Input("empty", "B0", zeros(size(G2, 1), 1))); 
   simDataTarget = systemTarget.simulate("t", tDomain.grid, systemTarget.stateName.pde(1), e0Tilde((indAgent - ind0)*mas.n + 1:(indAgent - ind0 + 1)*mas.n));
   pde2ode = [pde2ode; mas.E1.'*simDataTarget.plant.x.subs("z", 0)];
%    em0Targ = [em0Targ; mas.E1.'*simDataTarget.plant.x.subs("z", 0)];
%    em0Orig = [em0Orig; simData.("agent"+indAgent).error.measurement];
%    eTarg = [eTarg; simDataTarget.plant.x];
%    eOrig = [eOrig; simData.("agent" + indAgent).x - simData.("agent" + indAgent).observer.x];
%    emOrig = [emOrig; mas.E1.'*(simData.("agent" + indAgent).x - simData.("agent" + indAgent).observer.x)];
%    evOrig = [evOrig; simData.disturbance((indAgent - 1)*size(P, 1) + 1:indAgent*size(P, 1)) - simData.("agent"+indAgent).observer.disturbanceState];
end
ode = ss(A, -B, eye(length(ev0)), []);
evTarget = stateSpace.simulate(ode, [pde2ode.valueDiscrete], tDomain.grid, ev0);
% emTildeOrig = emOrig - kron(mas.H, mas.E1.'*distObserver.Omega)*evOrig;
% evOrig = stateSpace.simulate(ode, [emTildeOrig.subs("z", 0).valueDiscrete], tDomain.grid, ev0);
figure;
for ind = 1:size(evTarget, 2)
    ax = subplot(size(evTarget, 2), 1, ind);
    plot(tDomain.grid, evTarget(:, ind));
    hold on; 
    diff = dAgg(ind + (ind0 - 1)*length(v)) - simData.("agent" + (floor((ind - 1) / length(v)) + ind0)).observer.disturbanceState(mod(ind - 1, 2) + 1);
    plot(tDomain.grid, diff.valueDiscrete, '--g');
    hold off;
end
% figure;
% for ind = 1:size(ev, 2)
%     ax = subplot(size(ev, 2), 1, ind);
%     plot(tDomain.grid, evOrig(:, ind));
%     hold on; 
%     diff = dAgg(ind) - simData.("agent" + (floor((ind - 1) / length(v)) + 1)).observer.disturbanceState(mod(ind - 1, 2) + 1);
%     plot(tDomain.grid, diff.valueDiscrete, '--g');
%     hold off;
% end

% eOrigSoll = eTarg + kron(mas.H, distObserver.Omega)*quantity.Discrete(ev, tDomain);
% eOrig.plot();
% eOrigSoll.plot();

% figure;
% em0OrigDec = em0Orig - kron(mas.H, mas.E1.'*distObserver.Omega.at(0))*evOrig;
% for ind = 1:size(em0Orig, 1)
%     ax = subplot(size(em0Orig, 1), 1, ind);
%     plot(tDomain.grid, em0Targ(ind).valueDiscrete);
%     hold on; 
%     diff = em0OrigDec(ind);
%     plot(tDomain.grid, diff.valueDiscrete, '--g');
%     hold off;
% end

% simData.agent1.measurement.plot();
% hold on
% plot(tDomain.grid, simData.agent1.observer.measurement.valueDiscrete, '--g');
% hold off;

% em0TargRe = em0Targ + kron(mas.H, mas.E1.'*distObserver.Omega.at(0))*quantity.Discrete(evTarget, tDomain);
% observerMeasurementSoll1 = -(em0TargRe(1) - simData.agent1.measurement);
% observerMeasurementSoll1.plot();
% hold on;
% plot(tDomain.grid, simData.agent1.observer.measurement.valueDiscrete, '--r');
% hold off;
% 
% plot(-(em0TargRe(1) - simData.agent1.measurement) - simData.agent1.observer.measurement);
% 
% vObserverSoll = -(quantity.Discrete(evTarget(:, 1:2), tDomain) - w);
% vObserverSoll(1).plot();
% hold on;
% plot(tDomain.grid, simData.agent1.observer.disturbanceState(1).valueDiscrete, '--r');
