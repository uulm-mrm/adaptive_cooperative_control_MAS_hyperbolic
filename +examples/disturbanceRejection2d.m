zDomain = quantity.Domain("z", linspace(0, 1, 201));
tDomain = quantity.Domain("t", linspace(0, 8, 1601));
syms z t;
Lambda = quantity.Symbolic(diag([2, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 1; 0, 0], zDomain, "name", "A");
Q0 = 1;
Q1 = 1;

output = model.Output("controlOutput", "C0", [1 0]);

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 1 0 0];
% adjacencyMatrix = [0 0; 1 0];
G1 = [0, 1; 0, 0]*quantity.Symbolic(1 + z, zDomain, "name", "G1");
G2 = [0, 1];
G3 = [0, 0];
G4 = [1, 0];

mas = MultiAgent(adjacencyMatrix, Lambda, "A", A, "output", output, "Q0", Q0, "Q1", Q1,...
    "G1", G1, "G2", G2, "G3", G3, "G4", G4);

simulationSetting = {'t', tDomain.grid};
ic = {"agent1.x", quantity.Symbolic(zeros(2, 1), zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent3.x", quantity.Symbolic([0; 0], zDomain),...
    "agent4.x", quantity.Symbolic([0; 0], zDomain)};
S_r = [0 1; -1 0];
S = [0 1; -1 0];
P = eye(2);
[systemTarget, controller] = mas.feedforwardControl(S_r, [1 0], S, P);

%% 

w = quantity.Discrete([sin(tDomain.grid), cos(tDomain.grid)], tDomain);
v = quantity.Discrete.zeros([2, 1], tDomain);
% exogenousSignals = {'disturbanceState', [0; 0; 2*w; 0; 0; w], "disturbance", [0; 0; 2*w; 0; 0; w]};
aggregatedDisturbance = [v; v; v; w];
zeroDisturbance = [v; v; v; v];
exogenousSignals = {'disturbanceState', aggregatedDisturbance, "disturbance", aggregatedDisturbance};

simDataControlBoth = mas.simulate(simulationSetting, ic, controller, exogenousSignals);

for indAgent = 1:mas.N
    simDataTarget = systemTarget.simulate("t", tDomain.grid, ...
        systemTarget.stateName.pde, controller.backsteppingKernel.backsteppingTransformation(ic{2*indAgent}));
    simDataControlBoth.("agent"+indAgent).reference = simDataTarget.plant.controlOutput;
end

% mas.plotState(simDataControl);


% exogenousSignals = {'disturbanceState', zeroDisturbance, "disturbance", aggregatedDisturbance};
% 
% simDataControlDisturbance = mas.simulate(simulationSetting, ic, controller, exogenousSignals);
% 
% for indAgent = 1:mas.N
%     simDataTarget = systemTarget.simulate("t", tDomain.grid, ...
%         systemTarget.stateName.pde, controller.backsteppingKernel.backsteppingTransformation(ic{2*indAgent}));
%     simDataControlDisturbance.("agent"+indAgent).reference = simDataTarget.plant.controlOutput;
% end
% 
% 
% 
% exogenousSignals = {'disturbanceState', aggregatedDisturbance, "disturbance", zeroDisturbance};
% 
% simDataControlFeedforward = mas.simulate(simulationSetting, ic, controller, exogenousSignals);
% 
% for indAgent = 1:mas.N
%     simDataTarget = systemTarget.simulate("t", tDomain.grid, ...
%         systemTarget.stateName.pde, controller.backsteppingKernel.backsteppingTransformation(ic{2*indAgent}));
%     simDataControlFeedforward.("agent"+indAgent).reference = simDataTarget.plant.controlOutput;
% end
% 
% mas.plotOutput(simDataControlDisturbance);
% mas.plotOutput(simDataControlFeedforward);
mas.plotOutput(simDataControlBoth);

%Check if PI_v solves regulator equations
% shouldBeZero = mas.Lambda*diff(controller.PI_v)+controller.backsteppingKernel.getA0Target()*mas.E1.'*controller.PI_v.at(0)...
%     -controller.PI_v*S+controller.backsteppingKernel.backsteppingTransformation(mas.G1)*P;
% plot(shouldBeZero);
% mas.E2.'*controller.PI_v.at(0)-mas.Q0*mas.E1.'*controller.PI_v.at(0) - mas.G2*P
% outputTransformed = mas.output.backstepping(controller.backsteppingKernel, "inverse", true);
% outputTransformed.out(controller.PI_v) + mas.G4*P
% 
% K=controller.backsteppingKernel.getValue();
% shouldBeZero = mas.Lambda*K.diff("z")+diff(K*subs(mas.Lambda, "z", "zeta"), "zeta")-K*subs(mas.A, "z", "zeta");
% median(shouldBeZero)
%% 

% stableSystem = controller.backsteppingKernel.getTargetSystem(mas.network.agent(4));
% simDataStableRejection = stableSystem.simulate("t", tDomain.grid, "agent4.control", controller.disturbanceGain.agent4*aggregatedDisturbance, "disturbance", aggregatedDisturbance);
% % stableSystem.plotNorm(simDataStableRejection.agent4.controlOutput);
% % 
% distTransform = kron(mas.H, controller.PI_v);
% icTranformed = ic{8} - distTransform(end-1:end, :)*aggregatedDisturbance;
% simDataSystemTarget = systemTarget.simulate("t", tDomain.grid, systemTarget.stateName.pde, icTranformed);
% % systemTarget.plotNorm(simDataSystemTarget.plant.controlOutput);
% plot(subs(simDataSystemTarget.plant.x, "z", 0));
% xTilde = simDataSystemTarget.plant.x + distTransform(end-1:end, :)*aggregatedDisturbance;
% plot(subs(xTilde, "z", 0));
% xOrig = controller.backsteppingKernel.inverseBacksteppingTransformation(xTilde);
% plot(subs(xOrig, "z", 0));
% 
% % simDataAgent4 = mas.network.agent(4).simulate("t", tDomain.grid, "agent4.control", controller.disturbanceGain.agent4*aggregatedDisturbance, "disturbance", aggregatedDisturbance);
% % plot(simDataAgent4.agent4.x);
% % a4Tilde = controller.backsteppingKernel.backsteppingTransformation(simDataAgent4.agent4.x);
% % plot(a4Tilde);
% % a4Target = a4Tilde - distTransform(end-1:end, :)*aggregatedDisturbance;
% % plot(a4Target)
% %% 
% 
% A0 = controller.backsteppingKernel.getA0Target();
% G1dash = controller.backsteppingKernel.backsteppingTransformation(G1) * P;
% PI_v = tool.solveGenericBvpFirstOrderCheb(S, Lambda, quantity.Discrete.zeros(size(Lambda), mas.domain),...
%                 -G1dash, A0, mas.Q0, G2 * P, output, -G4 * P); 
% K = -1 * kron(mas.H, G3 * P - mas.E1.'*PI_v.at(1));
% G2tilde = kron(mas.H, G2 * P);
% G1tilde = kron(mas.H, G1dash);
% G4tilde = kron(mas.H, G4 * P);
% % systemTest = model.Transport(kron(eye(4), Lambda), "A0", kron(eye(4), A0), "Q0", kron(eye(4), mas.Q0), ...
% %                 "output", kron(eye(4), output), ...
% %                 "input", model.Input("disturbance", "B", kron(mas.H, G1dash), ...
% %                 "B0", kron(mas.H, G2 * P), "B1", kron(mas.H, G3 * P), "D", misc.Gain("disturbance", kron(mas.H, G4 * P), "outputType", "controlOutput")) + ...
% %                 model.Input("disturbanceState", "B1", K));
% % systemTest = model.Transport(Lambda, "A0", A0, "Q0", Q0, ...
% %                 "output", output, ...
% %                 "input", model.Input("disturbanceState", "B1", K(end, :)) + model.Input("disturbance", "B0", G2tilde(end, :)));
% systemTest = model.Transport(Lambda, "A0", A0, "Q0", Q0, ...
%                 "output", output, ...
%                 "input", model.Input("disturbanceState", "B1", K(end, :)) + model.Input("disturbance", "B0", G2tilde(end, :), "B", G1tilde(end-1:end, :)), ...
%                 "D", misc.Gain("disturbance", G4tilde(end, :)));
%            
% % zeroBackstepped = controller.backsteppingKernel.backsteppingTransformation(quantity.Symbolic.zeros([2, 1], zDomain));
% % icTest = kron(ones(4, 1), zeroBackstepped);
% 
% % simDataNoDistNoCont = systemTest.simulate("t", tDomain.grid, systemTest.stateName.pde, icTest);
% simDataDistNoCont = systemTest.simulate("t", tDomain.grid, systemTest.stateName.pde, quantity.Symbolic([0; 0], zDomain), "disturbance", aggregatedDisturbance);
% simDataDistCont = systemTest.simulate("t", tDomain.grid, systemTest.stateName.pde, quantity.Symbolic([0; 0], zDomain), "disturbance", aggregatedDisturbance, "disturbanceState", aggregatedDisturbance);
% simDataNoDistCont = systemTest.simulate("t", tDomain.grid, systemTest.stateName.pde, quantity.Symbolic([0; 0], zDomain), "disturbanceState", aggregatedDisturbance);
% % plot(simDataNoDistNoCont.plant.x);
% % plot(simDataDistNoCont.plant.x);
% % plot(simDataNoDistCont.plant.x);
% % plot(simDataDistCont.plant.x);
% % plot(subs(simDataDistCont.plant.x(4), "z", 0));
% % plot(subs(simDataDistCont.plant.x(8), "z", 0));
% plot(subs(simDataControlBoth.agent4.x, "z", 0));
% plot(subs(xTilde, "z", 0));
% plot(subs(simDataDistCont.plant.x, "z", 0));
% plot(subs(simDataStableRejection.agent4.x, "z", 0));