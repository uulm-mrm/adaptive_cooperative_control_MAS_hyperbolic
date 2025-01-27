zDomain = quantity.Domain("z", linspace(0, 1, 61));
syms z t;
doexport = true;

adjacencyMatrix = [0 0 0 0; 2 0 0 0; 0 1 0 0.5; 0 1 0.5 0];
Lambda = quantity.Symbolic(diag([1, -1]), zDomain, "name", "Lambda");
A = quantity.Symbolic([0, 0; 0, 0], zDomain, "name", "A");
Q0 = -1;
Q1 = 1;

g1 = quantity.Discrete(gaussmf(zDomain.grid, [0.1 0.25]), zDomain, "name", "G1")/1;
G1 = [1; 1]*g1;
G2 = [0];
G3 = [0];
G4 = [0];

controlInput = model.Input("control", "B1", eye(1));
disturbanceInput = model.Input("disturbance", "B", G1,...
	"B0", G2, "B1", G3, "D", misc.Gain("disturbance", G4, "outputType", "controlOutput"));
input = controlInput + disturbanceInput;
% output = model.Output("controlOutput", "C0", [1 0]/2);
% output = model.Output("controlOutput", "zk", 0.6, "Ck", [1 1]/2);
C = quantity.Discrete([ones(floor(length(zDomain.grid)/2), 1); zeros(floor(1*length(zDomain.grid)/2)+1, 1)]/2, zDomain);
% C = quantity.Discrete.ones(1, zDomain);
% C = quantity.Discrete([ones(floor(length(zDomain.grid)/2), 1); -ones(floor(length(zDomain.grid)/2)+1, 1)]/2, zDomain);

lambdas = [1, 1.2, 1.5, .9];
for idx = 1:4
	a = lambdas(idx);
	Lambdaa = quantity.Symbolic(diag([a, -a]), zDomain, "name", "Lambda"); 
	output = model.Output("controlOutput", "C", [C, -C]/a);
	agents{idx} = model.Transport(Lambdaa, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(idx));
end

simTime = 80;
tDomain = quantity.Domain("t", linspace(0, simTime, 16001));
simulationSetting = {'t', tDomain.grid};


%% Simulation Führungsverhalten + Störverhalten
splitDomain = tDomain.split(simTime/2);
Pr = [1 0]/2;
% freq = 1*pi;
% freq2 = 1.5*pi;
% freqd = 0.6*pi;
% freqd2 = 1*pi;

freq = 0.4*pi;
freq2 = 0.6*pi;
freqd = 0.6*pi;
freqd2 = 0.8*pi;

S = quantity.Piecewise({blkdiag([0 freq; -freq 0], [0 freqd; -freqd 0]) * quantity.Discrete.ones(1, splitDomain(1)), blkdiag([0 freq2; -freq2 0], [0 freqd2; -freqd2 0]) * quantity.Discrete.ones(1, splitDomain(2))});
S = quantity.Discrete(S);
% Sr = quantity.Piecewise({[0 freq; -freq 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 freq2; -freq2 0] * quantity.Discrete.ones(1, splitDomain(2))});
% Sr = quantity.Discrete(Sr);
w = quantity.Piecewise({quantity.Discrete([sin(freq*splitDomain(1).grid), cos(freq*splitDomain(1).grid)], splitDomain(1)), ...
    quantity.Discrete([sin(freq2*splitDomain(2).grid), cos(freq2*splitDomain(2).grid)], splitDomain(2))});
w = quantity.Discrete(w);
P1 = quantity.Piecewise({[6 0] * quantity.Discrete.ones(1, splitDomain(1)), [4 0] * quantity.Discrete.ones(1, splitDomain(2))});
P2 = quantity.Piecewise({[3 0] * quantity.Discrete.ones(1, splitDomain(1)), [5 0] * quantity.Discrete.ones(1, splitDomain(2))});
P3 = quantity.Piecewise({[5 0] * quantity.Discrete.ones(1, splitDomain(1)), [2 0] * quantity.Discrete.ones(1, splitDomain(2))});
P4 = quantity.Piecewise({[5 0] * quantity.Discrete.ones(1, splitDomain(1)), [6 0] * quantity.Discrete.ones(1, splitDomain(2))});
% P = quantity.Piecewise({[1 0] * quantity.Discrete.ones(1, splitDomain(1)), [1 0] * quantity.Discrete.ones(1, splitDomain(2))});
P1 = quantity.Discrete(P1);
P2 = quantity.Discrete(P2);
P3 = quantity.Discrete(P3);
P4 = quantity.Discrete(P4);
% S = quantity.Piecewise({[0 freq; -freq 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 freq2; -freq2 0] * quantity.Discrete.ones(1, splitDomain(2))});
% S = quantity.Discrete(S);
v = quantity.Piecewise({quantity.Discrete([sin(freqd*splitDomain(1).grid), cos(freqd*splitDomain(1).grid)], splitDomain(1)), ...
    quantity.Discrete([sin(freqd2*splitDomain(2).grid), cos(freqd2*splitDomain(2).grid)], splitDomain(2))});
v = quantity.Discrete(v);

b = {ones(length(S), 1)*4};
mu = 4;
L = 4*eye(length(S));
imc = InternalModelController(adjacencyMatrix, Lambda, b, mu, L, "A", A, "output", output, "Q0", Q0, "Q1", Q1, "nCoef", 20);
mas = MultiAgentHeterogenous(agents, adjacencyMatrix, "nCoefDist", 5, "nCoefRef", 5);


aggregatedDisturbance = [P1*v; P2*v; P3*v; P4*v];

ic = {"agent1.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent3.x", quantity.Symbolic([0; 0], zDomain),...
    "agent4.x", quantity.Symbolic([0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent3.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent4.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "S", [0 1; -1 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0],...
    "v", [0; 0; 0; 0; 0; 0; 0; 0]};
r = Pr*w;
exogenousInput = {"disturbance", aggregatedDisturbance, "S", S, "reference", r};
% simData = mas.adaptiveControlSimulation(simulationSetting, 2, 1, arf, "initialCondition", ic, "exogenousInput", exogenousInput, "mu_S", 2, ...
% 	"distCoef", 10, "adaptationDelayObserver", 10, "adaptationDelayController", 10);
c = 4;
R = 4;
Q = eye(length(S)*output.lengthOutput)*4;
adaptationDelay = 100;
simData = tool.simulateInternalModel(simulationSetting, mas, imc, "initialCondition", ic, "exogenousInput", exogenousInput, "c", c, "R", R, "Q", Q, "adaptationDelay", adaptationDelay);

for idx = 1:4
	simData.("agent" + idx).controlOutput.plot();
end

%% Export Führungsverhalten + Störverhalten
if doexport
tPlot = quantity.Domain("t", linspace(0, simTime, 1601));
zPlot = quantity.Domain("z", linspace(0, 1, 101));
basepath = "C:\Users\xzb84\Documents\tarikspapers\VortragGMA_2_15_2024\data";
% basepath = "D:\Eigene Dateien\Promotion\Vorträge\mrmbeamer-main\data";
% r = Pr*w;
for i = 1:mas.N
    omega = 2*(5-i)-1.5+0.5*[1 1]*simData.("agent"+i).x.subs("t","tau").int("tau", 0, "t");
    omega = omega.changeDomain([tPlot, zPlot]);
    
    dataSim = zPlot.grid;
    signalName = ["z"];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "omega" + string(tIt)];
        dataSim = cat(2, dataSim, omega.subs("t", tPlot.grid(tIt)).valueDiscrete);
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "simData"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);
end

for i = 1:mas.N
	dataSim = [];
	signalName = [];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "y" + string(tIt)];
        dataSim = cat(2, dataSim, 2*(5-i)-1.5+simData.("agent"+i).controlOutput.subs("t", tPlot.grid(tIt)));
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "outputData"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);
end

for i = 1:mas.N
	dataSim = [];
	dataSimNorm = [];
	signalName = [];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "y" + string(tIt)];
		newData = (simData.("agent"+i).controlOutput.subs("t", tPlot.grid(tIt))-r.subs("t", tPlot.grid(tIt)))...
			/max(abs(simData.("agent"+i).controlOutput - r));
		if abs(newData) > 0.1
			dataSim = cat(2, dataSim, newData);
			dataSimNorm = cat(2, dataSimNorm, nthroot(norm(newData), 2)/2);
		else
			dataSim = cat(2, dataSim, 0);
			dataSimNorm = cat(2, dataSimNorm, 0);
		end
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "normOutputErrorData"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);exportObj = export.dd(...
        "M", dataSimNorm, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "outputErrorData"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);
end

for i = 1:mas.N
	u = simData.("agent"+i).control/(2*lambdas(i))/2;
	dataSim = [];
	signalName = [];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "u" + string(tIt)];
        dataSim = cat(2, dataSim, u.subs("t", tPlot.grid(tIt)));
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "inputData"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);
end

% for i = 1:mas.N
% 	dataSim = [];
% 	signalName = [];
% 	for tIt = 1 : numel(tPlot.grid)
%         signalName = [signalName, "y" + string(tIt)];
% 		newData = nthroot(real(norm((simData.("agent"+i).controlOutput.subs("t", tPlot.grid(tIt))-r.subs("t", tPlot.grid(tIt))))),4)...
% 			/max(1,nthroot(real(max(norm(simData.("agent"+i).controlOutput - r))),4))*0.5;
% 		if abs(newData) > 0.22
% 			dataSim = cat(2, dataSim, newData);
% 		else
% 			dataSim = cat(2, dataSim, 0);
% 		end
%     end % for tIt = tSample.'
% %     for tIt = 1 : numel(tPlot.grid)
% %         signalName = [signalName, "y" + string(tIt)];
% % 		if tIt < 2
% % 			dataSim = cat(2, dataSim, 0);
% % 		elseif tIt < numel(tPlot.grid)/2
% % 			dataSim = cat(2, dataSim, nthroot(real(norm(simData.("agent"+i).controlOutput.subs("t", tPlot.grid(tIt))-r.subs("t", tPlot.grid(tIt)))),4)...
% % 				/max(1,nthroot(real(max(norm(simData.("agent"+i).controlOutput - r))),4))*0.6);
% % 		elseif tIt < numel(tPlot.grid)*3/4-5
% % 			dataSim = cat(2, dataSim, nthroot(norm(simData.("agent"+i).controlOutput.subs("t", tPlot.grid(tIt))-r.subs("t", tPlot.grid(tIt))),4)...
% % 				/max(1,nthroot(norm(max(simData.("agent"+i).controlOutput - r)),4))*0.6);
% % 		else
% % 			dataSim = cat(2, dataSim, 0);
% % 		end
% % 		elseif tIt < numel(tPlot.grid)/2
% % 			dataSim = cat(2, dataSim, nthroot(norm(simData.("agent"+i).S.subs("t", tPlot.grid(tIt)) - S.subs("t", tPlot.grid(tIt))),4)...
% % 				/max(1,nthroot(norm(simData.("agent"+i).S.at(0)-S.at(0)),4))*0.6);
% % 		elseif tIt < numel(tPlot.grid)*3/4-5
% % 			dataSim = cat(2, dataSim, nthroot(norm(simData.("agent"+i).S.subs("t", tPlot.grid(tIt)) - S.subs("t", tPlot.grid(tIt))),4)...
% % 				/max(1,nthroot(norm(simData.("agent"+i).S.at(10.01)-S.at(10.01)),4))*0.6);
% % 		else
% % 			dataSim = cat(2, dataSim, 0);
% % 		end
% %     end % for tIt = tSample.'
% 
% %     exportObj = export.dd(...
% %         "M", dataSim, ...
% %         "header", table2cell(array2table(signalName)).', ...
% %         "filename", "outputErrorData"+i, ...
% %         "basepath", basepath);
% %     exportObj.export("silent", true);
% end

dataSim = tPlot.grid;
signalName = "t";
for i = 1:mas.N
	dataSim = cat(2, dataSim, simData.("agent" + i).controlOutput.on(tPlot));
	signalName = cat(2, signalName, "y" + i);
end
exportObj = export.dd(...
    "M", dataSim, ...
    "header", signalName, ...
    "filename", "output", ...
    "basepath", basepath);
exportObj.export("silent", true);

for i = 1:mas.N
	d = aggregatedDisturbance((i-1)*length(P1*v)+1:i*length(P1*v));
	d_eff = g1*d/3;
    d_eff = d_eff.changeDomain([tPlot, zPlot]);
    
    dataSim = zPlot.grid;
    signalName = ["z"];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "d" + string(tIt)];
        dataSim = cat(2, dataSim, d_eff.subs("t", tPlot.grid(tIt)).valueDiscrete);
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "disturbance"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);
end

data.r = export.dd(...
    'M', [tDomain.grid, r.valueDiscrete], ...
    'header', {'t', 'r1'}, ...
    'filename', 'leader', ...
    'basepath', basepath, ...
	'N', 400 ...
    );

dataSim = [];
signalName = [];
for tIt = 1 : numel(tPlot.grid)
    signalName = [signalName, "r" + string(tIt)];
    dataSim = cat(2, dataSim, r.subs("t", tPlot.grid(tIt)));
end % for tIt = tSample.'
data.lvec = export.dd(...
    "M", dataSim, ...
    "header", table2cell(array2table(signalName)).', ...
    "filename", "leaderVec", ...
    "basepath", basepath);


dataSim = [];
signalName = [];
for tIt = 1 : numel(tPlot.grid)
    signalName = [signalName, "t" + string(tIt)];
    dataSim = cat(2, dataSim, floor(tPlot.grid(tIt)));
end % for tIt = tSample.'

data.t = export.dd(...
    "M", dataSim, ...
    "header", table2cell(array2table(signalName)).', ...
    "filename", "timeVec", ...
    "basepath", basepath);

        
export.Data.exportAll(data);
end