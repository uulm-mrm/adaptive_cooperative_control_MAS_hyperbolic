zDomain = quantity.Domain("z", linspace(0, 1, 61));
syms z t;

adjacencyMatrix = [0 0 0 0; 1 0 0 0; 0 1 0 1; 0 1 1 0];
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
C = quantity.Discrete([ones(floor(length(zDomain.grid)/2), 1); zeros(floor(1*length(zDomain.grid)/2)+1, 1)]/2, zDomain);
lambdas = [1, 1.5, 2, 0.9];
for idx = 1:4
	a = lambdas(idx);
	Lambda = quantity.Symbolic(diag([a, -a]), zDomain, "name", "Lambda"); 
	output = model.Output("controlOutput", "C", [C, -C]/a);
	agents{idx} = model.Transport(Lambda, "A", A, "input", input, "output", output, "Q0", Q0, "Q1", Q1, "prefix", "agent"+(idx));
end

simTime = 20;
tDomain = quantity.Domain("t", linspace(0, simTime, 2001));
simulationSetting = {'t', tDomain.grid};

%% Simulation Führungsverhalten
splitDomain = tDomain.split(simTime/2);
P = [1 0];
S = [0 1; -1 0];
Pr = [1 0]/2;
freq = 1*pi;
freq2 = 1.5*pi;
Sr = quantity.Piecewise({[0 freq; -freq 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 freq2; -freq2 0] * quantity.Discrete.ones(1, splitDomain(2))});
Sr = quantity.Discrete(Sr);
w = quantity.Piecewise({quantity.Discrete([sin(freq*splitDomain(1).grid), cos(freq*splitDomain(1).grid)], splitDomain(1)), ...
    quantity.Discrete([sin(freq2*splitDomain(2).grid), cos(freq2*splitDomain(2).grid)], splitDomain(2))});
w = quantity.Discrete(w);
v = quantity.Discrete([-sin(tDomain.grid), -cos(tDomain.grid)], tDomain)*0;

mas = MultiAgentHeterogenous(agents, adjacencyMatrix, "nCoefDist", 5, "nCoefRef", floor(max(freq, freq2)*5));


aggregatedDisturbance = [P*v*0; P*v; P*v*0; P*v/2];

obsGain = 2;
arf = ReferenceObserver(adjacencyMatrix, 2, obsGain, obsGain, obsGain, obsGain, 1, obsGain, obsGain);
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
exogenousInput = {"disturbance", aggregatedDisturbance, "S", S, "Sr", Sr, "w", w, "Pr", Pr, "P", P};
simData = mas.adaptiveControlSimulation(simulationSetting, 2, 1, arf, "initialCondition", ic, "exogenousInput", exogenousInput, "mu_S", 2, ...
	"distCoef", 10, "adaptationDelayObserver", 10, "adaptationDelayController", 10);

for idx = 1:4
	simData.("agent" + idx).controlOutput.plot();
end
%% Export Führungsverhalten
tPlot = quantity.Domain("t", linspace(0, simTime, 401));
zPlot = quantity.Domain("z", linspace(0, 1, 101));
basepath = "C:\Users\xzb84\Documents\Vorträge\GMA2022\mrmbeamer-main\data";
r = Pr*w;
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
	signalName = [];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "y" + string(tIt)];
		newData = (simData.("agent"+i).controlOutput.subs("t", tPlot.grid(tIt))-r.subs("t", tPlot.grid(tIt)))...
			/max(abs(simData.("agent"+i).controlOutput - r));
		if abs(newData) > 0.1
			dataSim = cat(2, dataSim, newData);
		else
			dataSim = cat(2, dataSim, 0);
		end
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "normOutputErrorData"+i, ...
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

for i = 1:mas.N
	dataSim = [];
	signalName = [];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "y" + string(tIt)];
		if tIt < 2
			dataSim = cat(2, dataSim, 0);
		elseif tIt < numel(tPlot.grid)/2
			dataSim = cat(2, dataSim, nthroot(norm(simData.("agent"+i).observer.Sr.subs("t", tPlot.grid(tIt)) - Sr.subs("t", tPlot.grid(tIt))),4)...
				/max(1,nthroot(norm(simData.("agent"+i).observer.Sr.at(0)-Sr.at(0)),4))*0.6);
		elseif tIt < numel(tPlot.grid)*3/4-5
			dataSim = cat(2, dataSim, nthroot(norm(simData.("agent"+i).observer.Sr.subs("t", tPlot.grid(tIt)) - Sr.subs("t", tPlot.grid(tIt))),4)...
				/max(1,nthroot(norm(simData.("agent"+i).observer.Sr.at(10.01)-Sr.at(10.01)),4))*0.6);
		else
			dataSim = cat(2, dataSim, 0);
		end
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "outputErrorData"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);
end


dataSim = [];
signalName = [];
for tIt = 1 : numel(tPlot.grid)
    signalName = [signalName, "t" + string(tIt)];
    dataSim = cat(2, dataSim, tPlot.grid(tIt));
end % for tIt = tSample.'

data.t = export.dd(...
    "M", dataSim, ...
    "header", table2cell(array2table(signalName)).', ...
    "filename", "timeVec", ...
    "basepath", basepath);

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
data.t = export.dd(...
    "M", dataSim, ...
    "header", table2cell(array2table(signalName)).', ...
    "filename", "leaderVec", ...
    "basepath", basepath);

dataSim = tPlot.grid;
signalName = "t";
for i = 1:mas.N
	dataSim = cat(2, dataSim, simData.("agent" + i).observer.Sr.norm().on(tPlot));
	signalName = cat(2, signalName, "S" + i);
end
exportObj = export.dd(...
    "M", dataSim, ...
    "header", table2cell(array2table(signalName)).', ...
    "filename", "normFObs", ...
    "basepath", basepath);
exportObj.export("silent", true);

dataSim = tPlot.grid;
signalName = "t";
for i = 1:mas.N
	dataSim = cat(2, dataSim, simData.("agent" + i).observer.Pr.norm().on(tPlot));
	signalName = cat(2, signalName, "P" + i);
end
exportObj = export.dd(...
    "M", dataSim, ...
    "header", table2cell(array2table(signalName)).', ...
    "filename", "normRObs", ...
    "basepath", basepath);
exportObj.export("silent", true);

dataSim = tPlot.grid;
signalName = "t";
for i = 1:mas.N
	dataSim = cat(2, dataSim, simData.("agent" + i).observer.w.on(tPlot));
	signalName = cat(2, signalName, "w1" + i);
	dataSim = cat(2, dataSim, simData.("agent" + i).observer.w.on(tPlot));
	signalName = cat(2, signalName, "w2" + i);
end
exportObj = export.dd(...
    "M", dataSim, ...
    "header", signalName, ...
    "filename", "wObs", ...
    "basepath", basepath);
exportObj.export("silent", true);

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

        
export.Data.exportAll(data);

%% Simulation Störverhalten
splitDomain = tDomain.split(simTime/2);
Pr = [1 0];
Sr = [0 1; -1 0];
w = quantity.Discrete.zeros([2, 1], tDomain);
P = quantity.Piecewise({[6 0] * quantity.Discrete.ones(1, splitDomain(1)), [4 0] * quantity.Discrete.ones(1, splitDomain(2))});
P = quantity.Discrete(P);
freq = 0.6*pi;
freq2 = 1*pi;
S = quantity.Piecewise({[0 freq; -freq 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 freq2; -freq2 0] * quantity.Discrete.ones(1, splitDomain(2))});
S = quantity.Discrete(S);
v = quantity.Piecewise({quantity.Discrete([sin(freq*splitDomain(1).grid), cos(freq*splitDomain(1).grid)], splitDomain(1)), ...
    quantity.Discrete([sin(freq2*splitDomain(2).grid), cos(freq2*splitDomain(2).grid)], splitDomain(2))});
v = quantity.Discrete(v);

mas = MultiAgentHeterogenous(agents, adjacencyMatrix, "nCoefDist", floor(max(freq, freq2)*5), "nCoefRef", 2);


aggregatedDisturbance = [P*v; P*v; P*v; P*v];

obsGain = 2;
arf = ReferenceObserver(adjacencyMatrix, 2, obsGain, obsGain, obsGain, obsGain, 1, obsGain, obsGain);
ic = {"agent1.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.x", quantity.Symbolic([0; 0], zDomain),...
    "agent3.x", quantity.Symbolic([0; 0], zDomain),...
    "agent4.x", quantity.Symbolic([0; 0], zDomain),...
    "agent1.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent2.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent3.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "agent4.observer.x", quantity.Symbolic([0; 0], zDomain),...
    "S", [0 freq; -freq 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0],...
    "v", [0; 0; 0; 0; 0; 0; 0; 0]};
exogenousInput = {"disturbance", aggregatedDisturbance, "S", S, "Sr", Sr, "w", w, "Pr", Pr, "P", P};
simData = mas.adaptiveControlSimulation(simulationSetting, 2, 1, arf, "initialCondition", ic, "exogenousInput", exogenousInput, "mu_S", 2, ...
	"distCoef", 11, "adaptationDelayObserver", 10, "adaptationDelayController", 10);

for idx = 1:4
	simData.("agent" + idx).controlOutput.plot();
end
%% Export Störverhalten
tPlot = quantity.Domain("t", linspace(0, simTime, 401));
zPlot = quantity.Domain("z", linspace(0, 1, 101));
basepath = "C:\Users\xzb84\Documents\Vorträge\GMA2022\mrmbeamer-main\dataDist";
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
	u = simData.("agent"+i).control/(2*lambdas(i));
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

for i = 1:mas.N
	dataSim = [];
	signalName = [];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "y" + string(tIt)];
		if tIt < 2
			dataSim = cat(2, dataSim, 0);
		elseif tIt < numel(tPlot.grid)/2
			dataSim = cat(2, dataSim, nthroot(norm(simData.("agent"+i).observer.S.subs("t", tPlot.grid(tIt)) - S.subs("t", tPlot.grid(tIt))),4)...
				/max(1,nthroot(norm(simData.("agent"+i).observer.S.at(0)-S.at(0)),4))*0.6);
		elseif tIt < numel(tPlot.grid)*3/4-5
			dataSim = cat(2, dataSim, nthroot(norm(simData.("agent"+i).observer.S.subs("t", tPlot.grid(tIt)) - S.subs("t", tPlot.grid(tIt))),4)...
				/max(1,nthroot(norm(simData.("agent"+i).observer.S.at(10.01)-S.at(10.01)),4))*0.6);
		else
			dataSim = cat(2, dataSim, 0);
		end
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "outputErrorData"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);
end

for i = 1:mas.N
	dataSim = [];
	signalName = [];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "y" + string(tIt)];
		newData = (simData.("agent"+i).controlOutput.subs("t", tPlot.grid(tIt)))...
			/max(abs(simData.("agent"+i).controlOutput));
		if abs(newData) > 0.1
			dataSim = cat(2, dataSim, newData);
		else
			dataSim = cat(2, dataSim, 0);
		end
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "normOutputErrorData"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);
end


dataSim = [];
signalName = [];
for tIt = 1 : numel(tPlot.grid)
    signalName = [signalName, "t" + string(tIt)];
    dataSim = cat(2, dataSim, tPlot.grid(tIt));
end % for tIt = tSample.'

data.t = export.dd(...
    "M", dataSim, ...
    "header", table2cell(array2table(signalName)).', ...
    "filename", "timeVec", ...
    "basepath", basepath);

for i = 1:mas.N
	d = aggregatedDisturbance((i-1)*length(P*v)+1:i*length(P*v));
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

data.g1 = export.dd(...
    "M", [zPlot.grid, g1.on(zPlot)], ...
    "header", ["z", "g1"], ...
    "filename", "g1", ...
    "basepath", basepath);

d = P*v;
data.dist = export.dd(...
    "M", [tPlot.grid, d.on(tPlot)], ...
    "header", ["t", "d"], ...
    "filename", "dist", ...
    "basepath", basepath,...
	"N", 400);

dataSim = tPlot.grid;
signalName = "t";
for i = 1:mas.N
	dataSim = cat(2, dataSim, simData.("agent" + i).observer.S.norm().on(tPlot));
	signalName = cat(2, signalName, "S" + i);
end
exportObj = export.dd(...
    "M", dataSim, ...
    "header", table2cell(array2table(signalName)).', ...
    "filename", "normSObs", ...
    "basepath", basepath);
exportObj.export("silent", true);

dataSim = tPlot.grid;
signalName = "t";
for i = 1:mas.N
	dataSim = cat(2, dataSim, simData.("agent" + i).observer.P.norm().on(tPlot));
	signalName = cat(2, signalName, "P" + i);
end
exportObj = export.dd(...
    "M", dataSim, ...
    "header", table2cell(array2table(signalName)).', ...
    "filename", "normPObs", ...
    "basepath", basepath);
exportObj.export("silent", true);

dataSim = tPlot.grid;
signalName = "t";
for i = 1:mas.N
	dataSim = cat(2, dataSim, simData.("agent" + i).observer.v(1).on(tPlot));
	signalName = cat(2, signalName, "v1" + i);
	dataSim = cat(2, dataSim, simData.("agent" + i).observer.v(2).on(tPlot));
	signalName = cat(2, signalName, "v2" + i);
end
exportObj = export.dd(...
    "M", dataSim, ...
    "header", signalName, ...
    "filename", "vObs", ...
    "basepath", basepath);
exportObj.export("silent", true);

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
        
export.Data.exportAll(data);

%% Simulation Führungsverhalten + Störverhalten
splitDomain = tDomain.split(simTime/2);
Pr = [1 0]/2;
freq = 1*pi;
freq2 = 1.5*pi;
Sr = quantity.Piecewise({[0 freq; -freq 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 freq2; -freq2 0] * quantity.Discrete.ones(1, splitDomain(2))});
Sr = quantity.Discrete(Sr);
w = quantity.Piecewise({quantity.Discrete([sin(freq*splitDomain(1).grid), cos(freq*splitDomain(1).grid)], splitDomain(1)), ...
    quantity.Discrete([sin(freq2*splitDomain(2).grid), cos(freq2*splitDomain(2).grid)], splitDomain(2))});
w = quantity.Discrete(w);
freq = 0.6*pi;
freq2 = 1*pi;
P = quantity.Piecewise({[6 0] * quantity.Discrete.ones(1, splitDomain(1)), [4 0] * quantity.Discrete.ones(1, splitDomain(2))});
P = quantity.Discrete(P);
S = quantity.Piecewise({[0 freq; -freq 0] * quantity.Discrete.ones(1, splitDomain(1)), [0 freq2; -freq2 0] * quantity.Discrete.ones(1, splitDomain(2))});
S = quantity.Discrete(S);
v = quantity.Piecewise({quantity.Discrete([sin(freq*splitDomain(1).grid), cos(freq*splitDomain(1).grid)], splitDomain(1)), ...
    quantity.Discrete([sin(freq2*splitDomain(2).grid), cos(freq2*splitDomain(2).grid)], splitDomain(2))});
v = quantity.Discrete(v);

mas = MultiAgentHeterogenous(agents, adjacencyMatrix, "nCoefDist", 23, "nCoefRef", 23);


aggregatedDisturbance = [P*v; P*v; P*v; P*v];

obsGain = 2;
arf = ReferenceObserver(adjacencyMatrix, 2, obsGain, obsGain, obsGain, obsGain, 1, obsGain, obsGain);
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
exogenousInput = {"disturbance", aggregatedDisturbance, "S", S, "Sr", Sr, "w", w, "Pr", Pr, "P", P};
simData = mas.adaptiveControlSimulation(simulationSetting, 2, 1, arf, "initialCondition", ic, "exogenousInput", exogenousInput, "mu_S", 2, ...
	"distCoef", 10, "adaptationDelayObserver", 10, "adaptationDelayController", 10);

for idx = 1:4
	simData.("agent" + idx).controlOutput.plot();
end

%% Export Führungsverhalten + Störverhalten
tPlot = quantity.Domain("t", linspace(0, simTime, 401));
zPlot = quantity.Domain("z", linspace(0, 1, 101));
basepath = "C:\Users\xzb84\Documents\Vorträge\GMA2022\mrmbeamer-main\dataCombi";
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
	signalName = [];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "y" + string(tIt)];
		newData = (simData.("agent"+i).controlOutput.subs("t", tPlot.grid(tIt))-r.subs("t", tPlot.grid(tIt)))...
			/max(abs(simData.("agent"+i).controlOutput - r));
		if abs(newData) > 0.1
			dataSim = cat(2, dataSim, newData);
		else
			dataSim = cat(2, dataSim, 0);
		end
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "normOutputErrorData"+i, ...
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

for i = 1:mas.N
	dataSim = [];
	signalName = [];
    for tIt = 1 : numel(tPlot.grid)
        signalName = [signalName, "y" + string(tIt)];
		if tIt < 2
			dataSim = cat(2, dataSim, 0);
		elseif tIt < numel(tPlot.grid)/2
			dataSim = cat(2, dataSim, nthroot(norm(simData.("agent"+i).observer.Sr.subs("t", tPlot.grid(tIt)) - Sr.subs("t", tPlot.grid(tIt))),4)...
				/max(1,nthroot(norm(simData.("agent"+i).observer.Sr.at(0)-Sr.at(0)),4))*0.6);
		elseif tIt < numel(tPlot.grid)*3/4-5
			dataSim = cat(2, dataSim, nthroot(norm(simData.("agent"+i).observer.Sr.subs("t", tPlot.grid(tIt)) - Sr.subs("t", tPlot.grid(tIt))),4)...
				/max(1,nthroot(norm(simData.("agent"+i).observer.Sr.at(10.01)-Sr.at(10.01)),4))*0.6);
		else
			dataSim = cat(2, dataSim, 0);
		end
    end % for tIt = tSample.'

    exportObj = export.dd(...
        "M", dataSim, ...
        "header", table2cell(array2table(signalName)).', ...
        "filename", "outputErrorData"+i, ...
        "basepath", basepath);
    exportObj.export("silent", true);
end

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

        
export.Data.exportAll(data);