function simData = simulateOde(adjacencyMatrix, stateSpace, tDomain, NameValue)
%SIMULATEODE Summary of this function goes here
%   Detailed explanation goes here
arguments
	adjacencyMatrix double {mustBe.quadraticMatrix(adjacencyMatrix)};
	stateSpace ss;
	tDomain quantity.Domain;
	NameValue.Delay = 0;
	NameValue.InitialCondition cell = {};
	NameValue.informedAgents (:,1) double = adjacencyMatrix(:,1);
end

%Compute the Laplacian and Leader-Follower Matrix
L =  diag(sum(adjacencyMatrix, 2)) - adjacencyMatrix;
H_ext = L + diag(NameValue.informedAgents);
%Initialize state variable
simData = struct();
for indAgent = 1:length(adjacencyMatrix)
	if ~isempty(NameValue.InitialCondition)
		simData.("agent"+indAgent) = NameValue.InitialCondition{indAgent};
	else
		simData.("agent"+indAgent) = zeros(length(stateSpace.A));
	end
end
y = stateSpace.C * x + stateSpace.D * u(:,1);

%Run simulation in a loop
tStep = tDomain.upper / (tDomain.n - 1);
for ind = 2:tDomain.n
	if ind > NameValue.Delay
		x = cat(2, x, tStep * (stateSpace.A*x + stateSpace.B*u(:, ind - NameValue.Delay - 1)) + x(:, end));
	else
		x = cat(2, x, tStep * (stateSpace.A*x) + x(:, end));
	end
end

end

