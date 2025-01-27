t_sim = 6;
t_n = 601;
tDomain = quantity.Domain("t", linspace(0, t_sim, t_n));

A = [1, 2, 1, 6; 0, 1, 0, 0; 3, 0, 5, 1; 2, 4, 0, 1];
b = [1; 4; 2; 1];

L1 = diag([-10; -11; -12; -13]);
V1 = [-4, -1, -1, -3; -1, -4, -4, -2; -3, -1, -3, -1; -1, -2, -1, -4];
M1 = inv(V1)*L1*V1;
V2 = [-3, -2, -3, -1; -2, -3, -1, -4; -4, -1, -6, -1; -3, -2, -1, -3];
L2 = diag([-35; -35; -38; -39]);
M2 = V2*L2/V1;
c1 = [1; 2; -1; 1];
eM1t = [];
eM2t = [];
for idx = 1:tDomain.n
	eM1t = cat(1, eM1t, reshape(expm(M1*(tDomain.grid(idx)-1)), [1, length(M1), length(M1)]));
	eM2t = cat(1, eM2t, (expm(M2*tDomain.grid(idx))*c1).');
end
% eA = -eye(length(A))/reshape(eM1t(51,:,:),[4, 4])*quantity.Discrete(eM1t, tDomain);
eA = eye(length(A))*quantity.Discrete(eM1t, tDomain);
eb = -quantity.Discrete(eM2t*10, tDomain);

A_hat = A - eA;
b_hat = b - eb;

% for idx = 1:tDomain.n
% 	if abs(det(A_hat.atIndex(idx))) < 10^(-8)
% 		idx
% 	end
% end

T = 25;
x_hat = [];
x_hatT = zeros(length(A), 1);
for idx = 1:tDomain.n
	if mod(idx-1, T) == 0 && abs(det(A_hat.atIndex(idx))) > 10^(-8)
		x_hatT = A_hat.atIndex(idx)\b_hat.atIndex(idx);
	end
	x_hat = cat(1, x_hat, x_hatT.');
end
x_hat = quantity.Discrete(x_hat, tDomain);
% x_hat.plot();

x = A\b;
ex = x - x_hat;
ex.plot();
ex.at(t_sim)

basepath = "C:\Users\xzb84\Documents\Überblick\AdaptiverStörbeobachter\zeitvarianteLineareGleichung";

% Export control output with reference
	header = {'t'};
	M = [tDomain.grid];
	for idx = 1:length(A)
		header = cat(2, header, ("y"+idx));
		M = cat(2, M, ex(idx).valueDiscrete);
	end
	data.controlOutputRef = export.dd(...
	    'M', M, ...
	    'header', header, ...
	    'filename', 'xError', ...
	    'basepath', basepath ...
	    );
export.Data.exportAll(data);