clear
clc
A = [];
P = [];
alpha = 1;

E = eig(P^-1*A);
flag = 0;
left = alpha.*abs(E).^2
right = 2*real(E)
for i=1:size(A)
	if  left(i) < right(i)
		flag = 1;
		fprintf("Richardson stazionario non converge\n");
	end
end

if flag == 0
	fprintf("Richardson stazionario converge\n");
end

%richardson(A,[1,1]',P,[1,1]',1e-6,100000,alpha)

% Se gli autovalori sono tutti reali basta verificare
% 0 < alpha*E < 2
if (0 < alpha*E) && (alpha*E < 2)
	fprintf("Check 2 (lambda reali) passato\n");
else
	fprintf("Check 2 (lambda reali) non passato\n");
end