function[K] = spectralcond(A,p)

%Calcola il numero di condizionamento della matrice A

% SDP check
if isequal(A,A') && all(eig(A))
		fprintf('\nA è SDP, K = K2\n');
else
	fprintf('\nA non è SDP, K != K2\n');
end


if nargin == 1
	K = max(abs(eig(A)))/min(abs(eig(A)));
else
	K = norm(A,p)*norm(inv(A),p);
end
end