function K = powercond(A,x0,nmax)
%
% POWERCOND Metodo per il calcolo del numero di condizionamento spettrale
% tramite metodo delle potenze dirette e inverse
% PER MATRICI SDP
% Parametri di ingresso:
%
% A: matrice di sistema
% nmax: numero massimo di iterazioni ammesse
%
% Parametri di uscita:
%
% K: numero di condizionamento spettrale di A
%
% Inizializzazione

if isequal(A,A') && all(eig(A))
else
	error('A non è SDP!');
end

if nargin < 3
	nmax = 100;
end

n = size(A,2);
iter = 0;
ymax = x0 / norm(x0);
ymin = ymax;
K = 1;
while iter < nmax
	iter = iter + 1;
	xmax = A * ymax;				% 2nˆ2−n ops
	ymax = xmax / norm(xmax);		% 2n+1 ops
	lambdamax = ymax' * A * ymax;	% 2nˆ2+n−1 ops
	xmin = A \ ymin;
	ymin = xmin / norm(xmin);				% 2n+1 ops
	lambdamin = ymin' * A * ymin;	% 2nˆ2+n−1 ops
	K_new = lambdamax / lambdamin;
	K = K_new;
end

fprintf("\nStima numero di condizionamento spettrale (nnmax = %d)", nmax);
K
end