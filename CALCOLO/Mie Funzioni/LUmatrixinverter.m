function[A_inv] = LUmatrixinverter(A)

% Inverte matrici nxn tramite fattorizzazione LU

[L,U] = lugauss(A);

[n,m] = size(A);

if (n ~= m)
    error('A non è una matrice quadrata'); 
end

B=eye(n);
A_inv=[];

for i=1:3
	y = fwsub(L,B(:,i));
	x = bksub(U,y);
	A_inv=[A_inv x];
end
end