clear
clc

% esercizio 1

%% punto 1 

A = [ 50 1 3; 1 6 0; 3 0 1 ]
B = [ 50 1 10; 3 20 1; 10 4 70 ]
C = [ 7 8 9; 5 4 3; 1 2 6 ]
 
% la matrice A e' simmetrica e definita positiva
if (A == A')
    disp('A e'' simmetrica')
    na = size(A, 1);
    for i=1:na
        if (det(A(1:i,1:i))>0)
            if (i==na) 
                disp('A e'' definita positiva') 
            end
        else
            error('A non e'' definita positiva')
        end   
    end
else
    disp('A non e'' simmetrica')
end

% la matrice B e' a dominanza diagonale stretta per colonne
db = diag(B);
nb = size(B,1);
for j = 1:nb
   r(j) = abs(db(j)) - sum(abs(B([1:j-1,j+1:nb], j)));
end
if (r > 0)
    disp('B e'' a dominanza diagonale per colonne')
end

% la matrice C non ha particolari proprieta'
% controllo che i suoi minori principali abbiano
% determinante diverso da zero

i = 1;
dc = 1;
nc = size(C,1);
while ( (i < nc) && (dc ~= 0) )
    dc = det(C(1:i, 1:i));
    i = i+1;        
end
if (dc == 0)
    disp('non posso applicare la decomposizione LU ad C')
else
    disp('posso applicare la decomposizione LU a C')
end  
  
%% punto 2
% vedi funzione lugauss.m
    
%% punto 3
% matrice A
[LA, UA] = lugauss(A)

% matrice B
[LB, UB] = lugauss(B)

% matrice B
[LC, UC] = lugauss(C)

%% punto 4
% vedi funzioni fwsub.m e bksub.m

%% punto 5
sol = ones(3,1);
b = A * sol;
y = fwsub(LA, b)
x = bksub(UA, y)

%% punto 6
err_rel = norm(sol - x)/norm(x)
res_nor = norm(A*x -b)/norm(b)