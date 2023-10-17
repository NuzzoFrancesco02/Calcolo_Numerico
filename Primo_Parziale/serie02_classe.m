A = [[50 1 3]; [1 6 0]; [3 0 1]];
B = [[50 1 10]; [3 20 1]; [10 4 70]];
C = [[7 8 9]; [5 4 3];[1 2 6]];

%Condizioni sufficienti
% if ( A == A')
%     n = size(A,1);
%     for i = 1 : n
%         if det(A(i:n,i:n))<=0
%             error('La matrice non è def positiva');
%         end
%     end
%     disp('La matrice è definita positiva');
% end
% 
% d_B = diag(B);
% n = size (B,1);
% for j = 1:n
%     r(j) = abs(d_B(j)) - sum( abs( B([1:j-1,j+1:n],j) ) );
% end
% if r > 0
%     disp('La matrice è a dominanza stretta per colonne');
% end
% 
% for i = 1 :n-1
%     if det(C(1 : i, 1 : i)==0)
%         error('La condizione nec e suff non è soddisfatta');
%     end
%     disp('La condizione nec e suff non è soddisfatta');
% end







[L_A, U_A] = lugauss(A);
[L_B, U_B] = lugauss(B);
[L_C, U_C] = lugauss(C);
x_es = [1 1 1]';
b = A*x_es;

%LU*x=b L*y = b U*x = y
y = fwsub(L_A,b);
x = bksub(U_A,y);
err_rel = norm(x_es-x)/norm(x_es);
r_norm = norm(b - A*x)/norm(b);

%% Es 1.3
clear
clc
n = 20;
I0=2;
R = 1;
b = [I0 zeros(1,n-1)]';
A = diag(-R*ones(1,n))+diag(R*ones(1,n-1),-1);
A(1,:) = 1;
A(2,1)=1e3;
[L_A, U_A, P]=lu(A);
if(P==eye(n))
    disp('non è stato utilizzato il pivoting');
else
    disp('è stato effettuato pivoting');
end
y = fwsub(L_A,P * b);
x = bksub(U_A,y);
i_ex = ones(n,1)*I0/n;

err_rel = norm(i_ex -x)/norm(i_ex);
res_norm = norm(b - A*x)/norm(b);
cond(A);
 %% Es 2.1
 n = 1000;

 A = hilb(n);
 x_ex = ones(n,1);
 b = A*x_ex;

 B = rand(n);
 y_ex = ones(n,1);
 c = B * y_ex;

 x = A \ b;
 disp('---Matrice A---');
 fprintf('cond(A) = %5.4e \n\n', cond(A));
 fprintf('errore relativo = %5.4e\n\n', norm(x-x_ex)/norm(x_ex));
 fprintf('residuo normalizzato = %5.4e\n\n', norm(b-A*x)/norm(b));

 y = B \ c;
  disp('---Matrice B---');
 fprintf('cond(B) = %5.4e \n\n', cond(B));
 fprintf('errore relativo = %5.4e\n\n', norm(y-y_ex)/norm(y_ex));
 fprintf('residuo normalizzato = %5.4e\n\n', norm(c-B*y)/norm(c));