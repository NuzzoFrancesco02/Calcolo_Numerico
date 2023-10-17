%% Autovalutazione 5
%% 1
A = [4 2 1; 2 4 1; 1 1 7];
b = [2 2 2]';
[x k r alpha] = richardson(A,b,eye(3),b,1e-6,1);
alpha
x
%% 2
A = hilb(3);
d = (cond(A)-1)/(cond(A)+1); %P = I == richardson non precondizionato se sdp
k = log(1/200)/log(d)
%% 3
clear
clc
g = sym('g', 'real');

A = [4 -1; -1 g];
b = [1, 1]';
x0 = zeros(2,1);
[x, k, p] = conjgrad_it(A,b,x0,1,1e-6);
(simplify(p))

%% 4
A = [8 -2 -2; -2 6 -1; -2 -1 9];
b = ones(3,1);
[x, flag, r, k] = pcg(A,b,1e-6,2,[],[],b);
x 
k
%% 5
%A*v = lambda*v
syms b;
A = [10 -2; -2 b];
v = [1 1]';

lambda = eig_approx(A,v);

simplify(lambda)

%%
M = diagonals([1 -8 14 -8 1],10)

