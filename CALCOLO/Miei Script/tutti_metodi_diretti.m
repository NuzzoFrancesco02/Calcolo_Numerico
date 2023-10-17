A = [];
b = [];

t = tic;
x_back = A \ b;
T (1, i) = toc (t);
t = tic;
x_back_2 = (-A)\(-b);
T (2, i) = toc (t);

t = tic;
[L, U, P] = lu (A);
x_lu = U \ (L \ (P*b));
T (3, i) = toc (t);

t = tic;
[L, U, x_thom] = thomas (A, b);
T (4, i) = toc (t);

A_chol = -A;
t_chol = -b;
t = tic;
H = chol (A_chol);
x_chol = H \ (H' \ t_chol);
T (5, i) = toc (t);