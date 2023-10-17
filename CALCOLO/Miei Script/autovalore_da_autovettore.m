syms a;
A = [];
x=[]';
y=x/norm(x);
lambda = y'*A^-1*y
simplify(lambda)