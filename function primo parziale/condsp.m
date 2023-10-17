function Ksp=condsp(A)
Ksp=(max(abs(eig(A))))/(min(abs(eig(A))));
end
