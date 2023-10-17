function alpha_opt=alphaopt(A,P)
alpha_opt=2/(max(abs(eig(P\A)))+min(abs(eig(P\A))));
end