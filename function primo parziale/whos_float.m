function [min,max,card,eps]=whos_float(beta,t,L,U)
min=beta^(L-1);
max=(1-beta^-t)*beta^U;
card=2*(beta-1)*beta^(t-1)*(U-L+1)+1;
eps=beta^(1-t);
end