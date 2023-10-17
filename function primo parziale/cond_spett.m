function [K]=cond_spett(A)
     K=max(abs(eig(A)))/min(abs(eig(A)));
end