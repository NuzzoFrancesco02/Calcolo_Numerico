function [alpha]=alpha_opt(A,P)
 if nargin==1
     P=eye(size(A,1));
 end
   alpha=2/(max(eig(inv(P)*A))+min(eig(inv(P)*A)));
end