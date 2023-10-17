function [alpha]=alpha(A,r,z)
  if (nargin==2)
      z=r;
  end
 alpha=(z'*r)/(z'*A*z);
end
