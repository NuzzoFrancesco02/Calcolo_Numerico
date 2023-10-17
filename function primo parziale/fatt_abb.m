function[d]=fatt_abb(A,P)
  if nargin==1
      P=eye(size(A,1));
  end
  d=(cond(inv(P)*A)-1)/(cond(inv(P)*A)+1);
end