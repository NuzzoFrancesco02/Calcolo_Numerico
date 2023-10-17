function[X]=inversa(A)
  if(size(A,1)~= size(A,2))
      error('matrice non quadrata')
  end
  if(det(A)==0)
      error('matrice non invertibile')
  end
  [L,U]=lugauss(A);
  n=size(A,1);
  B=eye(n);
  for i=1:n
      X(:,i)=bksub(U,fwsub(L,B(:,i)));
  end
end

