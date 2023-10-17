function []=unicita_LU(A)
  for i=1:size(A,1)
      if(det(A(1:i,1:i))==0)
          error('la fattorizzazione LU senza pivoting non è unica\n')
      end
  end
  disp('la fattorizzazione LU senza pivoting è unica')
end