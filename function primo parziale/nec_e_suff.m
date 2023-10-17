function[]=nec_e_suff(A,n)
 i=2;
  while(det(A(1:i,1:i))~=0 && i<n)
      i=i+1;
  end

  if(i==n)
      disp('la matrice rispetta la condizione necessaria e sufficiente')
  else
      disp('la matrice non rispetta la condizione necessaria e sufficiente')
  end
end

    
 