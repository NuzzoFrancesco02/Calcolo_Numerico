function[massimo]=matrix_max_lugauss(M)
 A=size(M);
 if(A(1)~=A(2))
     error('inserire una matrice quadrata')
 end
 n=A(1);
 L=eye(n);
 massimo=0;
 for k=1:n-1
     for i=k+1:n
         L(i,k)=M(i,k)./M(k,k);
          for j=k+1:n
              M(i,j)=M(i,j)-L(i,k).*M(k,j);
          end
     end
      if(abs(max(max(M)))>massimo)
      massimo=abs(max(max(M)));
      end
 end
end

