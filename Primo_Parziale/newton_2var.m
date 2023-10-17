%   function [xmat,it]=newton_2var(x0,nmax,toll,fun,jac)

%   ATTENZIONE!!! xmax parte da x(0)
function [xmat,it]=newton_2var(x0,nmax,toll,fun,jac)
  it=0;
  xmat=[x0];
  xv=x0;
  err=toll+1;
      while(err>toll && it<nmax)
        jacobiano=jac(xv(1),xv(2));
          if(det(jacobiano)==0)
                error('jacobiana singolare')
          end
        delta=jacobiano\(-1*fun(xv(1),xv(2)));
        xn=delta+xv;
        err=norm(xn-xv);
        xmat=[xmat xn];
        it=it+1;
        xv=xn;
      end
end
