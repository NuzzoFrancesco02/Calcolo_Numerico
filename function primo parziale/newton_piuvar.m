function [xmat,it]=newton_piuvar(x0,nmax,toll,fun,jac)
  it=0;
  xmat=[x0];
  xv=x0;
  err=toll+1;
      while(err>toll && it<nmax)
        jacobiano=jac(xv);
          if(det(jacobiano)==0)
                error('jacobiana singolare')
          end
        delta=jacobiano\(-fun(xv));
        xn=delta+xv;
        err=norm(xn-xv);
        xmat=[xmat xn];
        it=it+1;
        xv=xn;
      end
end
