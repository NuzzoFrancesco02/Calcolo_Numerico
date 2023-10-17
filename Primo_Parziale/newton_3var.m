function [xmat,it]=newton_3var(x0,nmax,toll,fun,jac)
  it=0;
  xmat=[x0];
  xv=x0;
  err=toll+1;
      while(err>toll && it<nmax)
        jacobiano=jac(xv(1),xv(2),xv(3));
          if(det(jacobiano)==0)
                error('jacobiana singolare')
          end
        delta=jacobiano\(-fun(xv(1),xv(2),xv(3)));
        xn=delta+xv;
        err=norm(xn-xv);
        xmat=[xmat xn];
        it=it+1;
        xv=xn;
      end
end
