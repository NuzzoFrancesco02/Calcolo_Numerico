function[xn,it]=newton_sislin(x0,nmax,tol,f,j)
   it=0;
   err=tol+1;
   xv=x0;
      while(err>tol && it<nmax)
          F=-1*f(xv);
          J=j(xv);
          delta=J\F;
          xn=xv+delta;
          err=norm(xn-xv);
          it=it+1;
          xv=xn;
      end
end