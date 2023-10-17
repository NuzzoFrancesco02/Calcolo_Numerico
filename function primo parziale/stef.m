function [x,it]=stef(x0,f,tol,nmax)
   x=x0;
   err=tol+1;
   it=0;
   while(err>tol && it<nmax)
       xv=x;
       x=xv-((f(xv)).^2)/(f(xv+f(xv))-f(xv));
       err=abs(f(x));
       it=it+1;
   end
end
