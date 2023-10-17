F=@(x,y,z) [sin(pi.*x)-y;
             y.^2-z;
             -x-y+z.^2];
J=@(x,y,z) [pi.*cos(pi.*x) -1 0;
            0 2.*y -1;
            -1 -1 2.*z];
x0=[ 1/5 1/5 1/5]';
if(det(J(x0(1),x0(2),x0(3)))==0)
    error('jacobiana singolare')
end

delta=J(x0(1),x0(2),x0(3))\(-F(x0(1),x0(2),x0(3)));
x1=delta+x0;
if(det(J(x1(1),x1(2),x1(3)))==0)
    error('jacobiana singolare')
end
delta=J(x1(1),x1(2),x1(3))\(-F(x1(1),x1(2),x1(3)));
x2=delta+x1;
if(det(J(x2(1),x2(2),x2(3)))==0)
    error('jacobiana singolare')
end
