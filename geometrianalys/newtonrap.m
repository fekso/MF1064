function  x = newtonrap(x,alpha)

L1=654;          %[mm]
L2=150;          %[mm]
L3=670;          %[mm]
L4x=132;         %[mm]
L4y=55;          %[mm]

beta = x(1);
gamma = x(2);

g=L1*cos(alpha)+L2*cos(gamma')-L3*cos(beta')-L4x;
h=L1*sin(alpha)-L2*sin(gamma')-L3*sin(beta')+L4y;
f=[g;h];

tol=10^-3; 
imax=1000; 
i=0;

while norm(f)>tol && i<imax
i = i+1;
A=[L3*sin(beta), -L2*sin(gamma);-L3*cos(beta),-L2*cos(gamma)];
x=x-A\f;
f1=L1*cos(alpha)+L2*cos(x(2))-L3*cos(x(1))-L4x;
f2=L1*sin(alpha)-L2*sin(x(2))-L3*sin(x(1))+L4y;
f=[f1;f2];
end

end
