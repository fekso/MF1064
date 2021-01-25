function  x = newtonrap(x,alpha)

L1=657;          %[mm]
L2=50;          %[mm]
L3=670;          %[mm]
L4x=35;         %[mm]
L4y=44;          %[mm]

beta = x(1);
gamma = x(2);

g=L1*cos(alpha)+L2*cos(gamma')-L3*cos(beta')-L4x;
h=L1*sin(alpha)-L2*sin(gamma')-L3*sin(beta')+L4y;
f=[g;h];

tol=10^-3; % Toleransen
imax=1000; % Max antal iterationer per värde, detta för att inte hamna i en oändlig loop
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
