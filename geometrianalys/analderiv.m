function xprim2 = analderiv(alpha,beta,gamma,L1,L2,L3,t,omega);

J11 = L3*sin(beta);
J12 = -L2*sin(gamma);
J21 = -L3*cos(beta);
J22 = -L2*cos(gamma);
N1 = L1*omega*sin(alpha); 
N2 = -L1*omega*cos(alpha);

J = [J11 J12;J21 J22];  
N = [N1;N2];

xprim2 = J\N;           

end