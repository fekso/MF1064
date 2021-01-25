function xprim=numderiv(beta,gamma,alpha,t,omega)

N = length(alpha);


for k = 2:(N-1)
    dt = (alpha(k+1)-alpha(k))/(omega);
    dbeta = beta(k+1)-beta(k-1);
    betaprim(k) = dbeta/(2*dt);
    dgamma = gamma(k+1)-gamma(k-1);
    gammaprim(k) = dgamma/(2*dt);
end

betaprim(1) = NaN;
gammaprim(1) = NaN;
betaprim(N) = NaN;
gammaprim(N) = NaN;

xprim = [betaprim' gammaprim'];