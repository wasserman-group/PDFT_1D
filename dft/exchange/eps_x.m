function out = eps_x(n,nup,ndn)
    kf = k_f(n);
    sigma = ksi_fun(n,nup,ndn)+eps;
    
    z = kf.*(1 + sigma);
    zp = kf.*(1 - sigma);   %z-prime
    f = f_fun(z);
    fp = f_fun(zp);

    a = ((1 + sigma).^2).*f.*n;
    b = ((1 - sigma).^2).*fp.*n;

    out = (-1/4)*(a + b);
end

function out = k_f(n)
    out = n*(pi/2);
end

function out = ksi_fun(n,nu,nd)
    out = (nu - nd)./n;
end

function out = f_fun(z)
    temp = MeijerG_2224(z.^2)./(4*z);
    out = temp;
end

function out = MeijerG_2224(z)
    out = meijerG([0.5,1.0],[],[0.5,0.5],[-0.5,0.0],z);
end