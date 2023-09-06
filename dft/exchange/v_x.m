function out = v_x(n,nup,ndn,str)
% [A/B] = integral_dx{n^2*(1[+/-]sigma)^2*f(kF*(1[+/-]sigma))}
    kf = k_f(n);
    ksi = ksi_fun(n,nup,ndn)+eps;       %>why? >B/C eps is an eps ans grid is a grid
    z = double(kf.*(1 + ksi));
    zp = double(kf.*(1 - ksi));   %z-prime
    f = f_fun(z);
    fp = f_fun(zp);
    
    dfdz = df_dz(z);
    dfdzp = df_dz(zp);      %df(z-prime)/d(z-prime):function of z-prime
    
    if (strcmp(str,'up'))
        temp = 4*n.*(1 + ksi).*f + (pi)*(n.^2).*((1 + ksi).^2).*dfdz;
    elseif (strcmp(str,'dn'))
        temp = 4*n.*(1 - ksi).*fp + (pi)*(n.^2).*((1 - ksi).^2).*dfdzp;
    else
    end
    out = -1/4*(temp);
end

function out = f_fun(z)%<-Yan
    temp = MeijerG_2224(z.^2)./(4*z);
    out = temp;
end

function out = df_dz(z)%<-Yan

    z2 = z.^2;
    temp = (MeijerG_2113(z2) - MeijerG_2224(z2))...
            ./(2*z2);

    out = temp;
end

function out = k_f(n)
    out = n*(pi/2);
end

function out = ksi_fun(n,nu,nd)
    out = (nu - nd)./n;
end

function out = MeijerG_2224(z)
    out = meijerG([0.5,1.0],[],[0.5,0.5],[-0.5,0.0],z);
end

function out = MeijerG_2113(z)
    out = meijerG(1.0,[],[0.5,0.5],0.0,z);
end