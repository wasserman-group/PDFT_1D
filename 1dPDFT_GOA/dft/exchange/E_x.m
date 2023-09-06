function out = E_x(n,nup,ndn,x)   
    eps = eps_x(n,nup,ndn);
    temp = trapz(x,n.*eps);
    out = temp;
end