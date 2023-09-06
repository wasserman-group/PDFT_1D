function out = E_c(n,nup,ndn,x)
    eps = eps_c(n,nup,ndn);   
    temp = trapz(x,n.*eps);
    out = temp;
end