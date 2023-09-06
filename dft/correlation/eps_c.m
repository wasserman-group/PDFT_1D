function out = eps_c(n,nup,ndn)
    rs = 1./(2*n);
    ksi = ksi_f(n,nup,ndn);
    e0 = e0_f(rs);
    e1 = e1_f(rs);
    out = e0 + (ksi.^2).*(e1 - e0);
end

function out = e0_f(rs)
    
    A = 18.40;
    B = 0.0;
    C = 7.501;
    D = 0.10185;
    E = 0.012827;
    a = 1.511;
    b = 0.258;
    m = 4.424;
    
    out = (-1/2)*(rs + E*rs.^2)./(A + B*rs + C*rs.^2 + D*rs.^3).*log(1 + a*rs + b*rs.^m);
end

function out = e1_f(rs)
    
    A = 5.24;
    B = 0.0;
    C = 1.568;
    D = 0.1286;
    E = 0.00320;
    a = 0.0538;
    b = 1.56e-5;
    m = 2.958;
    
    out = (-1/2)*(rs + E*rs.^2)./(A + B*rs + C*rs.^2 + D*rs.^3).*log(1 + a*rs + b*rs.^m);
end

function out = ksi_f(n,nup,ndn)
    out = (nup - ndn)./n;
end