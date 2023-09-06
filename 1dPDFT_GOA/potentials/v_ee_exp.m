function out = v_ee_exp(N,h)
    A = 1.071295;
    k = 1.0/2.385345;
    dindex = linspace(0,N-1,N);
    dindex = dindex';
    dr = dindex*h;
    out = A*exp(-k*abs(dr));
end  