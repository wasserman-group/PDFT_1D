function out = v_ee_DFT(N,h,lambda)
    c = 1.0;
    dindex = linspace(0,N-1,N);
    dindex = dindex';
    dr = dindex*h;
    out = lambda*1./sqrt(c + (dr).^2);
end