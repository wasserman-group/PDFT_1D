function out = v_ext_exp(x,R)
% exponential potential for H atom
    A = 1.071295;
    k = 1.0/2.385345;
    out = -1.0*A*exp(-k*abs(x-R));
end