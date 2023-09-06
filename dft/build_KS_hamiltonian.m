function out = build_KS_hamiltonian(N,h,arrayV)

    K = build_Kin6(N,1,h);

    H = K + spdiags(arrayV,0,N,N);

    out = H + eps;

end