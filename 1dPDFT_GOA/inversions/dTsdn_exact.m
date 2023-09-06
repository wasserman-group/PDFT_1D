function out = dTsdn_exact(N,h,n,Nup,Ndn,phiup,phidn,eup,edn)

K = build_Kin6(N,1,h);
temp = 0.0;

for i = 1:Nup
    temp = temp + (phiup(:,i).*(K*phiup(:,i)) - eup(i)*(phiup(:,i).^2));
end

for i = 1:Ndn
    temp = temp + (phidn(:,i).*(K*phidn(:,i)) - edn(i)*(phidn(:,i).^2));
end

out = (1./(n+1e-20)).*temp;

