function out = KS_densities(phiNor,Nele)
% expects normalized densities
    [a,~] = size(phiNor);
    den = zeros(a,1);
    for i = 1:Nele
        den = den + conj(phiNor(:,i)).*phiNor(:,i);
    end
    out = den;
end