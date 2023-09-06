function [outTsup,outTsdn,outEext,outEH,outEx,outEc] = ...
                            DFT_1D_energies(...
                            phiupNor,phidnNor,nup,ndn,Nup,Ndn,x,vH,vExt)
% expects normalized vectors
    % calculate Kin energies
    N = length(x);
    h = x(2) - x(1);
    K = build_Kin6(N,1,h);
    Tsup = 0.0;
    Tsdn = 0.0;
    
    if (Nup)
        for j = 1:Nup
            Tsup = Tsup + trapz( x,conj(phiupNor(:,j)).*(K*phiupNor(:,j)) );   
        end
    end
    
    if (Ndn)
        for j = 1:Ndn
            Tsdn = Tsdn + trapz( x,conj(phidnNor(:,j)).*(K*phidnNor(:,j)) );
        end
    end
    n = nup + ndn;
    
    Eext = trapz(x,n.*vExt);
    EH = (1/2)*trapz(x,n.*vH);
    Ex = E_x(n,nup,ndn,x);
    Ec = E_c(n,nup,ndn,x);
    
    outTsup = Tsup;
    outTsdn = Tsdn;
    outEext = Eext;
    outEH = EH;
    outEx = Ex;
    outEc = Ec;
end