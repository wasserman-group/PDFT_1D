function [Ep,Epkin,Epext,EpH,EpXC,totEext,totEH,totEXC] = decompose(vp,EAlpha,DensAlpha,totDens,totDensup,totDensdn,totTs, ...
                                                           x,h,RCell,vCell,veeMatrix)
                                                       
                                                       
Nfrag = size(DensAlpha,2);
                                                       
% Non-additive nuclear-electronic (external) energy Epext
vext = build_vext_partition(x,RCell,vCell,zeros(size(x)));
vext = cell2mat(vext);
vextMol = sum(vext,2);
EextMol = trapz(x,totDens.*vextMol);
EextfragAlpha = sum(vext.*DensAlpha,1)*h;
Epext = EextMol - sum(EextfragAlpha);

% Non-additive electron-electron repulsion (Hartree) energy EpH
vH = v_H(totDens,x,veeMatrix);
EHMol = (1/2)*trapz(x,totDens.*vH);
vHfrag = zeros(size(DensAlpha));
for i = 1:Nfrag
    vHfrag(:,i) = v_H(DensAlpha(:,i),x,veeMatrix);
end
EHfragAlpha = (1/2)*sum(vHfrag.*DensAlpha,1)*h;
EpH = EHMol - sum(EHfragAlpha);

% Non-additive exchange-correlation energy Epxc
ExMol = E_x(totDens,totDensup,totDensdn,x);
EcMol = E_c(totDens,totDensup,totDensdn,x);
EXCMol = ExMol + EcMol;
ExfragAlpha = zeros(1,Nfrag);
EcfragAlpha = zeros(1,Nfrag);
for i = 1:Nfrag
    ExfragAlpha(i) = E_x(DensAlpha(:,i),DensAlpha(:,i),zeros(size(x)),x);
    EcfragAlpha(i) = E_c(DensAlpha(:,i),DensAlpha(:,i),zeros(size(x)),x);
end
EXCfragAlpha = ExfragAlpha + EcfragAlpha;
EpXC = EXCMol - sum(EXCfragAlpha);

% Non-additive exchange-correlation energy Epxc
dE = NaN*ones(size(EAlpha));
for i = 1:length(EAlpha)
    dE(1,i) = dot(DensAlpha(:,i),vp)*h;
end
EfragAlpha = EAlpha - dE;
TsfragAlpha = EfragAlpha - EextfragAlpha - EHfragAlpha - EXCfragAlpha;
Epkin = totTs - sum(TsfragAlpha);

Ep = Epkin + Epext + EpH + EpXC;
totEext = EextMol;
totEH = EHMol;
totEXC = EXCMol;

end



