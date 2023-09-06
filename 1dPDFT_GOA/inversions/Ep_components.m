function [Ep_OA,Ep,Epkin,Epext,EpH,EpXC,totEext,totEH,totEXC] = Ep_components(vp,vpkin,EAlpha,TsAlpha, ...
                                                           DensAlpha,totDens,totDensup,totDensdn,totTs, ...
                                                           x,h,RCell,vCell,veeMatrix,S)
                                                       
                                                       
Nfrag = size(DensAlpha,2);


% Non-additive kinetic energy Ts^nad
dTs = NaN*ones(size(TsAlpha));
for i = 1:length(TsAlpha)
    dTs(1,i) = dot(DensAlpha(:,i),vpkin)*h;
end
TsfragAlpha = TsAlpha - dTs;
Epkin = totTs - sum(TsfragAlpha);
                                                       
% Non-additive nuclear-electronic (external) energy Eext^nad
vext = build_vext_partition(x,RCell,vCell,zeros(size(x)));
vext = cell2mat(vext);
vextMol = sum(vext,2);
EextMol = trapz(x,totDens.*vextMol);
EextfragAlpha = sum(vext.*DensAlpha,1)*h;
Epext = EextMol - sum(EextfragAlpha);

% Non-additive electron-electron repulsion (Hartree) energy EH^nad
vH = v_H(totDens,x,veeMatrix);
EHMol = (1/2)*trapz(x,totDens.*vH);
vHfrag = zeros(size(DensAlpha));
for i = 1:Nfrag
    vHfrag(:,i) = v_H(DensAlpha(:,i),x,veeMatrix);
end
EHfragAlpha = (1/2)*sum(vHfrag.*DensAlpha,1)*h;
EpH = EHMol - sum(EHfragAlpha);

% Non-additive exchange-correlation energy EXC^nad
ExMol = E_x(totDens,totDensup,totDensdn,x);
EcMol = E_c(totDens,totDensup,totDensdn,x);
EXCMol = ExMol + EcMol;
dE = NaN*ones(size(EAlpha));
for i = 1:length(EAlpha)
    dE(1,i) = dot(DensAlpha(:,i),vp)*h;
end
EfragAlpha = EAlpha - dE;
EXCfragAlpha = EfragAlpha - EextfragAlpha - EHfragAlpha - TsfragAlpha;
EpXC = EXCMol - sum(EXCfragAlpha);


Ep = Epkin + Epext + EpH + EpXC;
Ep_OA = Epkin + Epext + EpH + S*EpXC;
totEext = EextMol;
totEH = EHMol;
totEXC = EXCMol;

end


