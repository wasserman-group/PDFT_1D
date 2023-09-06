function [outEp,outEpkin,outEpext,outEpH,outEpxc] = calc_Ep(vpkin,vpext,vpH,vpXC,TsAlpha,EextAlpha,EHAlpha,ExcAlpha, ...
                                                           DensAlpha,totDens,totDensup,totDensdn,totTs, ...
                                                           x,h,RCell,vCell,veeMatrix)
% This functions calculates the partition energy: 
% Ep = Epkin + Epext + EpH + Epxc

% Exclude the contribution from vp in E_f
    dTs = NaN*ones(size(TsAlpha));
    dEext = NaN*ones(size(EextAlpha));
    dEH = NaN*ones(size(EHAlpha));
    dExc = NaN*ones(size(ExcAlpha));
    for i = 1:length(TsAlpha)
        dTs(1,i) = dot(DensAlpha(:,i),vpkin)*h;
        dEext(1,i) = dot(DensAlpha(:,i),vpext)*h;
        dEH(1,i) = dot(DensAlpha(:,i),vpH)*h;
        dExc(1,i) = dot(DensAlpha(:,i),vpXC)*h;
    end

% Non-additive kinetic energy Epkin
TsfragAlpha = TsAlpha - dTs;
Epkin = totTs - sum(TsfragAlpha);

% Non-additive nuclear-electronic (external) energy Epext
vext = build_vext(x,RCell,vCell,zeros(size(x)));
EextMol = trapz(x,totDens.*vext);
EextfragAlpha = EextAlpha - dEext;
Epext = EextMol - sum(EextfragAlpha);

% Non-additive electron-electron repulsion (Hartree) energy EpH
vH = v_H(totDens,x,veeMatrix);
EHMol = (1/2)*trapz(x,totDens.*vH);
EHfragAlpha = EHAlpha - dEH;
EpH = EHMol - sum(EHfragAlpha);

% Non-additive exchange-correlation energy Epxc
ExMol = E_x(totDens,totDensup,totDensdn,x);
EcMol = E_c(totDens,totDensup,totDensdn,x);
ExcMol = ExMol + EcMol;
ExcfragAlpha = ExcAlpha - dExc;
Epxc = ExcMol - sum(ExcfragAlpha);

% Output
outEpkin = Epkin;
outEpext = Epext;
outEpH = EpH;
outEpxc = Epxc;
outEp = Epkin + Epext + EpH + Epxc;
