function [vpnew,totTs] = vp_updater_vXC(x,h,N,RCell,vCell,totDens,vp)

vextMatrix = build_vext(x,RCell,vCell,zeros(size(x)));
totDensup = totDens/2;
totDensdn = totDens/2;
Nele = length(vCell);
[totTs,vXC_inv,vXC_fun] = wuyang_vXC(totDens,totDensup,totDensdn,x,h,N,Nele,vextMatrix);
vpnew = vp + (vXC_fun - vXC_inv);

end
