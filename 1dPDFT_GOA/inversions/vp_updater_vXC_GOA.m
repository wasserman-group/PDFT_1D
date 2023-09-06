function [vpnew,totTs] = vp_updater_vXC_GOA(x,h,N,RCell,vCell,totDens,DensAlpha,S,dSdnAlpha,vp,EpXC,eta)

vextMatrix = build_vext(x,RCell,vCell,zeros(size(x)));
totDensup = totDens/2;
totDensdn = totDens/2;
Nele = length(vCell);
[totTs,vXC_tot_inv,vXC_tot_fun] = wuyang_vXC(totDens,totDensup,totDensdn,x,h,N,Nele,vextMatrix);
clear('vextMatrix');

vextMatrix = build_vext_partition(x,RCell,vCell,zeros(size(x)));
Nfrag = size(DensAlpha,2);
vXC_frag_inv = zeros(size(DensAlpha));
vXC_frag_fun = zeros(size(DensAlpha));

for i = 1:Nfrag
    [~,vXC_frag_inv(:,i),vXC_frag_fun(:,i)] = wuyang_vXC(DensAlpha(:,i),DensAlpha(:,i)/2,DensAlpha(:,i)/2,x,h,N,Nele/Nfrag,vextMatrix{i});
end

vpnew = vp + eta*((1-S)*sum(vXC_frag_fun.*DensAlpha,2)./(totDens+1e-12) - vXC_tot_inv + S*vXC_tot_fun + dSdnAlpha./(totDens+1e-10)*(EpXC));
end
