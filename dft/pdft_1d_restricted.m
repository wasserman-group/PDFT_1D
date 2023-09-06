function [EAlpha,EfragAlpha,DensAlpha,totDens,totDensup,totDensdn,vp,vpH,vpext,vpXC,vpkin, ...
    Etot,totTs,totEext,totEH,totEXC,Ep,Epkin,Epext,EpH,EpXC,vpContr,optimality, ...
    DenspAlpha,DenspAlphaPlusOne,EpAlpha,EpAlphaPlusOne,muAlpha] = ...
pdft_1d_restricted(NAlpha,x,N,h,RCell,vCell,veeMatrix,tolvp,tolDFT,nGuess)

% Electron populations
Nfrag = length(vCell);
Nele = sum(NAlpha);
pAlpha = fix(NAlpha+1e-16);
pAlphaPlusOne = pAlpha + 1;
wAlpha = NAlpha - double(pAlpha);

% Calculation setup
vp = zeros(size(x));      % Initialization of vp

fprintf('\n');
fprintf('>> v_p optimization start\n');

maxIter = 8;             % max iterations for vp convergence
optimality = 1;
optimality_vp = 1;
iter = 1;
eta = 0.1;                 % vp optimization stepsize
err = zeros(maxIter,1);
err_vp = zeros(maxIter,1);
newDensf = zeros(size(x));
muAlpha = zeros(size(NAlpha));

% Fragment calculations
while (iter <= maxIter) && (optimality > tolvp) 
    fprintf('\n');
    fprintf('>> Fragment calculations for iteration %d: start\n',iter);
    [TsAlpha,EextAlpha,EHAlpha,EXCAlpha,EAlpha,DensAlpha,totDens,EpAlpha,EpAlphaPlusOne,DenspAlpha,DenspAlphaPlusOne,...
        phiupNorpAlpha,phidnNorpAlpha,phiupNorpAlphaPlusOne,phidnNorpAlphaPlusOne, ...
        euppAlpha,ednpAlpha,euppAlphaPlusOne,ednpAlphaPlusOne] = fragment_calculator(RCell,vCell,pAlpha,wAlpha,pAlphaPlusOne,x,vp,veeMatrix,nGuess,tolDFT);
    
    dn = newDensf - totDens;
    newDensf = totDens;
    err(iter) = (1/(Nele^2))*sum( (dn).^2 )*h;
    optimality = err(iter);

    [vpnew,totTs] = vp_updater_vXC(x,h,N,RCell,vCell,totDens,vp);
    err_vp(iter) = sum( (vp - vpnew).^2 )*h;
    optimality_vp = err_vp(iter);
    vp = vpnew;
    [vpkin,vpext,vpH,vpXC] = vp_components(Nfrag,x,N,h,RCell,vCell,veeMatrix,vp,totDens,DenspAlpha,DenspAlphaPlusOne,wAlpha);
    totDensup = totDens/2;
    totDensdn = totDens/2;
    [Ep,Epkin,Epext,EpH,EpXC,totEext,totEH,totEXC] = decompose(vp,EAlpha,DensAlpha,totDens,totDensup,totDensdn,totTs,x,h,RCell,vCell,veeMatrix);
    [EfragAlpha,Etot,vpContr] = Energies(EAlpha,Ep,vp,DensAlpha,h);

    fprintf('\n');
    fprintf('>> Fragment calculations for iteration %d: done\n',iter);
    fprintf('>> ***Current optimality of density = %5.2e***\n',optimality);
    fprintf('>> ***Current optimality of vp = %5.2e***\n',optimality_vp);
    fprintf('>> ***Current max of abs err = %5.2e***\n',max(abs(dn)));
    fprintf('\n');

    iter = iter + 1;

end

for i = 1:Nfrag
    if pAlpha(i) >= 1
        muAlpha(i) = EpAlphaPlusOne(i) - EpAlpha(i);
    else
        muAlpha(i) = EpAlphaPlusOne(i);
    end
end

fprintf('>> =-=-=-=-=-=-=-=-=-=-=-=\n');
fprintf('>> PDFT exact module: done\n');
fprintf('>> =-=-=-=-=-=-=-=-=-=-=-=\n');

end

function [EfragAlpha,Etot,vpContr] = Energies(EAlpha,Ep,vp,DensAlpha,h)
    dE = NaN*ones(size(EAlpha));
    for i = 1:length(EAlpha)
        dE(1,i) = dot(DensAlpha(:,i),vp)*h;
    end
    EfragAlpha = EAlpha - dE;              % exclude the contribution from vp to fragments
    Etot = sum(EfragAlpha) + Ep;
    vpContr = dE;
end
