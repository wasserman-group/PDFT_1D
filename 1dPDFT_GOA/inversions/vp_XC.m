function out = vp_XC(QpAlpha,QpAlphaPlusOne,...
    totDens,DenspAlpha,DenspAlphaPlusOne,wAlpha,Nfrag)


ntot = totDens;
ntotup = totDens/2;
ntotdn = totDens/2;

%
vpXCpAlpha = zeros(size(ntot));
vpXCpAlphaPlusOne = zeros(size(ntot));


% promolecular calculations
vXMol = (1/2)*(v_x(ntot,ntotup,ntotdn,'up') + v_x(ntot,ntotup,ntotdn,'dn'));
vCMol = (1/2)*(v_c(ntot,ntotup,ntotdn,'up') + v_c(ntot,ntotup,ntotdn,'dn'));
vXCMol = vXMol + vCMol;

% fragment calculations
for i = 1:Nfrag

    if any(DenspAlpha(:,i))
        vXpAlphaup = v_x(DenspAlpha(:,i),DenspAlpha(:,i)/2,DenspAlpha(:,i)/2,'up');
        vXpAlphadn = v_x(DenspAlpha(:,i),DenspAlpha(:,i)/2,DenspAlpha(:,i)/2,'dn');
        vXpAlpha = (1/2)*(vXpAlphaup + vXpAlphadn);
        vCpAlphaup = v_c(DenspAlpha(:,i),DenspAlpha(:,i)/2,DenspAlpha(:,i)/2,'up');
        vCpAlphadn = v_c(DenspAlpha(:,i),DenspAlpha(:,i)/2,DenspAlpha(:,i)/2,'dn');
        vCpAlpha = (1/2)*(vCpAlphaup + vCpAlphadn);
        vXCpAlpha = vXpAlpha + vCpAlpha;
        vpXCpAlpha = vpXCpAlpha + ...
            (1 - wAlpha(i)).*QpAlpha(:,i).*(vXCMol - vXCpAlpha);
    else
        fprintf('vp_XC: wrong input fragment density! \n');
    end

    if any(DenspAlphaPlusOne(:,i))
        vXpAlphaPlusOneup = v_x(DenspAlphaPlusOne(:,i),DenspAlphaPlusOne(:,i)/2,DenspAlphaPlusOne(:,i)/2,'up');
        vXpAlphaPlusOnedn = v_x(DenspAlphaPlusOne(:,i),DenspAlphaPlusOne(:,i)/2,DenspAlphaPlusOne(:,i)/2,'dn');
        vXpAlphaPlusOne = (1/2)*(vXpAlphaPlusOneup + vXpAlphaPlusOnedn);
        vCpAlphaPlusOneup = v_c(DenspAlphaPlusOne(:,i),DenspAlphaPlusOne(:,i)/2,DenspAlphaPlusOne(:,i)/2,'up');
        vCpAlphaPlusOnedn = v_c(DenspAlphaPlusOne(:,i),DenspAlphaPlusOne(:,i)/2,DenspAlphaPlusOne(:,i)/2,'dn');
        vCpAlphaPlusOne = (1/2)*(vCpAlphaPlusOneup + vCpAlphaPlusOnedn);
        vXCpAlphaPlusOne = vXpAlphaPlusOne + vCpAlphaPlusOne;
        vpXCpAlphaPlusOne = vpXCpAlphaPlusOne + ...
            wAlpha(i).*QpAlphaPlusOne(:,i).*(vXCMol - vXCpAlphaPlusOne);
    else
        fprintf('vp_XC: wrong input fragment density! \n');
    end
    
end

vpXC = vpXCpAlpha + vpXCpAlphaPlusOne;
out = vpXC;

end