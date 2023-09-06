function out = vp_kin(x,h,N,RCell,vCell,vp,totDens,DenspAlpha,DenspAlphaPlusOne,pAlpha,wAlpha,pAlphaPlusOne, ...
    QpAlpha,QpAlphaPlusOne,phiupNorpAlpha,phidnNorpAlpha,euppAlpha,ednpAlpha, ...
    phiupNorpAlphaPlusOne,phidnNorpAlphaPlusOne,euppAlphaPlusOne,ednpAlphaPlusOne)

% Calculate dTs[n_f]/dn_f
Nfrag = length(vCell);
Nele = Nfrag;    % Modify as needed
vextMatrix = build_vext(x,RCell,vCell,vp);
totDensup = totDens/2;
totDensdn = totDens/2;
[~,outvksMol,~] = wuyang(totDens,totDensup,totDensdn,x,h,N,0,Nele,vextMatrix);
dTsdnf = -outvksMol;

% Fragment calculations
dTsdnpAlpha = NaN*ones(N,length(vCell));
dTsdnpAlphaPlusOne = NaN*ones(N,length(vCell));
dTsdnAlpha = zeros(N,1);

for i = 1:Nfrag
    fprintf('Ts inversion for Fragment %d. \n',i);
    if any(DenspAlpha(:,i))
        if mod(pAlpha(i),2)
            %odd
            Nup = 1 + (pAlpha(i)-1)/2;
            Ndn = (pAlpha(i)-1)/2;
        else
            %even
            Nup = (pAlpha(i))/2;
            Ndn = (pAlpha(i))/2;
        end
        dTsdnpAlpha(:,i) = dTsdn_exact(N,h,DenspAlpha(:,i),Nup,Ndn, ...
            phiupNorpAlpha{i},phidnNorpAlpha{i},euppAlpha{i},ednpAlpha{i});
    else
        dTsdnpAlpha(:,i) = zeros(1,N);
    end
    
    if any(DenspAlphaPlusOne(:,i))
        if mod(pAlphaPlusOne(i),2)
            %odd
            Nup = 1 + (pAlphaPlusOne(i)-1)/2;
            Ndn = (pAlphaPlusOne(i)-1)/2;
        else
            %even
            Nup = (pAlphaPlusOne(i))/2;
            Ndn = (pAlphaPlusOne(i))/2;
        end
        dTsdnpAlphaPlusOne(:,i) = dTsdn_exact(N,h,DenspAlphaPlusOne(:,i),Nup,Ndn, ...
            phiupNorpAlphaPlusOne{i},phidnNorpAlphaPlusOne{i},euppAlphaPlusOne{i},ednpAlphaPlusOne{i});
    else
        dTsdnpAlphaPlusOne(:,i) = zeros(1,N);
    end
    
    dTsdnAlpha = dTsdnAlpha + (1-wAlpha(i))*dTsdnpAlpha(:,i).*QpAlpha(:,i) + wAlpha(i)*dTsdnpAlphaPlusOne(:,i).*QpAlphaPlusOne(:,i);
end

out = dTsdnf - dTsdnAlpha;
end
    
    
        
        
        
        