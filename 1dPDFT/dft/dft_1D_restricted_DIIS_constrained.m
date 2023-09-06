function  [outTsup,outTsdn,outEext,outEH,outEx,outEc,outEtot,outnup,outndn,...
                        outvxup,outvxdn,outvcup,outvcdn,...
                        outvH,outvext,...
                        outphiupNor,outphidnNor,...
                        outphiupNorVirt,outphidnNorVirt,outepsup,outepsdn,outepsupVirt,outepsdnVirt]...
                        = dft_1D_restricted_DIIS_constrained(n0,Nele,veeMatrix,vextMatrix,x,tol)                
    fprintf('1-D DFT restricted constrained DIIS accelerated \n');
    %                
    if (Nele)
        %% 
        N = length(x);
        h = x(2) - x(1);
        if (tol >= 1)
            tol = 1e-4;
        end
        %
        optimality = 1;
        Nitermax = 150;
        %
        if mod(Nele,2)
            %odd
            Nup = 1 + (Nele-1)/2;
            Ndn = (Nele-1)/2;
            Norb = Nup;
        else
            %even
            Nup = (Nele)/2;
            Ndn = (Nele)/2;
            Norb = Nup;
        end 
        % Sorting out the initial guess
        [~,bn0] = size(n0);
        if (bn0 == 1)
            nup0 = Nup*n0/Nele;
            ndn0 = Ndn*n0/Nele;
        elseif (bn0 == 2)
            nup0 = n0(:,1);
            ndn0 = n0(:,2);
        else
            fprintf('>> dft_1D: iniitial density guess is unclear');
        end
        %
        nup = nup0;
        ndn = ndn0;
        n = nup0 + ndn0;     
        %
        Ni = 1;
        %
        rhoupArray(:,Ni) = nup;
        rhodnArray(:,Ni) = ndn;
        rhoArray(:,Ni) = n;
        %
        B = ones(1,1);
        Bup = B;
        Bdn = B;
        %
        while (optimality >= tol) && (Ni <= Nitermax)
            %
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            fprintf('>> dft_1D: iteration = %d  \n',Ni);
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            %
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            fprintf('First call to regular update \n');
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            %
            [Nrhoup,Nrhodn,Nrho,~] = ...
                                Reg_updater(rhoupArray(:,Ni), ...
                                            rhodnArray(:,Ni), ...
                                            rhoArray(:,Ni), ...
                                            Nup,Ndn,Norb, ...
                                            veeMatrix,vextMatrix, ...
                                            x,N,h);
            %
            YArray(:,Ni) = Nrho;
            YupArray(:,Ni) = Nrhoup;
            YdnArray(:,Ni) = Nrhodn;
            %
            clear('Nrhoup');
            clear('Nrhodn');
            clear('Nrho');
            %
            if (Nup)
                [nup,Bup,cup] = DIIS_updater_n(rhoupArray,YupArray,Bup,Ni,x);
            end
            if (Ndn)
                [ndn,Bdn,cdn] = DIIS_updater_n(rhodnArray,YdnArray,Bdn,Ni,x);
            end
            %
            n = nup + ndn;
            %
            clear('phiupNor');
            clear('phidnNor');
            %
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            fprintf('Second call to regular update \n');
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');                      
            %
            [rhoupArray(:,Ni+1),rhodnArray(:,Ni+1),rhoArray(:,Ni+1), ...
                phiNor,phiNorVirt,eps,epsVirt] = ...
                        Reg_updater(nup,ndn,n, ...
                        Nup,Ndn,Norb,veeMatrix,vextMatrix,x,N,h);
            %
            optimality1 = rhoArray(:,Ni+1) - rhoArray(:,Ni);
            optimality2 = trapz(x,optimality1.*optimality1);
            optimality3 = trapz(x,rhoArray(:,Ni).*rhoArray(:,Ni));
            optimality = optimality2/optimality3;
            %
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            fprintf('>> dft_1D: iteration = %d \n',Ni);
            fprintf('>> dft_1D: optimality = %e \n',optimality);
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            %
            %
            clear('cup');
            clear('cdn');
            %
            if (Nup) clear('nup'), end
            if (Ndn) clear('ndn'), end
            clear('n');
            %
            Ni = Ni + 1;
        end
        %
        fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
        fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
        fprintf('>> dft_1D: \n');
        fprintf('>> optimization done \n');
        fprintf('>> calculating final values \n');
        fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
        %
        nConv = rhoArray(:,Ni);
        nConvup = rhoupArray(:,Ni);
        nConvdn = rhodnArray(:,Ni);
        vxup = v_x(nConv,nConvup,nConvdn,'up');
        vxdn = v_x(nConv,nConvup,nConvdn,'dn');
        vcup = v_c(nConv,nConvup,nConvdn,'up');
        vcdn = v_c(nConv,nConvup,nConvdn,'dn');
        vH = v_H(nConv,x,veeMatrix);
        vext = v_ext(vextMatrix,zeros(size(x)));
        %
        [Tsup,Tsdn,Eext,EH,Ex,Ec] = ...
                                DFT_1D_energies( ...
                                phiNor(:,1:Nup),phiNor(:,1:Ndn), ...
                                nConvup,nConvdn,Nup,Ndn,x,vH,vext);
        %
        Etot = Tsup + Tsdn + Eext + EH + Ex + Ec;   
    else
        %
        nConvup = zeros(size(x));
        nConvdn = zeros(size(x));
        %
        vxup = zeros(size(x));
        vxdn = zeros(size(x));
        vcup = zeros(size(x));
        vcdn = zeros(size(x));
        vH = zeros(size(x));
        vext = v_ext(vextMatrix,zeros(size(x)));
        phiNor = zeros(length(x));
        phiNorVirt = zeros(length(x));
        eps = zeros(length(x));
        epsVirt = zeros(length(x));
        %
        Tsup = 0.0;
        Tsdn = 0.0;
        Eext = 0.0;
        EH = 0.0;
        Ex = 0.0;
        Ec = 0.0;
        Etot = Tsup + Tsdn + Eext + EH + Ex + Ec;
    end
    %    
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
    fprintf('1-D S-DFT energy breakdown\n');
    fprintf('Ts-up = %d\n',Tsup);
    fprintf('Ts-down = %d\n',Tsdn);
    fprintf('E-ext = %d\n',Eext);
    fprintf('E-Hartree = %d\n',EH);
    fprintf('Ex = %d\n',Ex);
    fprintf('Ec = %d\n',Ec);
    fprintf('Etotal = %d\n',Etot);
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
    %
    outnup = nConvup;
    outndn = nConvdn;
    outTsup = Tsup;
    outTsdn = Tsdn;
    outEext = Eext;
    outEH = EH;
    outEx = Ex;
    outEc = Ec;
    outEtot = Etot;
    %
    outvxup = vxup;
    outvxdn = vxdn;
    outvcup = vcup;
    outvcdn = vcdn;
    outvH = vH;
    outvext = vext;
    %
    outphiupNor = phiNor;   %Yan: Test orbitals!!
    outphidnNor = phiNor;
    %
    outphiupNorVirt = phiNorVirt;
    outphidnNorVirt = phiNorVirt;
    outepsup = eps;
    outepsdn = eps;
    outepsupVirt = epsVirt;
    outepsdnVirt = epsVirt;
end

function [nupNew,ndnNew,nNew,phiNor,phiNorVirt,eps,epsVirt] = ...
    Reg_updater(nupOld,ndnOld,nOld,Nup,Ndn,Norb,veeMatrix,vextMatrix,x,N,h)
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
    fprintf('DFT optimization: calculating regular update\n');
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
%     if Nup == Ndn
%         vxup = v_x(nOld,nupOld,ndnOld,'up');
%         vxdn = v_x(nOld,nupOld,ndnOld,'dn');
%         vcup = v_c(nOld,nupOld,ndnOld,'up');
%         vcdn = v_c(nOld,nupOld,ndnOld,'dn');
%         vc = (1/2)*(vcup + vcdn);
%         vx = (1/2)*(vxup + vxdn);
%     elseif Nup == (Nup+Ndn)
%         vx = v_x(nOld/2,nupOld/2,zeros(size(x)),'up');
%         vc = v_c(nOld/2,nupOld/2,zeros(size(x)),'up');
%     else
%         vxup = v_x(nOld,nupOld,ndnOld,'up');
%         vxdn = v_x(nOld,nupOld,ndnOld,'dn');
%         vcup = v_c(nOld,nupOld,ndnOld,'up');
%         vcdn = v_c(nOld,nupOld,ndnOld,'dn');
%         vc = (Nup/(Nup+Ndn))*vcup + (Ndn/(Nup+Ndn))*vcdn;
%         vx = (Nup/(Nup+Ndn))*vxup + (Ndn/(Nup+Ndn))*vxdn;
%     end
    vxup = v_x(nOld,nupOld,ndnOld,'up');
    vxdn = v_x(nOld,nupOld,ndnOld,'dn');
    vcup = v_c(nOld,nupOld,ndnOld,'up');
    vcdn = v_c(nOld,nupOld,ndnOld,'dn');
    vc = (1/2)*(vcup + vcdn);
    vx = (1/2)*(vxup + vxdn);
    vH = v_H(nOld,x,veeMatrix);
    vext = v_ext(vextMatrix,zeros(size(x)));

    arrayV = vext + vH + vx + vc;

    H = build_KS_hamiltonian(N,h,arrayV);

    [phi,e] = eig(full(H));
    e = diag(e);
    [e,ind] = sort(e);
    phi = phi(:,ind);
    clear('ind');

    % normalize occupied vectors
    %Yan: Debug!!
    %phiNor = NaN*ones(length(x),Norb);
    %for j=1:Norb
    %    Aup = sqrt( 1/trapz(x,conj(phi(:,j)).*phi(:,j)) );
    %    fprintf('>> dft_1D_restricted_DIIS_constrained: Aup: %d \n',Aup)
    %    phiNor(:,j) = Aup*phi(:,j);
    %    AupTemp = Aup;
    %    clear('Aup');
    %end
    
    Aup = trapz(x,conj(phi).*phi,1);
    Aup = sqrt(Aup);
    %fprintf('>> dft_1D_restricted_DIIS_constrained v2: Aup: %d \n',1./Aup)
    %phiNorTemp = phiNor;
    %clear('phiNor');
    phiNor = phi(:,1:Norb)./Aup(1:Norb);
    phiNorVirt = phi(:,Norb+1:end)./Aup(Norb+1:end);
    %clear('Aup');

    nupNew = KS_densities(phiNor,Nup);
    ndnNew = KS_densities(phiNor,Ndn);
    nNew = nupNew + ndnNew;
    
    eps = e(1:Norb);
    epsVirt = e(Norb+1:end);
    
end

function [outn,outB,outc] = DIIS_updater_n(rhoArray,NArray,B,i,x)
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
    fprintf('DFT optimization: DIIS step update\n');
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
    for j = 1:i
        igrd1 = rhoArray(:,i) - NArray(:,i);
        igrd2 = rhoArray(:,j) - NArray(:,j);
        igrl = trapz(x,igrd1.*igrd2);
        B(i,j) = igrl;
        if (i ~= j)
            B(j,i) = igrl;
        end
        clear('igrl');
    end

    c = constr_min_DIIS(B);
    
    c = repmat(c',size(rhoArray,1),1);
    
    n = sum(c.*rhoArray,2); 
    
    outn = n;
    outB = B;
    outc = c;
end

function out = constr_min_DIIS(B)
    [l,~] = size(B);
    A = -1*eye(l);
    Aeq = ones(1,l);
    b = zeros(l,1);
    beq = 1;
    
    fun = @(q) DIIS_G(q,B);
    
    x0 = (1/l)*ones(l,1);
    
    [x,fval,exitflag,~] = fmincon(fun,x0,A,b,Aeq,beq);
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
    fprintf('constrained DIIS optimization done\n');
    fprintf('exitflag = %d \n',exitflag);
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
    
    c = x;
    
    out = c;
end

function out = DIIS_G(c,B)
    temp = (c')*B*c;
    out = temp;
end