function [outTsup,outTsdn,outEext,outEH,outEx,outEc,outEtot,outnup,outndn,...
                        outvxup,outvxdn,outvcup,outvcdn,...
                        outvH,outvext,...
                        outphiupNor,outphidnNor,...
                        outphiupNorVirt,outphidnNorVirt,outepsup,outepsdn,outepsupVirt,outepsdnVirt]...
                        = dft_1D_DIIS_constrained(n0,Nup,Ndn,veeMatrix,vextMatrix,x,tol)
    fprintf('*** 1-D DFT unrestricted constrained DIIS accelerated v0.1 ***\n');
    %                
    Nele = Nup + Ndn;   
    N = length(x);
    h = x(2) - x(1);
    %
    if (tol >= 1)
        tol = 1e-4;
    end
    %
    optimality = 1;
    Nitermax = 150;
    %
    if (Nele)
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
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            fprintf('>> dft_1D: iteration = %d  \n',Ni);
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            %
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            fprintf('First call to regular update \n');
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            %
            [Nrhoup,Nrhodn,Nrho,~,~] = ...
                Reg_updater(rhoupArray(:,Ni), ...
                            rhodnArray(:,Ni), ...
                            rhoArray(:,Ni), ...
                            Nup,Ndn,veeMatrix,vextMatrix,x,N,h);
            %
            YArray(:,Ni) = Nrho;    %Yan: Should actually fix it one day!
            YupArray(:,Ni) = Nrhoup;
            YdnArray(:,Ni) = Nrhodn;
            %
            clear('Nrhoup');
            clear('Nrhodn');
            clear('Nrho');
            %
            if (Nup)
                [nup,Bup,cup] = ...
                    DIIS_updater_n(rhoupArray,YupArray,Bup,Ni,x);
            end
            if (Ndn)
                [ndn,Bdn,cdn] = ...
                    DIIS_updater_n(rhodnArray,YdnArray,Bdn,Ni,x);
            end
            % nup or ndn is zero for 1-ele case
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
                phiupNor,phidnNor,...
                phiupNorVirt,phidnNorVirt,epsup,epsdn,epsupVirt,epsdnVirt] = ...
                                        Reg_updater(nup,ndn,n, ...
                                        Nup,Ndn,veeMatrix,vextMatrix,x,N,h);
            %
            optimality1 = rhoArray(:,Ni+1) - rhoArray(:,Ni);
            optimality2 = trapz(x,optimality1.*optimality1);
            optimality3 = trapz(x,rhoArray(:,Ni).*rhoArray(:,Ni));
            optimality = optimality2/optimality3;
            %
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
            fprintf('>> dft_1D: iteration = %d DONE \n',Ni);
            fprintf('>> dft_1D: optimality = %e \n',optimality);
            fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=\n');
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
        % Sorting out the external matrix
        [~,bvext] = size(vextMatrix);
        if (bvext == 1)
            vextup = v_ext(vextMatrix,zeros(size(x)));
            vextdn = v_ext(vextMatrix,zeros(size(x)));
        elseif (bvext == 2)
            vextup = v_ext(vextMatrix(:,1),zeros(size(x)));
            vextdn = v_ext(vextMatrix(:,2),zeros(size(x)));
        else
            fprintf('>> Reg_updater: vextMatrix is unclear');
        end
        %vext = v_ext(vextMatrix,zeros(size(x)));
        %
        [Tsup,Tsdn,Eext,EH,Ex,Ec] = ...
                                DFT_1D_energies(...
                                phiupNor,phidnNor,nConvup,nConvdn,Nup,Ndn,x,vH,v_ext(vextMatrix,zeros(size(x))));
        %                        phiupNor,phidnNor,nConvup,nConvdn,Nup,Ndn,x,vH,vext);
        Etot = Tsup + Tsdn + Eext + EH + Ex + Ec;
    else
        nConvup = zeros(size(x));
        nConvdn = zeros(size(x));
        %
        vxup = zeros(size(x));
        vxdn = zeros(size(x));
        vcup = zeros(size(x));
        vcdn = zeros(size(x));
        vH = zeros(size(x));
        % Sorting out the external matrix
        [~,bvext] = size(vextMatrix);
        if (bvext == 1)
            vextup = v_ext(vextMatrix,zeros(size(x)));
            vextdn = v_ext(vextMatrix,zeros(size(x)));
        elseif (bvext == 2)
            vextup = v_ext(vextMatrix(:,1),zeros(size(x)));
            vextdn = v_ext(vextMatrix(:,2),zeros(size(x)));
        else
            fprintf('>> Reg_updater: vextMatrix is unclear');
        end
        %
        %vext = v_ext(vextMatrix,zeros(size(x)));
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
    %
    outvext(:,1) = vextup;
    outvext(:,2) = vextdn;
    %outvext = vext;
    %
    outphiupNor = phiupNor; %Yan: Test orbitals!!--Huh?
    outphidnNor = phidnNor;
    %
    outphiupNorVirt = phiupNorVirt;
    outphidnNorVirt = phidnNorVirt;
    outepsup = epsup;
    outepsdn = epsdn;
    outepsupVirt = epsupVirt;
    outepsdnVirt = epsdnVirt;
end

function [nupNew,ndnNew,nNew,phiupNor,phidnNor,...
          phiupNorVirt,phidnNorVirt,epsup,epsdn,epsupVirt,epsdnVirt] = ...
    Reg_updater(nupOld,ndnOld,nOld,Nup,Ndn,veeMatrix,vextMatrix,x,N,h)
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
    fprintf('DFT optimization: calculating regular update\n');
    fprintf('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n');
%This is updater function that returns n_{i+1} density from n_i density
    vxup = v_x(nOld,nupOld,ndnOld,'up');
    vxdn = v_x(nOld,nupOld,ndnOld,'dn');
    vcup = v_c(nOld,nupOld,ndnOld,'up');
    vcdn = v_c(nOld,nupOld,ndnOld,'dn');
    vH = v_H(nOld,x,veeMatrix);
    % Sorting out the external matrix, for embedding vextup and vextdn may
    % be diffirent
    [~,bvext] = size(vextMatrix);
    if (bvext == 1)
        vextup = v_ext(vextMatrix,zeros(size(x)));
        vextdn = v_ext(vextMatrix,zeros(size(x)));
    elseif (bvext == 2)
        vextup = v_ext(vextMatrix(:,1),zeros(size(x)));
        vextdn = v_ext(vextMatrix(:,2),zeros(size(x)));
    else
        fprintf('>> Reg_updater: vextMatrix is unclear');
    end
    %vext = v_ext(vextMatrix,zeros(size(x)));

    arrayVup = vextup + vH + vxup + vcup;
    arrayVdn = vextdn + vH + vxdn + vcdn;
    %arrayVup = vext + vH + vxup + vcup;
    %arrayVdn = vext + vH + vxdn + vcdn;

    Hup = build_KS_hamiltonian(N,h,arrayVup);
    Hdn = build_KS_hamiltonian(N,h,arrayVdn);

    [phiup,eup] = eig(full(Hup));
    eup = diag(eup);
    [eup,ind] = sort(eup);
    phiup = phiup(:,ind);
    clear('ind');
    [phidn,edn] = eig(full(Hdn));
    edn = diag(edn);
    [edn,ind] = sort(edn);
    phidn = phidn(:,ind);
    clear('ind');

    % normalize occupied vectors
    %phiupNor = NaN*ones(length(x),Nup);
    %for j=1:Nup
    %    Aup = sqrt(1/trapz(x,conj(phiup(:,j)).*phiup(:,j)));
    %    phiupNor(:,j) = Aup*phiup(:,j);
    %    clear('Aup');
    %end
    Aup = trapz(x,conj(phiup).*phiup,1);
    Aup = sqrt(Aup);
    phiupNor = phiup(:,1:Nup)./Aup(1:Nup);
    phiupNorVirt = phiup(:,Nup+1:end)./Aup(Nup+1:end);
    
    %phidnNor = NaN*ones(length(x),Ndn);
    %for j=1:Ndn
    %    Adn = sqrt(1/trapz(x,conj(phidn(:,j)).*phidn(:,j))); 
    %    phidnNor(:,j) = Adn*phidn(:,j);
    %    clear('Adn');
    %end
    Adn = trapz(x,conj(phidn).*phidn,1);
    Adn = sqrt(Adn);
    phidnNor = phidn(:,1:Ndn)./Adn(1:Ndn);
    phidnNorVirt = phidn(:,Ndn+1:end)./Adn(Ndn+1:end);

    nupNew = KS_densities(phiupNor,Nup);
    ndnNew = KS_densities(phidnNor,Ndn);

    nNew = nupNew + ndnNew;
    
    epsup = eup(1:Nup);
    epsupVirt = eup(Nup+1:end);
    
    epsdn = edn(1:Ndn);
    epsdnVirt = edn(Ndn+1:end);
end

function [outn,outB,outc,outNArray,outrhoArray] = ...
                                    DIIS_updater_n(rhoArray,NArray,B,i,x)
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
    outNArray = NArray;
    outrhoArray = rhoArray;
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