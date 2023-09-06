function [TsAlpha,EextAlpha,EHAlpha,EXCAlpha,EAlpha,DensAlpha,totDens,EpAlpha,EpAlphaPlusOne,DenspAlpha,DenspAlphaPlusOne,...
            phiupNorpAlpha,phidnNorpAlpha,phiupNorpAlphaPlusOne,phidnNorpAlphaPlusOne, ...
            euppAlpha,ednpAlpha,euppAlphaPlusOne,ednpAlphaPlusOne] = fragment_calculator(RCell,vCell,pAlpha,wAlpha,pAlphaPlusOne,x,vp,veeMatrix,nGuess,tolDFT);
        
    % Initialization of fragment properties
    Nfrag = length(vCell);                                               % number of fragments
    TsAlpha = NaN*ones(1,Nfrag);                                 % non-interacting kinetic energy of fragment alpha
    EextAlpha = NaN*ones(1,Nfrag);                              % external energy of fragment alpha
    EHAlpha = NaN*ones(1,Nfrag);                               % Hartree energy of fragment alpha
    EXCAlpha = NaN*ones(1,Nfrag);                              % exchange-correlation (XC) energy of fragment alpha
    TspAlpha = NaN*ones(1,Nfrag);                               
    EextpAlpha = NaN*ones(1,Nfrag);                              
    EHpAlpha = NaN*ones(1,Nfrag);
    EXCpAlpha = NaN*ones(1,Nfrag);
    TspAlphaPlusOne = NaN*ones(1,Nfrag);
    EextpAlphaPlusOne = NaN*ones(1,Nfrag);
    EHpAlphaPlusOne = NaN*ones(1,Nfrag);
    EXCpAlphaPlusOne = NaN*ones(1,Nfrag);
    EAlpha = NaN*ones(1,Nfrag);                                  % total energy of fragment alpha
    EpAlpha = NaN*ones(1,Nfrag);                                  
    EpAlphaPlusOne = NaN*ones(1,Nfrag);                                
    DensAlpha = NaN*ones(length(vp),Nfrag);             % density of fragment alpha
    DenspAlpha = NaN*ones(length(vp),Nfrag);            
    DenspAlphaPlusOne = NaN*ones(length(vp),Nfrag);             
    phiupNorpAlpha = cell(1,Nfrag);                            % spin-up orbitals for ensemble components
    phidnNorpAlpha = cell(1,Nfrag);                            % spin-down orbitals for ensemble components
    euppAlpha = cell(1,Nfrag);                                     % spin-up eigenvalues for ensemble components
    ednpAlpha = cell(1,Nfrag);                                    % spin-down eigenvalues for ensemble components
    phiupNorpAlphaPlusOne = cell(1,Nfrag);                    
    phidnNorpAlphaPlusOne = cell(1,Nfrag);                          
    euppAlphaPlusOne = cell(1,Nfrag);                                 
    ednpAlphaPlusOne = cell(1,Nfrag);                                    
    
    % vp added to vext
    vextMatrix = build_vext_partition(x,RCell,vCell,vp);
    
    % Fragment calculations
    for i = 1:Nfrag
        fprintf('>> \n');
        fprintf('>> fragment_calculator_restricted: calculation for fragment %d start\n',i);
        fprintf('>> p = %d \n', pAlpha(i));
        fprintf('>> p+1 = %d \n', pAlphaPlusOne(i));
        fprintf('>> w = %d \n', wAlpha(i));
        fprintf('>> \n');
        
        % The first term in PPLB eq.
        [Tsup1,Tsdn1,Eext1,EH1,EX1,EC1,Etot1,nup1,ndn1,~,~,~,~,~,~,phiupNor1,phidnNor1,~,~,eup1,edn1,~,~]...
            = dft_1D_restricted_DIIS_constrained(nGuess{i},double(pAlpha(i)),veeMatrix,vextMatrix{i},x,tolDFT);
        EpAlpha(i) = Etot1;
        DenspAlpha(:,i) = nup1 + ndn1;
        phiupNorpAlpha{i} = phiupNor1;
        phidnNorpAlpha{i} = phidnNor1;
        euppAlpha{i} = eup1;
        ednpAlpha{i} = edn1;
        TspAlpha(i) = Tsup1 + Tsdn1;
        EextpAlpha(i) = Eext1;
        EHpAlpha(i) = EH1;
        EXCpAlpha(i) = EX1 + EC1;
        
        % The second term in PPLB eq.
        [Tsup2,Tsdn2,Eext2,EH2,EX2,EC2,Etot2,nup2,ndn2,~,~,~,~,~,~,phiupNor2,phidnNor2,~,~,eup2,edn2,~,~]...
            = dft_1D_restricted_DIIS_constrained(2*nGuess{i},double(pAlphaPlusOne(i)),veeMatrix,vextMatrix{i},x,tolDFT);
        EpAlphaPlusOne(i) = Etot2;
        DenspAlphaPlusOne(:,i) = nup2 + ndn2;
        phiupNorpAlphaPlusOne{i} = phiupNor2;
        phidnNorpAlphaPlusOne{i} = phidnNor2;
        euppAlphaPlusOne{i} = eup2;
        ednpAlphaPlusOne{i} = edn2;
        TspAlphaPlusOne(i) = Tsup2 + Tsdn2;
        EextpAlphaPlusOne(i) = Eext2;
        EHpAlphaPlusOne(i) = EH2;
        EXCpAlphaPlusOne(i) = EX2 + EC2;
        
        % Fragment properties via ensemble
        TsAlpha(i) = (1-wAlpha(i))*TspAlpha(i) + wAlpha(i)*TspAlphaPlusOne(i);
        EextAlpha(i) = (1-wAlpha(i))*EextpAlpha(i) + wAlpha(i)*EextpAlphaPlusOne(i);
        EHAlpha(i) = (1-wAlpha(i))*EHpAlpha(i) + wAlpha(i)*EHpAlphaPlusOne(i);
        EXCAlpha(i) = (1-wAlpha(i))*EXCpAlpha(i) + wAlpha(i)*EXCpAlphaPlusOne(i);
        EAlpha(i) = (1-wAlpha(i))*EpAlpha(i) + wAlpha(i)*EpAlphaPlusOne(i);
        DensAlpha(:,i) = (1-wAlpha(i))*DenspAlpha(:,i)  + wAlpha(i)*DenspAlphaPlusOne(:,i);
    end
    totDens = sum(DensAlpha,2);
end
    
    
