function [outE,outB] = dimer_geom(v_atomic,vCell,veeMatrix,Nele,B0,testScan,x)
% Function that calculates optimal bond length for a 1D dimer
    %% Define some stuff
    box_l = x(end);
    h = x(2) - x(1);
    N = length(x);
    MaxIter = 10;
    i = 1;
    optimality = 1;
    convTol = 1e-8;
    numEig = 3;
    numReq = 2;
    %% Test Scan
    if (testScan)
        BT = NaN*ones(1,10);
        ET = NaN*ones(1,10);
        for l=1:10
            BT(l) = l*1.0;
            %[ET(l),~] = get_E(BT(l),x,numEig,symm,fun_v_ext,fun_v_ee,fun_v_Atomic);
            [R1,R2] = dimer_coord(BT(l),box_l);
            RCell{1} = R1;
            RCell{2} = R2;
            ET(l) = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
            clear('RCell');
            
            fprintf('***SCAN***\n');
            fprintf('>>Geometry opt: scan run\n');
            fprintf('>>Geometry opt: scan energy %.4f\n',ET(l));
            fprintf('>>Geometry opt: scan bond length %.2f\n',BT(l));
            fprintf('*********\n');
        end
        plot(BT,ET)
        fprintf('\n');
        inB = input('>>Geometry opt: enter new bond guess or click enter\n');
        if (inB)
            B0 = inB;
        end
        fprintf('\n');
        %% Calculate E at the initial guess
        %[E0,~] = get_E(B0,x,numEig,symm,fun_v_ext,fun_v_ee,fun_v_Atomic);
        [R10,R20] = dimer_coord(B0,box_l);
        RCell{1} = R10;
        RCell{2} = R20;
        E0 = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
        clear('RCell')
        
        fprintf('\n');
        fprintf('>>Geometry opt: pre-cycle run\n');
        fprintf('>>Geometry opt: current energy %.4f\n',E0);
        fprintf('>>Geometry opt: current bond length %.2f\n',B0);
        fprintf('>>Geometry opt: current optimality %5.2e\n',optimality);
        fprintf('\n');
    else
        [R10,R20] = dimer_coord(B0,box_l);
        RCell{1} = R10;
        RCell{2} = R20;
        E0 = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
        clear('RCell');
    end
    %% Optimization  
    scale = 10;   
    while (i <= MaxIter) && (optimality > convTol)
        step = 1/(i*scale);

        BR0 = B0;
        BL0 = B0;
        ER0 = E0;
        EL0 = E0;
        
        BRT = BR0 + step;
        BLT = BL0 - step;
        
        %[ERT,~] = get_E(BRT,x,numEig,symm,fun_v_ext,fun_v_ee,fun_v_Atomic);
        [R1RT,R2RT] = dimer_coord(BRT,box_l);
        RCell{1} = R1RT;
        RCell{2} = R2RT;
        ERT = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
        clear('RCell')
        
        %[ELT,~] = get_E(BLT,x,numEig,symm,fun_v_ext,fun_v_ee,fun_v_Atomic);
        [R1LT,R2LT] = dimer_coord(BLT,box_l);
        RCell{1} = R1LT;
        RCell{2} = R2LT;
        ELT = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
        clear('RCell')
        
        if (ELT < E0) && (ERT > E0)
            EL0 = ELT;
            BL0 = BLT;
            while ELT <= EL0
                EL0 = ELT;
                BL0 = BLT;
                BLT = BLT - step;
                %[ELT,~] = get_E(BLT,x,numEig,symm,fun_v_ext,fun_v_ee,fun_v_Atomic);
                [R1LT,R2LT] = dimer_coord(BLT,box_l);
                RCell{1} = R1LT;
                RCell{2} = R2LT;
                ELT = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
                clear('RCell')
                
                fprintf('>>Geometry opt: walking to the left\n');
                fprintf('>>Geometry opt: iteration %d current BOND LENGTH = %.4f ENERGY = %.4f\n',i,BLT,ELT);
                BMIN = fitPar([B0,BL0,BLT],[E0,EL0,ELT]);
                %EMIN = get_E(BMIN,x,numEig,symm,fun_v_ext,fun_v_ee,fun_v_Atomic);
                [R1MIN,R2MIN] = dimer_coord(BMIN,box_l);
                RCell{1} = R1MIN;
                RCell{2} = R2MIN;
                EMIN = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
                clear('RCell')      
            end
        elseif (ELT > E0) && (ERT < E0)
            ER0 = ERT;
            BR0 = BRT;
            while ERT <= ER0
                ER0 = ERT;
                BR0 = BRT;
                BRT = BRT + step;
                %[ERT,~] = get_E(BRT,x,numEig,symm,fun_v_ext,fun_v_ee,fun_v_Atomic);
                [R1RT,R2RT] = dimer_coord(BRT,box_l);
                RCell{1} = R1RT;
                RCell{2} = R2RT;
                ERT = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
                clear('RCell')
                
                fprintf('>>Geometry opt: walking to the right\n');
                fprintf('>>Geometry opt: iteration %d current BOND LENGTH = %.4f ENERGY = %.4f\n',i,BRT,ERT);
                BMIN = fitPar([B0,BR0,BRT],[E0,ER0,ERT]);
                %EMIN = get_E(BMIN,x,numEig,symm,fun_v_ext,fun_v_ee,fun_v_Atomic);
                [R1MIN,R2MIN] = dimer_coord(BMIN,box_l);
                RCell{1} = R1MIN;
                RCell{2} = R2MIN;
                EMIN = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
                clear('RCell')
            end
        elseif (ELT > E0) && (ERT > E0)
            % fit parabola using E0, EL0, and ER0
            BMIN = fitPar([B0,BRT,BLT],[E0,ERT,ELT]);
            %EMIN = get_E(BMIN,x,numEig,symm,fun_v_ext,fun_v_ee,fun_v_Atomic);
            [R1MIN,R2MIN] = dimer_coord(BMIN,box_l);
            RCell{1} = R1MIN;
            RCell{2} = R2MIN;
            EMIN = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq);
            clear('RCell')
        else
            ;% come up with simething...
            fprintf('\n');
            fprintf('>>>>Geometry opt: Multiple minima\n');
            fprintf('>>>>Geometry opt: Try different initial guess or adjust step size\n');
            fprintf('\n');
            outE = E0;
            outB = B0;
            return
        end

        optimality = abs(E0 - EMIN);
        
        fprintf('\n');
        fprintf('>>Geometry opt: %d iteration complete\n',i);
        if (E0 - EMIN) < 0
            fprintf('>>Geometry opt: Energy went up after this iteration\n');
        end
        fprintf('>>Geometry opt: current BOND LENGTH = %.4f ENERGY = %.4f\n',BMIN,EMIN);
        fprintf('>>Geometry opt: previous BOND LENGTH = %.4f ENERGY = %.4f\n',B0,E0);
        fprintf('>>Geometry opt: current optimality %5.2e\n',optimality);
%         fprintf('>>Geometry opt: click to contine\n');
        fprintf('\n');
%         pause;
        
        i = i + 1;
        if (E0 - EMIN) > 0
            E0 = EMIN;
            B0 = BMIN;
        end
    end
    
    outE = E0;
    outB = B0;
    
end

function out = fitPar(x,y)
%returns x_min        
        x1 = x(1);
        y1 = y(1);
        x2 = x(2);
        y2 = y(2);
        x3 = x(3);
        y3 = y(3);
        
        A = NaN*ones(3,3);
        d = NaN*ones(3,1);
        A(1,1) = x2^2;
        A(1,2) = x2;
        A(1,3) = 1;
        
        A(2,1) = x3^2;
        A(2,2) = x3;
        A(2,3) = 1;
        
        A(3,1) = x1^2;
        A(3,2) = x1;
        A(3,3) = 1;
        
        d(1,1) = y2;
        d(2,1) = y3;
        d(3,1) = y1;
        
        c = A\d;
        
        xmin = -c(2)/(2*c(1));
        out = xmin;
end

function out = get_E(v_atomic,veeMatrix,RCell,vCell,x,h,N,Nele,numEig,numReq)
    vextMatrix = build_vext(x,RCell,vCell,zeros(size(x)));
    H = build_hamiltonian(N,Nele,h,vextMatrix,veeMatrix);
    [~,eigenE,~] = calc_Eig(H,N,Nele,numEig,numReq);
    out = eigenE(numReq) + v_atomic(abs(RCell{1} - RCell{2}));
end