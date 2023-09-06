function [Ts,vXC_inv,vXC_fun] = wuyang_vXC(n0,n0up,n0dn,x,h,N,Nele,vextMatrix)
% This function uses matlab's built in optimization function
% fminunc to maximize Ts with respect to the potential

options = optimoptions(@fminunc,'Algorithm','quasi-newton',...
    'SpecifyObjectiveGradient',true,...
    'CheckGradients',false,...
    'OptimalityTolerance',1e-7, ...
    'MaxIterations',1500,...
    'Display','iter');

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

vxup = v_x(n0,n0up,n0dn,'up');
vxdn = v_x(n0,n0up,n0dn,'dn');
vcup = v_c(n0,n0up,n0dn,'up');
vcdn = v_c(n0,n0up,n0dn,'dn');
vc = (1/2)*(vcup + vcdn);
vx = (1/2)*(vxup + vxdn);
vXC_fun = vx + vc;
veeHandle = @v_ee_DFT;
veeMatrix = build_vee(veeHandle,N,h,1);
vH = v_H(n0,x,veeMatrix);
vext = v_ext(vextMatrix,zeros(size(x)));

v0 = zeros(size(x));

[v_opt,Gval,exitflag,output,gradient]= fminunc(@G,v0,options);
vs = vext + vH + +vXC_fun + v_opt;
vXC_inv = vs - vext - vH;

    function [G,gradient] = G(v_opt)
        % This function calculates G its gradient and its hessian
        % for a given potential.

        v_guide = vH + vXC_fun;
        vs = vext + v_guide + v_opt;

        % orbitals
        H = build_KS_hamiltonian(N,h,vs);
        [phi,e] = eig(full(H));
        [e,ind] = sort(diag(e));
        e = diag(e);
        phi = phi(:,ind);
        clear('ind'); 
        A = sqrt(trapz(x,conj(phi).*phi,1));
        phiNor = phi(:,1:Norb)./A(1:Norb);
        phiupNor = phiNor(:,1:Nup);
        phidnNor = phiNor(:,1:Ndn);

        % density
        nup = KS_densities(phiNor,Nup);
        ndn = KS_densities(phiNor,Ndn);
        n = nup + ndn;
        
        % Kinetic energy
        Tsup = 0.0;
        Tsdn = 0.0;
        K = build_Kin6(N,1,h);
        if (Nup)
            for j = 1:Nup
                Tsup = Tsup + trapz( x,conj(phiupNor(:,j)).*(K*phiupNor(:,j)) );
            end
        end
        
        if (Ndn)
            for j = 1:Ndn
                Tsdn = Tsdn + trapz( x,conj(phidnNor(:,j)).*(K*phidnNor(:,j)) );
            end
        end
        
        Ts = Tsup + Tsdn;


        % difference between fragment sum and target.
        d_nf = (n-n0);

        % calculate objective function and its gradient
        G = -(Ts + sum((vs).*(d_nf))*h);
        gradient = -(d_nf)*h;

    end
end