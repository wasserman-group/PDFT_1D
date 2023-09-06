function out = build_Kin6(N,Nele,h)
%Builds kinetic energy operator on a real grid of N points for 1D systeme of
%Nele electrons ~O(h^6).
    K = spalloc(N^Nele,N^Nele,(2*Nele + 1)*N^Nele);
    if (Nele == 1)
        diag1 = (3/2)*ones(N-1,1);
        diag2 = (-3/20)*ones(N-1,1);
        diag3 = (1/90)*ones(N-2,1);
        K = spdiags(diag1,-1,N,N);
        K = K + spdiags(diag2,-2,N,N);
        K = K + spdiags(diag3,-3,N,N);
        K = K + K';
        %
        K = K + spdiags(-(49/18)*(ones(N,1)),0,N,N);
    elseif (Nele == 2)
        diag1 = (3/2)*ones(N-1,1);
        diag2 = (-3/20)*ones(N-1,1);
        diag3 = (1/90)*ones(N-2,1);
        Kgen = spdiags(diag1,-1,N,N);
        Kgen = Kgen + spdiags(diag2,-2,N,N);
        Kgen = Kgen + spdiags(diag3,-3,N,N);
        Kgen = Kgen + Kgen';
        %
        Igen = spdiags(ones(N,1),0,N,N);
        %
        K = kron(Igen,Kgen) + kron(Kgen,Igen);
        K = K + spdiags(-(49/9)*(ones(N^2,1)),0,N^2,N^2);
    elseif (Nele == 3)
        % if it's wrong, blame Laura
        diag1 = (3/2)*ones(N-1,1);
        diag2 = (-3/20)*ones(N-1,1);
        diag3 = (1/90)*ones(N-2,1);
        Kgen = spdiags(diag1,-1,N,N);
        Kgen = Kgen + spdiags(diag2,-2,N,N);
        Kgen = Kgen + spdiags(diag3,-3,N,N);
        Kgen = Kgen + Kgen';
        %
        Igen = spdiags(ones(N,1),0,N,N);
        %
        K = kron(Igen,kron(Igen,Kgen)) + ...
            kron(Igen,kron(Kgen,Igen)) + ...
            kron(Kgen,kron(Igen,Igen));
        K = K + spdiags(-(49/6)*(ones(N^3,1)),0,N^3,N^3);
    elseif (Nele == 0)
        K = 0;
    else
        fprintf('>> Kinetic to O(h^6): cannot build K\n');
    end
    %
    out = (-1/2)*(1/(h^2))*K;
end