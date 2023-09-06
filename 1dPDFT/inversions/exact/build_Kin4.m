function out = build_Kin4(N,Nele,h)
%Builds kinetic energy operator on a real grid of N points for 1D systeme of
%Nele electrons ~O(h^4).
    tic;
    K = spalloc(N^Nele,N^Nele,(2*Nele + 1)*N^Nele);
    if (Nele == 1)
        diag1 = 16*ones(N-1,1);
        diag2 = -1*ones(N-1,1);
        K = spdiags(diag1,-1,N,N);
        K = K + spdiags(diag2,-2,N,N);
        K = K + K';
        
        K = K + spdiags(-30*(ones(N,1)),0,N,N);
    elseif (Nele == 2)
        diag1 = 16*ones(N-1,1);
        diag2 = -1*ones(N-1,1);
        Kgen = spdiags(diag1,-1,N,N);
        Kgen = Kgen + spdiags(diag2,-2,N,N);
        
        
        Kgen = Kgen + Kgen';
        Igen = spdiags(ones(N,1),0,N,N);
        
        K = kron(Igen,Kgen) + kron(Kgen,Igen);
        K = K + spdiags(-60*(ones(N^2,1)),0,N^2,N^2);
    elseif (Nele == 0)
        K = 0;
    else
        fprintf('>> Kinetic to O(h^4): cannot build K\n');
    end
    out = (-1/2)*(1/(12*h^2))*K;
    timeK = toc;
    fprintf('>> Kinetic to O(h^4) part build time: %2.2e sec\n',timeK);
end