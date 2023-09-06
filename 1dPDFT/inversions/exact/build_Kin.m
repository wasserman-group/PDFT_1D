function out = build_Kin(N,Nele,h)
%Builds kinetic energy operator on a real grid of N points for 1D systeme of
%Nele electrons ~O(h^2).
    tic;
    K = spalloc(N^Nele,N^Nele,(2*Nele + 1)*N^Nele);
    for i = 1:Nele
        szDiag = N^Nele - N^(i-1) + 1 - 1;
        szPieces = (N - 1)*N^(i-1);
        numPieces = N^(Nele - 1 - (i-1));

        currDiag = zeros(1,N^Nele);
        currDiag(1:szDiag) = 1;
        szBlank = (szDiag - numPieces*szPieces)/(numPieces - 1);

        if (~isnan(szBlank))
            for j = 1:szBlank
                currDiag(szPieces+j:szPieces+szBlank:end) = 0;
            end
        end

        K = K + spdiags(currDiag',-N^(i-1),N^Nele,N^Nele);
    end

    K = K + K' + spdiags(-2*Nele*ones(N^Nele,1),0,N^Nele,N^Nele);
    out = (-1/2)*(1/(h^2))*K;
    timeK = toc;
    fprintf('>> Kinetic part build time: %2.2e sec\n',timeK);
end