function out = build_hamiltonian(N,Nele,h,vextMatrix,veeMatrix)
    if (Nele)
        %K = build_Kin(N,Nele,h);
        %K = build_Kin4(N,Nele,h);
        K = build_Kin6(N,Nele,h);
        %
        %arrayV = build_V(Nele,N,veeMatrix,vextMatrix);
        arrayV = fast_build_V(Nele,N,veeMatrix,vextMatrix);
        %
	H = K + spdiags(arrayV,0,N^Nele,N^Nele);
	fprintf('>> Number of elements in H: %d \n',nnz(H));
	fprintf('>> Size of H:%d \n',numel(H));
	% 
	out = H;
    else
        out = 0.0;
    end
end

function out = fast_build_V(Nele,N,veeMatrix,vextMatrix)
    tic
    g = 1:1:N;
    C = cell(Nele,1);
    [C{:}] = ndgrid(g);
    indexext = cellfun(@(g){g(:)}, C);
    indexext = [indexext{:}];
    clear('C');
    if (Nele > 1)
        % "small" double loop here -- I think it's the only way to do it
        k = 1;
        indexee = NaN*ones(N^Nele,nchoosek(Nele,2));
        for i = 1:Nele
            for j = (i+1):Nele
                indexee(:,k) = abs(indexext(:,i) - indexext(:,j));
                k = k + 1;
            end
        end
        %
        vee = sum(veeMatrix(indexee+1),2);
    end
    %
    vext = sum(vextMatrix(indexext),2);
    timeV = toc;
    fprintf('>> time to build V: %2.2e sec\n',timeV);
    %
    if (Nele > 1)
        out = vee + vext;
    else
        out = vext;
    end
end

function out = build_V(Nele,N,veeMatrix,vextMatrix)
%Same as fast_build_V but works very slowly. The function is more readable
%and testable
    tic;
    iMatrix = NaN*ones(N^Nele,Nele);
    vDiagonal = NaN*ones(N^Nele,1);
%     fullVext = NaN*ones(N^Nele,1);
%     fullVee = NaN*ones(N^Nele,1);
    for k = 1:N^Nele   
        temp = k - 1;
        for i = 1:Nele
            iMatrix(k,i) = rem(fix(temp/(N^(i-1))),N) + 1;
        end
        
        clear('temp');
        vee = 0;
        vext = 0;
        for i = 1:Nele
            for j = (i+1):Nele
                vee = vee + veeMatrix(abs(iMatrix(k,j) - iMatrix(k,i)) + 1);
            end
            vext = vext + vextMatrix(iMatrix(k,i));
        end

        fullVext(k,1) = vext;
        fullVee(k,1) =  vee;
        vDiagonal(k,1) = vext + vee;
        clear('vext');
        clear('vee');        
    end
    
    clear('iMatrix');
    clear('veeMatrix');
    clear('vextMatrix');
    timeV = toc;
    fprintf('>> time to build V: %2.2e sec\n',timeV);
    out = vDiagonal;
end
