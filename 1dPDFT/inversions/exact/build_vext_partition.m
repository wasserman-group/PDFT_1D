function out = build_vext_partition(x,RCell,vCell,vemb)
% Function returns external potential for fragments of a hydrogen chain.
% x -> is the grid vector. 
% RCell is the cell which elements are the 1d arrays with coordinates of the atoms.
    Nfrag = length(vCell);
    
    if (length(RCell) ~= Nfrag)
        fprintf('>> build_vext: Warning! Wrong number of atoms \n');
        return
    end
    
    vextMatrix = cell(1,Nfrag);

    for i = 1:Nfrag
        vextMatrix{i} = vCell{i}(x,RCell{i}) + vemb;
    end
    
    out = vextMatrix;
end