function out = build_vext(x,RCell,vCell,vOther)
%Function returns external potential of the system. x -> is the grid vector
%RCell is the cell which elements are the 1d arrays with coordinates of the
%atoms and vCell is the cell of the external potential function handles  
    Nfrag = length(vCell);
    %
    if (length(RCell) ~= Nfrag)
        fprintf('>> build_vext: Warning! Wrong number of atoms \n');
        return
    end
    %
    vextMatrix = zeros(size(x));
    %
    for i = 1:Nfrag
        vextMatrix = vextMatrix + vCell{i}(x,RCell{i});
    end
    %
    %lift = vextMatrix(1)/2 + vextMatrix(end)/2;
    lift = 0;
    %
    out = vextMatrix + vOther + (-lift);
end