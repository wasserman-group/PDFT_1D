function out = HChain_coord(B0,boxL,Nfrag)

R = cell(1,Nfrag);

for i = 1:Nfrag
    R{i} = boxL/2.0 - B0/2 - (Nfrag/2 - i).*B0;
end

out = R;

end