function [outQpAlpha,outQpAlphaPlusOne] = calc_local_Q(Nfrag,N,totDens,DenspAlpha,DenspAlphaPlusOne)
QpAlpha = NaN*ones(N,Nfrag);
QpAlphaPlusOne = NaN*ones(N,Nfrag);

for i = 1:Nfrag
    QpAlpha(:,i) = DenspAlpha(:,i)./totDens;
    QpAlphaPlusOne(:,i) = DenspAlphaPlusOne(:,i)./totDens;
end

outQpAlpha = QpAlpha;
outQpAlphaPlusOne = QpAlphaPlusOne;