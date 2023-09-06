function out = overlap(c,x,DensAlpha)

temp = 0;
Nfrag = size(DensAlpha,2);
for i = 1:(Nfrag-1)
    temp = temp + trapz(x,sqrt(DensAlpha(:,i).*DensAlpha(:,i+1)));
end
temp = temp/(Nfrag-1);
out = erf(c*temp);

end