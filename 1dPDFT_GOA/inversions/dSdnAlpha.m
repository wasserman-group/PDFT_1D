function out = dSdnAlpha(c,x,DensAlpha)

p = 0;
temp = 0;
Nfrag = size(DensAlpha,2);
for i = 1:(Nfrag-1)
    p = p + trapz(x,sqrt(DensAlpha(:,i).*DensAlpha(:,i+1)));
end
p = p/(Nfrag-1);

for i = 1:(Nfrag-1)
    temp = temp + sqrt(DensAlpha(:,i).*DensAlpha(:,i+1));
end
temp = temp/(Nfrag-1);

out = 2*c/sqrt(pi)*temp*exp(-(c*p)^2);

end