function out = vp_H(QpAlpha,QpAlphaPlusOne, wAlpha, ...
                totDens,DenspAlpha,DenspAlphaPlusOne,veeMatrix, ...
                N,h,Nfrag)
    %
    vpHartree = zeros(size(totDens));
    %
    for i =1:Nfrag    
        vpHartreepAlpha = zeros(size(totDens));
        vpHartreepAlphaPlusOne = zeros(size(totDens));

        for j = 1:N
            int1 = 0;
            int2 = 0;

                for k = 1:N
                    int1 = int1 + h*(totDens(k) - DenspAlpha(k,i)).*veeMatrix(abs(j-k)+1);
                    int2 = int2 + h*(totDens(k) - DenspAlphaPlusOne(k,i)).*veeMatrix(abs(j-k)+1);
                end

            vpHartreepAlpha(j)        = vpHartreepAlpha(j)        + (1-wAlpha(i))*int1;
            vpHartreepAlphaPlusOne(j) = vpHartreepAlphaPlusOne(j) + wAlpha(i)*int2;
        end
        %
        vpHartree = vpHartree + vpHartreepAlpha.*QpAlpha(:,i) + vpHartreepAlphaPlusOne.*QpAlphaPlusOne(:,i);
        clear('vpHartreepAlpha');
        clear('vpHartreepAlphaPlusOne');
    end    
    %
    out = vpHartree;
end