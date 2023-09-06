function out = v_H(n,x,veeMatrix)
% function v_H(n,h,veeMatrix)
% n <- density
% h <- dx increment
% veeMatrix <- array generated for given 1D interaction and grid.
% veeMatrix(index) should return interaction potential between 2 electrons
% located (index-1) grid points away from each other. For example:
% veeMatrix(1) = v_ee(|xi-xj| = 0); veeMatrix(7) = v_ee(|xi-xj| = 6dx).
    temp = zeros(size(n));
    N = length(x);
    h = x(2) - x(1);
    currentIndex = 1:1:(N);
    % better "MatLab" way of doing this ??
%     h = x(2) - x(1);
    for i = 1:length(temp)
        currentIndex = currentIndex - 1;
%         temp(i) = trapz( x,n.*veeMatrix(int16(abs(currentIndex)+1)) );
        temp(i) = sum(n.*veeMatrix(int16(abs(currentIndex)+1)))*h;
    end
    
    out = temp;
end
