function out = build_nGuess(N,n0,R)
% This function generates initial guess of densities for fragments in PDFT.
% N: number of grid points
% n0: isolated fragment density
% R: nuclear coordinates

Nc = 10*R + 1;
N1 = int16(Nc - (length(n0) + 1)/2);
N2 = int16(Nc + (length(n0) - 1)/2);

temp = zeros(N,1);

if (1 <= N1) && (N1 < N2) && (N2 <= N)
    
    for i = 1:N1
        temp(i) = 1e-13;
    end
    
    for i = (N1+1):N2
        temp(i) = n0(i-N1);
    end
    
    for i = (N2+1):N
        temp(i) = 1e-13;
    end
    
else
    fprintf('Incorrect fragment density or nuclear coordinates \n');
end

out = temp;
end
      