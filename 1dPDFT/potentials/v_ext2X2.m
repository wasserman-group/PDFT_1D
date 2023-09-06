function out = v_ext2X2(x,R2)
    b = 1.0;
    out = -1.00./sqrt(b + (x - R2).^2);
end