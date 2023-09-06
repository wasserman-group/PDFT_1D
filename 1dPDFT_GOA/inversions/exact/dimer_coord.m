function [out1,out2] = dimer_coord(B,box_l)
    R1 = 0.0 + (box_l - B)/2;
    R2 = box_l - (box_l - B)/2;
    out1 = R1;
    out2 = R2;
end