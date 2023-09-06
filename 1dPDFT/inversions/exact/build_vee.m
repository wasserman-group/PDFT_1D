function out = build_vee(v_ee,N,h,lambda)
%Function returns electron-electron potential matrix for 2 electrons. The
%index of the matrix is the distance (in h = dx units) between two
%electrons. v_ee -> function handle, h -> dx, N -> number of grid points
    veeMatrix = v_ee(N,h,lambda);
    %
    out = veeMatrix;
end