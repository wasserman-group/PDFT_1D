function [vpkin,vpext,vpH,vpXC] = vp_components(Nfrag,x,N,h,RCell,vCell,veeMatrix,vp,totDens,DenspAlpha,DenspAlphaPlusOne,wAlpha)

[QpAlpha,QpAlphaPlusOne] = calc_local_Q(Nfrag,N,totDens,DenspAlpha,DenspAlphaPlusOne);

vpext = vp_ext(QpAlpha,QpAlphaPlusOne,x,RCell,vCell,wAlpha);

vpH = vp_H(QpAlpha,QpAlphaPlusOne, wAlpha, ...
                totDens,DenspAlpha,DenspAlphaPlusOne,veeMatrix, ...
                N,h,Nfrag);

vpXC = vp_XC(QpAlpha,QpAlphaPlusOne,...
    totDens,DenspAlpha,DenspAlphaPlusOne,wAlpha,Nfrag);

vpkin = vp - (vpext + vpH + vpXC);
