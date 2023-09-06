addpath(genpath('./potentials'));
addpath(genpath('./inversions'));
addpath(genpath('./dft'));
addpath(genpath('./temp_mat_files'));
% Set up:
boxL = 25;
N = 251;
tolvp = 1e-6;
tolDFT = 1e-6;
Nele = int16(2.0); % make sure it is int format when it is integer
vCell = cell(1,Nele);

for i = 1:Nele
    vCell{i} = @v_ext1X2;
end

Nfrag = length(vCell);
B0 = 5.0;
x = linspace(0,boxL,N);
x= x';
h = x(2) - x(1);
veeHandle = @v_ee_DFT;
veeMatrix = build_vee(veeHandle,N,h,1);
RCell = HChain_coord(B0,boxL,Nfrag);

% Print details:
fprintf('>> current functions are: \n');
fprintf('>> %s \n',func2str(vCell{1}));
fprintf('>> %s \n',func2str(veeHandle));
fprintf('>> current number of grid points : %d\n',N);
fprintf('>> current box length : %d\n',boxL);
fprintf('>> current internuclear separation : %d\n',B0);
fprintf('>> current number of electrons : %d\n', double(Nele));
fprintf('>> vp convergence tolerance: %e\n',tolvp);
fprintf('>> S-DFT convergence tolerance: %e\n',tolDFT);

% Download full-molecular density:
load('H2_partition_R5.mat');

% Sort initial guess of density
vextMatrix = build_vext_partition(x,RCell,vCell,zeros(size(x)));
vextMatrix = cell2mat(vextMatrix);
w = zeros(N,Nfrag);
for i = 1:Nfrag
    w(:,i) = vextMatrix(:,i)./sum(vextMatrix,2);
end

nGuess = cell(1,Nfrag);
for i = 1:Nfrag
    nGuess{i} = densMol.*w(:,i);
end

NAlpha = [1,1];

[EAlpha,EfragAlpha,DensAlpha,totDens,totDensup,totDensdn,vp,vpH,vpext,vpXC,vpkin, ...
    Etot,totTs,totEext,totEH,totEXC,Ep_OA,Ep,Epkin,Epext,EpH,EpXC,vpContr,optimality, ...
    DenspAlpha,DenspAlphaPlusOne,EpAlpha,EpAlphaPlusOne,muAlpha,S,dSdn_Alpha,n_temp,vp_temp] = ...
    pdft_1d_restricted_GOA(NAlpha,x,N,h,RCell,vCell,veeMatrix,tolvp,tolDFT,nGuess);


save('output_H2_R5_GOA');
