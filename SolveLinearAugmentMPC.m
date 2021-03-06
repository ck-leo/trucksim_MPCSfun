function [delta_control] = SolveLinearAugmentMPC(matrix_a,matrix_b,matrix_c,pre_matrix_c,matrix_q,matrix_r,...
    matrix_lower,matrix_upper,matrix_state,pre_matrix_state,reference,pre_reference,control)
% [F,Phi,Phid,A_e, B_e,C_e] = mpcgain(Ackesi,Bckesi,Bdkesi,Cckesi,Nc,Np)
% obtain F,Phi,Phid and A_e,B_e,C_e, from parameter
% Ackesi,Bckesi,Bdkesi,Cckesi,Np,Nc of the augmented discretization model
% the discretization model:
% kesi(k+1) = Ackesi*kesi(k) + Bckesi*u(k) + Bdkesi*d(k)
%     y(k)  = Cckesi*kesi(k)
% augmented model:
% [delta_kesi(k+1)] [Ackesi         0ckesi'] [delta_kesi(k)]
% [         y(k+1)]=[Cckesi*Ackesi     I   ]*[         y(k)]+
%                          A_e(A)
%                   [Bckesi       ]           [Bdkesi       ]
%                   [Cckesi*Bckesi]*deltau(k)+[Cckesi*Bdkesi]*d(k)
%                          B_e(B)                  Bd_e(Bd)
%                              [delta_kesi(k)]
% [         y(k)]=[0ckesi   I]*[         y(k)]
%                     C_e(C)
%   [CA   ]
%   [CA^2 ]
% F=[...  ]
%   [CA^Np]
% narginchk(5,7);
% if (size(matrix_a,1) ~= size(matrix_a,2) ||
%       size(matrix_b,1) ~= size(matrix_a,1) ||
%       size(matrix_lower,1) ~= size(matrix_upper,1))
%     error('One or more matrices have incompatible dimensions. Aborting.');
% end
[nkesi,nu] = size(matrix_b);                             % nkesi:状态变量个数；nu:控制输入参数个数；
ny = nkesi;
nd = size(matrix_c,1);
Ackesi = matrix_a;
Bckesi = matrix_b;
Bdkesi = eye(nd);
Cckesi = eye(nkesi);
Nc = size(control,1)/nu;
dis = matrix_c;
Np = size(reference,1)/nkesi;
% [ny,nkesi] = size(Cckesi);                             % ny:number of output; nkesi:number of state;
% [nkesi,nu] = size(Bckesi);                             % nkesi:number of state; nu:number of control;
% [nkesi,nd] = size(Bdkesi);                             % nkesi:number of state;nd:number of disturbance;
%% calculate A_e,B_e,Bd_e,C_e
% A_e
A_e = eye(nkesi+ny,nkesi+ny);
A_e(1:nkesi,1:nkesi) = Ackesi;
A_e(nkesi+1:nkesi+ny,1:nkesi) = Cckesi*Ackesi;
% B_e
B_e = zeros(nkesi+ny,nu);
B_e(1:nkesi,:) = Bckesi;
B_e(nkesi+1:nkesi+ny,:)=Cckesi*Bckesi;
% Bd_e
Bd_e=zeros(nkesi+ny,nd);
Bd_e(1:nkesi,:)=Bdkesi;
Bd_e(nkesi+1:nkesi+ny,:)=Cckesi*Bdkesi;
% C_e
C_e = zeros(ny,nkesi+ny);
C_e(:,nkesi+1:nkesi+ny) = eye(ny,ny);
%% calculate F,Phi,Phid
% F
F = zeros(ny*Np,nkesi+ny);
F(1:ny,:) = C_e*A_e;
for i=2:Np
    F((i-1)*ny+1:i*ny,:) = F((i-2)*ny+1:(i-1)*ny,:)*A_e;
end
% Phi
Phi = zeros(ny*Np,nu*Nc);
Phi(:,1:nu) = [C_e;F(1:(Np-1)*ny,:)]*B_e;
for i = 2:Nc
    Phi(:,(i-1)*nu+1:i*nu) = [zeros((i-1)*ny,nu);Phi(1:(Np-i+1)*ny,1:nu)];
end
% Phid
Phid = [C_e;F(1:(Np-1)*ny,:)]*Bd_e;
%% initialize Yref,matrix_q,matrix_r,etc.
Yref = [reference-pre_reference;reference];
Xf = [matrix_state-pre_matrix_state;matrix_state];
Delta_dis = matrix_c-pre_matrix_c;
%% calculate Q,R
% Q
Q=zeros(Np*ny,Np*ny);
for j=1:Np
    Q((j-1)*ny+1:j*ny,(j-1)*ny+1:j*ny)=matrix_q;
end
% R
R=eye(Nc*nu,Nc*nu);
for j=1:Nc
    R((j-1)*nu+1:j*nu,(j-1)*nu+1:j*nu)=matrix_r;
end
% umin,umax,u0
umin = matrix_lower;
umax = matrix_upper;
u0 = control(1);
%% calculate QP parameter
[H,G,A_cons,b] = ToQP(F,Phi,Phid,Q,R,nu,Nc,Delta_dis,Xf,Yref,umin,umax,u0);
% caluculate QP problem
delta_control = quadprog(H,G,A_cons,b);
end

