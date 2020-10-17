function control = SolveLinearMPC(matrix_a,matrix_b,matrix_c,matrix_q,matrix_r,matrix_lower,matrix_upper,matrix_initial_state,reference,eps,max_iter,control)
% [F,Phi,Phid,A_e, B_e,C_e] = mpcgain(Ackesi,Bckesi,Bdkesi,Cckesi,Nc,Np)
% 根据离散状态方程（非增量模型）中的Ackesi,Bckesi,Bdkesi,Cckesi以及预测周期Np和控制周期Nc，
% 获取F、Phi、Phid以及增量模型中的A、B、C
% the discretization model:
% kesi(k+1) = Ackesi*kesi(k) + Bckesi*u(k) + Bdkesi*d(k)
%     y(k)  = Cckesi*kesi(k)
% augmented model:
% [delta_kesi(k+1)] [Ackesi         0ckesi'] [delta_kesi(k)]
% [         y(k+1)]=[Cckesi*Ackesi     I   ]*[         y(k)]+
%                          A_e（A）
%                   [Bckesi       ]           [Bdkesi       ]
%                   [Cckesi*Bckesi]*deltau(k)+[Cckesi*Bdkesi]*d(k)
%                          B_e（B）                  Bd_e（Bd）
%                              [delta_kesi(k)]
% [         y(k)]=[0ckesi   I]*[         y(k)]
%                     C_e(C)
%   [CA   ]
%   [CA^2 ]
% F=[...  ]
%   [CA^Np]
% narginchk(5,7);
if (size(matrix_a,1) ~= size(matrix_a,2) ||
      size(matrix_b,1) ~= size(matrix_a,1) ||
      size(matrix_lower,1) ~= size(matrix_upper,1))
    error('One or more matrices have incompatible dimensions. Aborting.');
end
[nkesi,nu] = size(matrix_b);                             % nkesi:状态变量个数；nu:控制输入参数个数；
ny = nkesi;
nd = size(matrix_c,2);
Ackesi = matrix_a;
Bckesi = matrix_b;
Bdkesi = matrix_c;
Cckesi = eye(nkesi);
Nc = size(control,1)/nu;
Np = size(reference,1)/nkesi;
% [ny,nkesi] = size(Cckesi);                             % ny:输出参数个数；nkesi:状态变量个数；
% [nkesi,nu] = size(Bckesi);                             % nkesi:状态变量个数；nu:控制输入参数个数；
% [nkesi,nd] = size(Bdkesi);                             % nkesi:状态变量个数；nd:干扰变量参数个数；
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