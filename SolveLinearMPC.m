function [control] = SolveLinearMPC(matrix_a,matrix_b,matrix_c,pre_matrix_c,matrix_q,matrix_r,...
    matrix_lower,matrix_upper,matrix_state_initial,pre_matrix_state_initial,reference,pre_reference,control)
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
[nkesi,nu] = size(matrix_b);                             % nkesi:number of state; nu:number of control;
ny = nkesi;
nd = size(matrix_c,1);
Ackesi = matrix_a;
Bckesi = matrix_b;
Bdkesi = eye(nd);
Cckesi = eye(nkesi);
% Nc = size(control,1)/nu;
dis = matrix_c;
Np = size(reference,1)/nkesi;
% [ny,nkesi] = size(Cckesi);                             % ny:number of output; nkesi:number of state;
% [nkesi,nu] = size(Bckesi);                             % nkesi:number of state; nu:number of control;
% [nkesi,nd] = size(Bdkesi);                             % nkesi:number of state;nd:number of disturbance;
%% Update augment reference matrix_t
matrix_t = reference;

%% Update augment control matrix_v
matrix_v = control;

matrix_a_power_cell = cell(Np,1);
for i = 1:Np
    matrix_a_power_cell{i,1} = matrix_a^i;
end

matrix_k_cell = cell(Np,Np);
for i = 1:1:Np
    for k = 1:1:Np
        if k<=i
            matrix_k_cell{i,k} = matrix_a^(i-k)*matrix_b;
        else
            matrix_k_cell{i,k} = zeros(nkesi,nu);
        end
    end
end
matrix_k = cell2mat(matrix_k_cell);

%% Initialize matrix_k, matrix_m, matrix_t and matrix_v, matrix_qq, matrix_rr, vector of matrix A power
matrix_qq_cell = cell(Np,Np);
for i = 1:Np
    for j = 1:Np
        matrix_qq_cell{i,j} = zeros(nkesi,nkesi);
    end
end
matrix_rr_cell = cell(Np,Np);
for i = 1:Np
    for j = 1:Np
        matrix_rr_cell{i,j} = zeros(nu,nu);
    end
end
matrix_ll_cell = cell(Np,1);
matrix_uu_cell = cell(Np,1);
matrix_cc_cell = cell(Np,1);
matrix_aa_cell = cell(Np,1);
% matrix_aa = [A;A^2;A^3;...,A^Np];
for i = 1:Np
    matrix_aa_cell{i,1} = matrix_a^i;
end
matrix_aa = cell2mat(matrix_aa_cell);
matrix_m = matrix_aa*matrix_state_initial;
% matrix_ll = [matrix_lower;matrix_lower;...;matrix_lower]
% matrix_uu = [matrix_upper;matrix_upper;...;matrix_upper]
% matrix_qq = diag(matrix_q,matrix_q,...,matrix_q)
% matrix_rr = diag(matrix_r,matrix_r,...,matrix_r)
% matrix_cc = [C;C+AC;C+AC+A^2*C;...;C+AC+A^2*C+...+A^(Np-1)*C];
matrix_cc_cell{1,1} = matrix_c;
for i = 1:Np
    matrix_ll_cell{i,1} = matrix_lower;
    matrix_uu_cell{i,1} = matrix_upper;
    matrix_qq_cell{i,i} = matrix_q;
    matrix_rr_cell{i,i} = matrix_r;
    if (i>1)
        matrix_cc_cell{i,1} = matrix_a*matrix_cc_cell{i-1,1}+matrix_c;
    end
end
matrix_ll = cell2mat(matrix_ll_cell);
matrix_uu = cell2mat(matrix_uu_cell);
matrix_qq = cell2mat(matrix_qq_cell);
matrix_rr = cell2mat(matrix_rr_cell);
matrix_cc = cell2mat(matrix_cc_cell);
%% Update matrix_m1,matrix_m2,convert MPC to QP
% min_x : q(x) = 0.5 * x' * matrix_m1 * x + x' * matrix_m2
matrix_m1 = matrix_k' * matrix_qq * matrix_k + matrix_rr;
matrix_m1 = (matrix_m1+matrix_m1')/2;
matrix_m2 = matrix_k' * matrix_qq * (matrix_m + matrix_cc - matrix_t);

%% Update matrix_m0,matrix_ctrl_gain,matrix_add_gain,obtain the analytical
%   control gain matrix, corresponding to the unconstrained QP problem
matrix_m0 = -matrix_m1'*matrix_k'*matrix_qq;
matrix_ctrl_gain = matrix_m0*matrix_aa;
matrix_add_gain = matrix_m0*matrix_cc;

%% Format in qp_solver
% 
%           min_x  : q(x) = 0.5 * x^T * Q * x  + x^T c
%           with respect to:  A * x = b (equality constraint)
%                             C * x <= d (inequality constraint)
%

% TODO(QiL) : change qp solver to box constraint or substitute QPOASES
% Method 1: qp
matrix_inequality_constrain_ll = eye(size(matrix_ll,1));
matrix_inequality_constrain_uu = eye(size(matrix_uu,1));
matrix_inequality_constrain = [-matrix_inequality_constrain_ll;matrix_inequality_constrain_uu];
matrix_inequality_boundary = [-matrix_ll;matrix_uu];
% caluculate QP problem
control = quadprog(matrix_m1,matrix_m2,matrix_inequality_constrain,matrix_inequality_boundary);
end
