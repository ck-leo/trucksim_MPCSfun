function [A_open] = getAopen(matrix_a_,matrix_a_coeff_,matrix_b_,matrix_q_,matrix_r_,Np,ts_,linear_v)
% get the matrix A-B*Kmpc
% 
[nkesi,nu] = size(matrix_b_);                             % nkesi:number of state; nu:number of control;
matrix_a_(2, 2) = matrix_a_coeff_(2, 2) / linear_v;
matrix_a_(2, 4) = matrix_a_coeff_(2, 4) / linear_v;
matrix_a_(4, 2) = matrix_a_coeff_(4, 2) / linear_v;
matrix_a_(4, 4) = matrix_a_coeff_(4, 4) / linear_v;
matrix_ad_ = inv(eye(nkesi) - ts_ * 0.5 * matrix_a_) * (eye(nkesi) + ts_ * 0.5 * matrix_a_);
matrix_bd_ = matrix_b_ * ts_;

% Phi
matrix_k_cell = cell(Np,Np);
for i = 1:1:Np
    for k = 1:1:Np
        if k<=i
            matrix_k_cell{i,k} = matrix_ad_^(i-k)*matrix_bd_;
        else
            matrix_k_cell{i,k} = zeros(nkesi,nu);
        end
    end
end
matrix_k = cell2mat(matrix_k_cell);
%% Initialize matrix_k, matrix_m, matrix_t and matrix_v, matrix_qq, matrix_rr, vector of matrix A power
matrix_qq_cell = cell(Np,Np);
matrix_rr_cell = cell(Np,Np);
matrix_aa_cell = cell(Np,1);
% matrix_aa = [A;A^2;A^3;...,A^Np];
% matrix_qq = diag(matrix_q,matrix_q,...,matrix_q)
% matrix_rr = diag(matrix_r,matrix_r,...,matrix_r)
for i = 1:Np
    for j = 1:Np
        matrix_qq_cell{i,j} = zeros(nkesi,nkesi);
        matrix_rr_cell{i,j} = zeros(nu,nu);
    end
end
for i = 1:Np
    matrix_aa_cell{i,1} = matrix_ad_^i;
    matrix_qq_cell{i,i} = matrix_q_;
    matrix_rr_cell{i,i} = matrix_r_;
end
matrix_aa = cell2mat(matrix_aa_cell);
matrix_qq = cell2mat(matrix_qq_cell);
matrix_rr = cell2mat(matrix_rr_cell);
%% Update matrix_m1,matrix_m2,convert MPC to QP
% min_x : q(x) = 0.5 * x' * matrix_m1 * x + x' * matrix_m2
matrix_m1 = matrix_k' * matrix_qq * matrix_k + matrix_rr;
KMPC = inv(matrix_m1)*(matrix_k'*matrix_qq*matrix_aa);
Kmpc = KMPC(1,:);
A_open = matrix_ad_-matrix_bd_*Kmpc;
end

