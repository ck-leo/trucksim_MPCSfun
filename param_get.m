clc;
%% configure vehicle parameter
basic_state_size_ = 4;% number of state:lateral error,lateral error rate,heading error, heading error rate
controls_ = 1;% number of controls:delta_f
horizon_ = 20;  % Np
vertical_ = horizon_;% Nc
M_SU = 4455;M_US1 = 570;M_US2 = 735;
mass_ = M_SU+M_US1+M_US2;
lf_ = 1110/1000;
lr_ = 2790/1000;
cf_ = 2 * (20164.4-15677.2)/(2*pi/180);
cr_ = cf_;
iz_ = 34802.6;
%% continuous model:dx/dt = A * x + B * u + C
matrix_a_ = zeros(basic_state_size_,basic_state_size_);
matrix_a_coeff_ = zeros(basic_state_size_,basic_state_size_);
matrix_b_ = zeros(basic_state_size_,controls_);
matrix_d_ = zeros(basic_state_size_,1);
pre_matrix_d_ = zeros(basic_state_size_,1);
matrix_state_ = zeros(basic_state_size_,1);
pre_matrix_state_ = zeros(basic_state_size_,1);
matrix_a_(1, 2) = 1.0;
matrix_a_(2, 3) = (cf_ + cr_) / mass_;
matrix_a_(3, 4) = 1.0;
matrix_a_(4, 3) = (lf_ * cf_ - lr_ * cr_) / iz_;
matrix_a_coeff_(2, 2) = -(cf_ + cr_) / mass_;
matrix_a_coeff_(2, 4) = (lr_ * cr_ - lf_ * cf_) / mass_;
matrix_a_coeff_(3, 4) = 1.0;
matrix_a_coeff_(4, 2) = (lr_ * cr_ - lf_ * cf_) / iz_;
matrix_a_coeff_(4, 4) = -1.0 * (lf_ * lf_ * cf_ + lr_ * lr_ * cr_) / iz_;
matrix_b_(2, 1) = cf_ / mass_;
matrix_b_(4, 1) = lf_ * cf_ / iz_;
% weight param
% test
% weight_lateral_error = 1;
% weight_lateral_error_rate = 0;
% weight_heading_error = 1;
% weight_heading_error_rate = 0;
% weight_steer = 2;
% Np = 20;
% dt = 0.01;
% linear_v = 70/3.6;
% A_open = getAopen(matrix_a_,matrix_a_coeff_,matrix_b_,matrix_q_,matrix_r_,Np,dt,linear_v);
% eig_real = real(eig(A_open));

weight_lateral_error = [1:10];
weight_lateral_error_rate = [1:10];
weight_heading_error = [1:10];
weight_heading_error_rate = [1:10];
weight_steer = [1:10];
Np = [1:2:20];
dt = [0.01:0.01:0.1];
linear_v = 70/3.6;
m = 1;
time_toc = zeros(1,3);
for i1 = 1:length(weight_lateral_error)
    for i2 = 1:length(weight_lateral_error_rate)
        for i3 = 1:length(weight_heading_error)
            for i4 = 1:length(weight_heading_error_rate)
                for iu = 1:length(weight_steer)
                    for iNp = 1:length(Np)
                        for idt = 1:length(dt)
                            tic;
                            matrix_q_ = diag([weight_lateral_error(i1),weight_lateral_error_rate(i2),weight_heading_error(i3),...
                                weight_heading_error_rate(i4)]);
                            matrix_r_ = weight_steer(iu);
                            A_open = getAopen(matrix_a_,matrix_a_coeff_,matrix_b_,matrix_q_,matrix_r_,Np(iNp),dt(idt),linear_v);
                            eig_real = real(eig(A_open));
                            if(sum(eig_real<0)==4 || norm(eig_real) <= 1e-3)
                                 weight_all(m,:) = [weight_lateral_error(i1),weight_lateral_error_rate(i2),weight_heading_error(i3),...
                                     weight_heading_error_rate(i4),weight_steer(iu),Np(iNp),dt(idt)];
                                 eig_real_all(m,:) = eig(A_open);
                                 m = m+1;
                            end
                            time_toc(1) = time_toc(1)+1;
                            time_toc(2) = time_toc(2) + toc;
                            time_toc(3) = time_toc(2)/time_toc(1);
                        end
                    end
                end
            end
        end
    end
end
save('weight_and_eig.mat','weight_all','eig_real_all');
time = time_toc(2)

