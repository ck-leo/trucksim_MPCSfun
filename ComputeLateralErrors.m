function debug = ComputeLateralErrors(u,trajectory_analyzer_)
x0 = u(1);
y0 = u(2);
theta0 = u(3);
linear_v = u(4);
anglar_v = u(5);
X_ref = trajectory_analyzer_(:,1);
Y_ref = trajectory_analyzer_(:,2);
theta = trajectory_analyzer_(:,3);
kappa = trajectory_analyzer_(:,4);
N_ref = length(X_ref);
index_min = 1;
dis_min = sqrt((x0-X_ref(1))^2+(y0-Y_ref(1))^2);
% find the nearest point
for i = 1:N_ref
    dis_temp = sqrt((x0-X_ref(i))^2+(y0-Y_ref(i))^2);
    if(dis_temp<dis_min)
        index_min = i;
        dis_min = dis_temp;
    end
end
%% 路径点较密集
debug = zeros(4,1);
dx = x0 - X_ref(index_min);
dy = y0 - Y_ref(index_min);
delta_theta = NormalizeAngle(theta0-theta(index_min));
debug(1) = dy*cos(theta(index_min)) - dx*sin(theta(index_min)); % lateral_error
debug(2) = linear_v*sin(delta_theta);             % lateral_error_rate
debug(3) = delta_theta;                                        % heading_error
debug(4) = anglar_v - linear_v*kappa(index_min);               % heading_error_rate
%% 如果路径点较稀疏
% if index_min==N_ref
%     x1 = X_ref(index_min-1);y1 = Y_ref(index_min-1);
%     x2 = X_ref(index_min);y2 = Y_ref(index_min);
% else
%     x1 = X_ref(index_min);y1 = Y_ref(index_min);
%     x2 = X_ref(index_min+1);y2 = Y_ref(index_min+1);
% end
% % 求u到x1x2直线的垂足
% A = y2-y1;B = x1-x2;C = x2*y1-x1*y2;
% debug(1) = (B*B*x0 - A*B*y0 - A*C) / (A*A + B*B);
% debug(2) = (-A*B*x0 + A*A*y0 - B*C) / (A*A + B*B);
% debug(3) = atan2(A,-B);