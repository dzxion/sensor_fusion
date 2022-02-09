close all;
clear all;
clc;

global ToRad ToDeg;
ToDeg = 180/pi;ToRad = pi/180;
GravityAcc = 980.665;

%% 导入数据
load('航线飞行验证纯惯导and数据集_data.mat');

%% 初始化
deltat_imu = 0.004;deltat_gps = 0.2;

gravity = [0 0 -1];
P_acc = 0.005;P_yaw = 0.1;
pitch_acc_filter = zeros(length(imu_sys_time), 1);roll_acc_filter = zeros(length(imu_sys_time), 1);
pitch_acc_compensate = zeros(length(imu_sys_time), 1);roll_acc_compensate = zeros(length(imu_sys_time), 1);
acc_motion_gps = zeros(length(imu_sys_time), 3);
history = 72;delay_acc = 72;delay_yaw_gps = 21;

gps_location_lat = gps_location_lat * 1e-7;
gps_location_lon = gps_location_lon * 1e-7;
gps_location_alt = gps_location_alt * 1e-1;
gps_head = gps_head * 1e-2;
[gps_x,gps_y,gps_z] = map_projection_project( gps_location_lat, gps_location_lon, gps_location_alt);

Qwb_his = zeros(history, 4);
Qwb_hat = zeros(length(imu_sys_time), 4);

% 参考欧拉角
roll = roll * 1e-4 * ToDeg;
Pitch = Pitch * 1e-4 * ToDeg;
yaw = yaw * 1e-4 * ToDeg;

%% 加速度滤波
fs=250;fc=5;
[b,a] = butter(2,fc/(fs/2));
acc_filter_x = filter(b,a,acc_x);
acc_filter_y = filter(b,a,acc_y);
acc_filter_z = filter(b,a,acc_z);

%% 加速度算roll pitch
for t=1:length(imu_sys_time)
    acc_filter_b = [acc_filter_x(t) acc_filter_y(t) acc_filter_z(t)];
    acc_filter_b = acc_filter_b / norm(acc_filter_b);
    
    pitch_acc_filter(t) = asin(acc_filter_b(1)) * ToDeg;
    roll_acc_filter(t) = atan2(-acc_filter_b(2), -acc_filter_b(3)) * ToDeg;
end

%% gps速度首次更新时间计算
gps_count = 1;
while 1
    gps_count = gps_count + 1;
    if gps_vitow(gps_count) ~= gps_vitow(gps_count-1) 
         break;
    end
end
gps_vel_first_update_time = gps_sys_time(gps_count);

%% gps运动加速度延时补偿
Qwb = [1 0 0 0];
gps_update = false;imu_count = 0;
vel_gps_last = [gps_vel_n(1),gps_vel_e(1),gps_vel_d(1)];vel_gps = [gps_vel_n(1),gps_vel_e(1),gps_vel_d(1)];
acc_motion_w = [0,0,0];acc_motion_w_last = [0,0,0];delta_acc_w = [0,0,0];
acceleration_available = false;gps_yaw_available = false;
error_rp = [];error_yaw = [];
% 融合算法
for k=1:length(imu_sys_time)
    % 信息处理
    imu_count = imu_count + 1;
    % gps速度首次更新
    if imu_sys_time(k) == gps_vel_first_update_time  
        imu_count = 0;
        vel_gps = [gps_vel_n(gps_count),gps_vel_e(gps_count),gps_vel_d(gps_count)];
        pos_gps = [gps_x(gps_count),gps_y(gps_count),gps_z(gps_count)];
        head_gps = gps_head(gps_count)*ToRad;
        acc_motion_w = (vel_gps - vel_gps_last) / deltat_gps;
        delta_acc_w = acc_motion_w - acc_motion_w_last;
        gps_yaw_available = true;
%         gps_update = true;
    end
    if imu_count == 50
        while 1
            gps_count = gps_count + 1;
            if gps_vitow(gps_count) ~= gps_vitow(gps_count-1) 
                break;
            end
        end
        vel_gps_last = vel_gps;
        vel_gps = [gps_vel_n(gps_count),gps_vel_e(gps_count),gps_vel_d(gps_count)];
        pos_gps = [gps_x(gps_count),gps_y(gps_count),gps_z(gps_count)];
        head_gps = gps_head(gps_count)*ToRad;
        acc_motion_w_last = acc_motion_w;
        acc_motion_w = (vel_gps - vel_gps_last) / deltat_gps;
        delta_acc_w = acc_motion_w - acc_motion_w_last;
        imu_count = 0;
        gps_yaw_available = true;
%         acceleration_available = true;
%         gps_update = true;
    end
    acc_motion_w_interp = acc_motion_w_last + (imu_count/50) * delta_acc_w;
    if k == delay_acc+1
        acceleration_available = true;
    end
    
    % 信息融合
    % 陀螺预测
    Qwb = Qwb + 0.5 * quaternProd(Qwb, [0 gyro_x(k)*deltat_imu gyro_y(k)*deltat_imu gyro_z(k)*deltat_imu]); 
    Qwb = Qwb / norm(Qwb);
    
    % 加速度修正
    if acceleration_available
        acc_b_filter = [acc_filter_x(k-delay_acc) acc_filter_y(k-delay_acc) acc_filter_z(k-delay_acc)]';
        Rwb = quatern2rotMat(Qwb_his(1,:));
        acc_w_filter = Rwb * acc_b_filter;
        acc_ref = acc_motion_w_interp + gravity * GravityAcc;
        acc_ref = acc_ref / norm(acc_ref);
        [angle_error,axis] = get_included_angle_from_unit_vector(acc_w_filter.', acc_ref);
        error_rp = [error_rp,angle_error*ToDeg];
        angle_error = P_acc * angle_error;
        q_correction = axisAngle2quatern(axis, angle_error);
        q_correction = q_correction / norm(q_correction);
        for i=1:history
            Qwb_his(i,:) = quaternProd(q_correction, Qwb_his(i,:));
            Qwb_his(i,:) = Qwb_his(i,:) / norm(Qwb_his(i,:));
        end
        Qwb = quaternProd(q_correction, Qwb);
        Qwb = Qwb / norm(Qwb);
%         Qwb_his(history) = Qwb;
    end
    
    % gps偏航修正
    if gps_yaw_available
        index = history-delay_yaw_gps+1;
        euler_pred = quatern2euler(Qwb_his(index,:));
        yaw_pred = euler_pred(3);
        % 角度误差化为-180~180
        angle_error = head_gps - yaw_pred;
        while angle_error < -pi
            angle_error = angle_error + 2*pi;
        end
        while angle_error > pi
            angle_error = angle_error - 2*pi;
        end
        error_yaw = [error_yaw,angle_error*ToDeg];
        angle_error = P_yaw * angle_error;
        axis = [0,0,1];
        q_correction = axisAngle2quatern(axis, angle_error);
        q_correction = q_correction / norm(q_correction);
        for i=index:history
            Qwb_his(i,:) = quaternProd(q_correction, Qwb_his(i,:));
            Qwb_his(i,:) = Qwb_his(i,:) / norm(Qwb_his(i,:));
        end
        Qwb = quaternProd(q_correction, Qwb);
        Qwb = Qwb / norm(Qwb);
        gps_yaw_available = false;
        %          gps_update = false;
    end
    
    % 历史状态滚动
    for i=1:history-1
        Qwb_his(i,:) = Qwb_his(i+1,:);
    end
    Qwb_his(history,:) = Qwb;
    Qwb_hat(k,:) = Qwb;
end
euler = quatern2euler(Qwb_hat) * ToDeg; % 四元数转欧拉角

%% 加速度算roll pitch（加入运动加速度补偿）
for t=1:length(imu_sys_time)
    acc_b = [acc_x(t) acc_y(t) acc_z(t)];
    acc_b = acc_b / norm(acc_b);
    acc_filter_b = [acc_filter_x(t) acc_filter_y(t) acc_filter_z(t)];
    
    Qwb = Qwb_hat(t,:);
    Qbw = quaternConj(Qwb);
    Rbw = quatern2rotMat(Qbw);
    acc_motion = acc_motion_gps(t,:);
    acc_motion = acc_motion.';
    acc_motion_b = Rbw * acc_motion;
    acc_motion_b = acc_motion_b.';
    acc_compensate = acc_filter_b * GravityAcc - acc_motion_b;
    acc_compensate = acc_compensate/norm(acc_compensate);
    
    pitch_acc_compensate(t) = asin(acc_compensate(1)) * ToDeg;
    roll_acc_compensate(t) = atan2(-acc_compensate(2), -acc_compensate(3)) * ToDeg;
end

imu_sys_time = imu_sys_time * 1e-3;
gps_sys_time = gps_sys_time * 1e-3;

%% 姿态角对比
figure
subplot(311)
plot(imu_sys_time, euler(:,1));
grid on;hold on;
% plot(imu_sys_time, roll_acc_filter);
% grid on;hold on;
% plot(imu_sys_time, roll_acc_compensate);
% grid on;hold on;
plot(imu_sys_time, roll);
xlabel('time (s)');
ylabel('angle (deg)');
legend('roll','roll ref');
title('Roll');

subplot(312)
plot(imu_sys_time, euler(:,2));
grid on;hold on;
% plot(imu_sys_time, pitch_acc_filter);
% grid on;hold on;
% plot(imu_sys_time, pitch_acc_compensate);
% grid on;hold on;
plot(imu_sys_time, Pitch);
xlabel('time (s)');
ylabel('angle (deg)');
legend('pitch','pitch ref');
title('Pitch');

subplot(313)
plot(imu_sys_time, euler(:,3));
grid on;hold on;
plot(imu_sys_time, yaw);
xlabel('time (s)');
ylabel('angle (deg)');
legend('yaw','yaw ref');
title('Yaw');
