close all;
clear all;
clc;

global ToRad ToDeg;
ToDeg = 180/pi;ToRad = pi/180;

%% 导入数据
load('航线飞行验证纯惯导and数据集_data.mat');

%% 初始化
deltat_imu = 0.004;deltat_gps = 0.2;
imu_sys_time = imu_sys_time * 1e-3;
gps_sys_time = gps_sys_time * 1e-3;

gravity = [0 0 -1];
P_acc = 0.005;GravityAcc = 980.665;

gps_location_lat = gps_location_lat * 1e-7;
gps_location_lon = gps_location_lon * 1e-7;
gps_location_alt = gps_location_alt * 1e-1;
gps_head = gps_head * 1e-2;
[gps_x,gps_y,gps_z] = map_projection_project( gps_location_lat, gps_location_lon, gps_location_alt);

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
    acc_b = [acc_x(t) acc_y(t) acc_z(t)];
    acc_b = acc_b / norm(acc_b);
    acc_filter_b = [acc_filter_x(t) acc_filter_y(t) acc_filter_z(t)];
    acc_filter_b = acc_filter_b / norm(acc_filter_b);
    
    pitch_acc(t) = asin(acc_b(1)) * ToDeg;
    roll_acc(t) = atan2(-acc_b(2), -acc_b(3)) * ToDeg;
    pitch_acc_filter(t) = asin(acc_filter_b(1)) * ToDeg;
    roll_acc_filter(t) = atan2(-acc_filter_b(2), -acc_filter_b(3)) * ToDeg;
end

%% gps运动加速度补偿
Qwb = [1 0 0 0];
gps_count = 1;gps_update = false;imu_count = 0;
acc_motion_w = [0,0,0];vel_gps_prev = [0,0,0];vel_gps = [0,0,0];
for t=1:length(imu_sys_time)
    imu_count = imu_count + 1;
    % gps首次更新
    if imu_sys_time(t) == gps_sys_time(1) 
        imu_count = 0;
        gps_count = 1;
        vel_gps = [gps_vel_n(1),gps_vel_e(1),gps_vel_d(1)];
        pos_gps = [gps_x(1),gps_y(1),gps_z(1)];
        head_gps = gps_head(1)*ToRad;
        gps_update = true;
    end
    % gps更新
    if imu_count == 50
        while 1
            gps_count = gps_count + 1;
            if gps_vitow(gps_count) ~= gps_vitow(gps_count-1) 
                break;
            end
        end
        vel_gps_prev = vel_gps;
        vel_gps = [gps_vel_n(gps_count),gps_vel_e(gps_count),gps_vel_d(gps_count)];
        pos_gps = [gps_x(gps_count),gps_y(gps_count),gps_z(gps_count)];
        head_gps = gps_head(gps_count)*ToRad;
        imu_count = 0;
        gps_update = true;
    end
    
    % 陀螺积分预测姿态
    Qwb = Qwb + 0.5 * quaternProd(Qwb, [0 gyro_x(t)*deltat_imu gyro_y(t)*deltat_imu gyro_z(t)*deltat_imu]); 
    Qwb = Qwb / norm(Qwb);
    % 加速度求误差角
    acc_b_filter = [acc_filter_x(t) acc_filter_y(t) acc_filter_z(t)]';
    Rwb = quatern2rotMat(Qwb);
    acc_w_filter = Rwb * acc_b_filter;
    acc_w_filter = acc_w_filter / norm(acc_w_filter);
    % gps速度差分运动加速度
    if gps_update
        acc_motion_w = (vel_gps - vel_gps_prev) / deltat_gps;
    end
    acc_motion_gps(t,:) = acc_motion_w;
    acc_ref = acc_motion_w + gravity * GravityAcc;
    acc_ref = acc_ref / norm(acc_ref);
    [angle,axis] = get_included_angle_from_unit_vector(acc_w_filter.', acc_ref);
    % 求误差四元数
    angle = P_acc * angle;
    q_err = axisAngle2quatern(axis, angle);
    q_err = q_err / norm(q_err);
    % 修正姿态四元数
    Qwb = quaternProd(q_err, Qwb);
    Qwb = Qwb / norm(Qwb);
    % gps修yaw
    if gps_update
         % 四元数分解
         euler_pred = quatern2euler(Qwb);
         yaw_pred = euler_pred(3);
         q_yaw = axisAngle2quatern([0,0,1], yaw_pred);
         q_yaw = quaternConj(q_yaw);
         q_pr = quaternProd(q_yaw, Qwb);
         % 使用gps航向
         q_yaw_hat = axisAngle2quatern([0,0,1], head_gps);
         Qwb = quaternProd(q_yaw_hat, q_pr);
         gps_update = false;
    end
    
    Qwb_hat(t,:) = Qwb;
end
euler = quatern2euler(Qwb_hat) * ToDeg;

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

%% 姿态角对比
figure
subplot(311)
plot(imu_sys_time, euler(:,1));
grid on;hold on;
plot(imu_sys_time, roll_acc_filter);
grid on;hold on;
plot(imu_sys_time, roll_acc_compensate);
grid on;hold on;
plot(imu_sys_time, roll);
xlabel('time (s)');
ylabel('angle (deg)');
legend('roll','roll acc','roll acc compensate','roll ref');
title('Roll');

subplot(312)
plot(imu_sys_time, euler(:,2));
grid on;hold on;
plot(imu_sys_time, pitch_acc_filter);
grid on;hold on;
plot(imu_sys_time, pitch_acc_compensate);
grid on;hold on;
plot(imu_sys_time, Pitch);
xlabel('time (s)');
ylabel('angle (deg)');
legend('pitch','pitch acc','pitch acc compensate','pitch ref');
title('Pitch');

subplot(313)
plot(imu_sys_time, euler(:,3));
grid on;hold on;
plot(imu_sys_time, yaw);
xlabel('time (s)');
ylabel('angle (deg)');
legend('yaw','yaw ref');
title('Yaw');
