close all;
clear all;
clc;

global ToRad ToDeg;
ToDeg = 180/pi;ToRad = pi/180;

%% ��������
load('���߷�����֤���ߵ�and���ݼ�_data.mat');

%% ��ʼ��
deltat_imu = 0.004;
imu_sys_time = imu_sys_time * 1e-3;
gravity = [0 0 -1];P_acc = 0.001;
% �ο�ŷ����
roll = roll * 1e-4 * ToDeg;
Pitch = Pitch * 1e-4 * ToDeg;
yaw = yaw * 1e-4 * ToDeg;

%% ���ٶ��������޵�ͨ�˲���
Qwb = [1 0 0 0];
for t=1:length(imu_sys_time)
    % ���ݻ���Ԥ����̬
    Qwb = Qwb + 0.5 * quaternProd(Qwb, [0 gyro_x(t)*deltat_imu gyro_y(t)*deltat_imu gyro_z(t)*deltat_imu]); 
    Qwb = Qwb / norm(Qwb);
    % ���ٶ�������
    acc_b = [acc_x(t) acc_y(t) acc_z(t)]';
    Rwb = quatern2rotMat(Qwb);
    acc_w = Rwb * acc_b;
    acc_w = acc_w / norm(acc_w);
    [angle,axis] = get_included_angle_from_unit_vector(acc_w.', gravity);
    % �������Ԫ��
    angle = P_acc * angle;
    q_err = axisAngle2quatern(axis, angle);
    q_err = q_err / norm(q_err);
    % ������̬��Ԫ��
    Qwb = quaternProd(q_err, Qwb);
    Qwb = Qwb / norm(Qwb);
    
    Qwb_hat(t,:) = Qwb;
end
euler = quatern2euler(Qwb_hat) * ToDeg;

%% ��̬�ǶԱ�
figure
subplot(311)
plot(imu_sys_time, euler(:,1));
grid on;hold on;
plot(imu_sys_time, roll);
xlabel('time (s)');
ylabel('angle (deg)');
legend('roll','roll ref');
title('Roll');

subplot(312)
plot(imu_sys_time, euler(:,2));
grid on;hold on;
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
