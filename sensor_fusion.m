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
% �ο�ŷ����
roll = roll * 1e-4 * ToDeg;
Pitch = Pitch * 1e-4 * ToDeg;
yaw = yaw * 1e-4 * ToDeg;

%% ��������̬
Qwb = [1 0 0 0];
for t=1:length(imu_sys_time)
     Qwb = Qwb + 0.5 * quaternProd(Qwb, [0 gyro_x(t)*deltat_imu gyro_y(t)*deltat_imu gyro_z(t)*deltat_imu]); 
     Qwb = Qwb / norm(Qwb);
     Qwb_gyro(t,:) = Qwb;
end
euler_gyro = quatern2euler(Qwb_gyro) * ToDeg;

%% ��̬�ǶԱ�
figure
subplot(311)
plot(imu_sys_time, euler_gyro(:,1));
grid on;hold on;
plot(imu_sys_time, roll);
xlabel('time (s)');
ylabel('angle (deg)');
legend('roll gyro','roll ref');
title('Roll');

subplot(312)
plot(imu_sys_time, euler_gyro(:,2));
grid on;hold on;
plot(imu_sys_time, Pitch);
xlabel('time (s)');
ylabel('angle (deg)');
legend('pitch gyro','pitch ref');
title('Pitch');

subplot(313)
plot(imu_sys_time, euler_gyro(:,3));
grid on;hold on;
plot(imu_sys_time, yaw);
xlabel('time (s)');
ylabel('angle (deg)');
legend('yaw gyro','yaw ref');
title('Yaw');
