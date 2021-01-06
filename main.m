clear variables; clc;
Re = 6371*1000; %������ �����
h = 400*1000; %������ ��� ������� ����
mu = 398600e+9;           %�������������� ���������� �����
%%
%�������� ������
sma = h+Re; %������� ������� ������
ecc = 0.0004; %��������������
inc = 51.64*pi/180; %����������
raan = 10.752*pi/180; %������� ����������� ����
aop = 75.14*pi/180; %�������� �������
ta = 115.95*pi/180; %�������� ��������
w0 = sqrt(mu/(sma)^3); %������� ��������
%%
global orbit_vec
    orbit_vec = [sma ecc inc raan aop ta w0]; %������ ��������� ������
%%
%����������
h = 0.001;%��� �� �����
N = 2;
t = 0:h:N;
x = zeros(length(t));
y = zeros(7,length(t));
%%��������� �������
y0 = [1 0 0 0 0 0 0]'; x0 = 0;
y(:,1) = y0; x(1) = x0;
%%�����������
[x,y] = ode4(@motionEquations,h,[0,N],[x0;y0]);
for i = 1:length(t)
    [a(i) b(i) c(i)] = quat2angle(y(1:4,i)','XYX');
    q_norm(i) = sqrt(y(1,i)^2+y(2,i)^2+y(3,i)^2+y(4,i)^2);
end
figure;
subplot(2,3,1); plot(x,a);xlabel("Time,s"); ylabel("roll, rad");grid;
subplot(2,3,2); plot(x,b);xlabel("Time,s"); ylabel("pitch, rad");grid;
subplot(2,3,3); plot(x,c);xlabel("Time,s"); ylabel("yaw, rad");grid;

subplot(2,3,4); plot(x, y(5,:)); xlabel("Time,s"); ylabel("w_x, rad/s"); grid;
subplot(2,3,5); plot(x, y(6,:)); xlabel("Time,s"); ylabel("w_y, rad/s"); grid;
subplot(2,3,6); plot(x, y(7,:)); xlabel("Time,s"); ylabel("w_z, rad/s"); grid;

figure; plot(x,q_norm); xlabel("Time,s"); ylabel("quaternion norm"); grid;