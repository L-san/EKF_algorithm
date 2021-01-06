clear variables; clc;
Re = 6371*1000; %радиус Земли
h = 400*1000; %высота над уровнем моря
mu = 398600e+9;           %гравитационная постоянная Земли
%%
%элементы орбиты
sma = h+Re; %большая полуось орбиты
ecc = 0.0004; %эксцентриситет
inc = 51.64*pi/180; %наклонение
raan = 10.752*pi/180; %долгота восходящего узла
aop = 75.14*pi/180; %аргумент перигея
ta = 115.95*pi/180; %истинная аномалия
w0 = sqrt(mu/(sma)^3); %среднее движение
%%
global orbit_vec
    orbit_vec = [sma ecc inc raan aop ta w0]; %вектор элементов орбиты
%%
%интегратор
h = 0.001;%шаг по сетке
N = 2;
t = 0:h:N;
x = zeros(length(t));
y = zeros(7,length(t));
%%начальные условия
y0 = [1 0 0 0 0 0 0]'; x0 = 0;
y(:,1) = y0; x(1) = x0;
%%интегрируем
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