clear variables; clc;
load("data.mat");
global mu
        mu = 398600e+9; %gravitational constant
%%
%orbital elemets
%%
global orbit_vec
       orbit_vec.Re = 6371*1000; %earth radius
       orbit_vec.h = 400*1000; %height above sea level
       orbit_vec.sma =  orbit_vec.h+orbit_vec.Re; %semimajor axis
       orbit_vec.ecc = 0.0004; %eccentricity
       orbit_vec.inc = 51.64*pi/180; %inclination
       orbit_vec.raan = 10.752*pi/180; %right ascension of the ascending node
       orbit_vec.aop = 75.14*pi/180; %argument of perigee
       orbit_vec.ta = 115.95*pi/180; %truth anomaly
       orbit_vec.w0 = sqrt(mu/(orbit_vec.sma)^3); %mean motion

global satellite
       satellite.c0 = 2.2; %drag coefficient
       satellite.Sm = 0.1*0.3; %drag area
       satellite.ra = [0; 0.001; 0]; % center of mass shift
       satellite.I = diag([0.01, 0.01, 0.02]);
%%
R = 100;
Q = 0;

q = y(1:4,:);
w = y(5:7,:);
r = y(8:10,:);
V = y(11:13,:);
dt = h;

sensor_data = [B_b; wb];
x_hat = zeros(7,length(x));

Jfx = zeros(7,length(x));
Jhx = zeros(7,length(x));
fx = zeros(7,length(x));
hx = zeros(7,length(x));
P = eye(7);
x_hat(:,1) = [q(:,1); w(:,1)];

for i = 2:length(x)
    %x_hat(:,i-1) = y(1:7,i-1);
    Ain2orb = in2orb(r(:,i-1),V(:,i-1));
    B_I = magneticField(r(:,i-1));
    Borb = Ain2orb*B_I;
    %
    JF = stateTransitionMatrix([x_hat(:,i-1); r(:,i-1); V(:,i-1)],dt);
    JH = observationMatrix(x_hat(:,i-1), Borb);
    f = getAttitudeVector([x_hat(:,i-1); r(:,i-1); V(:,i-1)], dt);
    h = getHiddenState([x_hat(:,i-1); r(:,i-1); V(:,i-1)], Borb);
    Jfx(:,i)= JF*x_hat(:,i-1);
    Jhx(:,i)= JH*x_hat(:,i-1);
    fx(:,i) = f;
    hx(:,i) = h;
    z = [0; sensor_data(:,i)];
    [x_hat(:,i), P] = kalmanReal(R, Q, z, P, JF, JH, h, f);
end

figure;
title("fx, Jfx");
subplot(2,4,1); plot(x,fx(1,:), x, Jfx(1,:), x, y(1,:)); grid; legend("fx","Jfx","y");
subplot(2,4,2); plot(x,fx(2,:), x, Jfx(2,:), x, y(2,:)); grid; legend("fx","Jfx","y");
subplot(2,4,3); plot(x,fx(3,:), x, Jfx(3,:), x, y(3,:)); grid; legend("fx","Jfx","y");
subplot(2,4,4); plot(x,fx(4,:), x, Jfx(4,:), x, y(4,:)); grid; legend("fx","Jfx","y");
subplot(2,4,5); plot(x,fx(5,:), x, Jfx(5,:), x, y(5,:)); grid; legend("fx","Jfx","y");
subplot(2,4,6); plot(x,fx(6,:), x, Jfx(6,:), x, y(6,:)); grid; legend("fx","Jfx","y");
subplot(2,4,7); plot(x,fx(7,:), x, Jfx(7,:), x, y(7,:)); grid; legend("fx","Jfx","y");

figure;
title("hx, Jhx, z");
subplot(2,4,1); plot(x,hx(1,:), x, Jhx(1,:), x, zeros(1,length(sensor_data(1,:)))); legend("hx","Jhx","z"); grid; 
subplot(2,4,2); plot(x,hx(2,:), x, Jhx(2,:), x, sensor_data(1,:)); legend("hx","Jhx","z"); grid;
subplot(2,4,3); plot(x,hx(3,:), x, Jhx(3,:), x, sensor_data(2,:)); legend("hx","Jhx","z"); grid; 
subplot(2,4,4); plot(x,hx(4,:), x, Jhx(4,:), x, sensor_data(3,:)); legend("hx","Jhx","z"); grid; 
subplot(2,4,5); plot(x,hx(5,:), x, Jhx(5,:), x, sensor_data(4,:)); grid; legend("hx","Jhx","z");
subplot(2,4,6); plot(x,hx(6,:), x, Jhx(6,:), x, sensor_data(5,:)); grid; legend("hx","Jhx","z");
subplot(2,4,7); plot(x,hx(7,:), x, Jhx(7,:), x, sensor_data(6,:)); grid; legend("hx","Jhx","z");



figure;
subplot(2,4,1); plot(x,y(1,:), x, x_hat(1,:)); xlabel("Время, с"); ylabel("q0"); legend("данные","фильтр"); grid;
subplot(2,4,2); plot(x,y(2,:), x, x_hat(2,:)); xlabel("Время, с"); ylabel("q1"); legend("данные","фильтр"); grid;
subplot(2,4,3); plot(x,y(3,:), x, x_hat(3,:)); xlabel("Время, с"); ylabel("q2"); legend("данные","фильтр"); grid;
subplot(2,4,4); plot(x,y(4,:), x, x_hat(4,:)); xlabel("Время, с"); ylabel("q3"); legend("данные","фильтр"); grid; 
subplot(2,4,5); plot(x, y(5,:)*180/pi, x, x_hat(5,:)*180/pi); xlabel("Время, с"); ylabel("w_x,  градусы/c"); grid; legend("данные","фильтр");
subplot(2,4,6); plot(x, y(6,:)*180/pi, x, x_hat(5,:)*180/pi); xlabel("Время, с"); ylabel("w_y,  градусы/c"); grid; legend("данные","фильтр");
subplot(2,4,7); plot(x, y(7,:)*180/pi, x, x_hat(7,:)*180/pi); xlabel("Время, с"); ylabel("w_z,  градусы/c"); grid; legend("данные","фильтр");
subplot(2,4,8); plot(x,vecnorm(x_hat(1:4,:)));legend("норма кватерниона");grid; 

                     
                     
 