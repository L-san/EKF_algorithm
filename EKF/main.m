clear variables; clc;
load('data.mat');
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
       satellite.Sm = 0.1*0.2; %drag area
       satellite.ra = [0; 0.001; 0]; % center of mass shift
       satellite.I = diag([0.01, 0.01, 0.02]);
delta_T = 0.01;          %time parameter
simT = 5;
simConstT = simT/delta_T;
Rc = 0.01;
Qc = 1;
[x_hat, time] = calcualtePosition(h,x,Rc,Qc,y);
subplot(1,3,1); plot(time,y(5,:),time, x_hat(5,:)); xlabel("Time,s"); ylabel("w_x, rad/s"); grid;
subplot(1,3,2); plot(time,y(6,:),time, x_hat(6,:)); xlabel("Time,s"); ylabel("w_y, rad/s"); grid;
subplot(1,3,3); plot(time,y(7,:),time, x_hat(7,:)); xlabel("Time,s"); ylabel("w_z, rad/s"); grid;