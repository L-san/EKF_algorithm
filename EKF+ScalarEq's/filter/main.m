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
Q = 0.1;

q = y(1:4,:);
w = y(5:7,:);
r = y(8:10,:);
V = y(11:13,:);
dt = h;
sensor_data = [B_b; wb];
x_hat = zeros(7,length(x));
h = zeros(7,length(x));
P = eye(7);
x_hat(:,1) = [q(:,1); w(:,1)];
for i = 1:length(x)
    
    qw = [q(:,i); w(:,i)];
    %
    Ain2orb = in2orb(r(:,i),V(:,i));
    B_I = magneticField(r(:,i));
    Borb = Ain2orb*B_I;
    %
    JF = stateTransitionMatrix(qw,dt);
    JH = observationMatrix(qw, Borb);
    f = getAttitudeVector([qw; r(:,i); V(:,i)], dt);
    h = getHiddenState([qw; r(:,i); V(:,i)], Borb);
    
    z = [0;sensor_data(:,i)];
    
    [x_hat(:,i), P] = kalmanReal(R, Q, z, P, JF, JH, h, f);
end
figure;
subplot(1,3,1); plot(x, wb(1,:), x, x_hat(5,:)); xlabel("Время, с"); ylabel("w_x,  рад/c"); grid; legend("raw","filtered");
subplot(1,3,2); plot(x, wb(2,:), x, x_hat(5,:)); xlabel("Время, с"); ylabel("w_y,  рад/c"); grid; legend("raw","filtered");
subplot(1,3,3); plot(x, wb(3,:), x, x_hat(7,:)); xlabel("Время, с"); ylabel("w_z,  рад/c"); grid; legend("raw","filtered");

figure;
subplot(2,2,1); plot(x,y(1,:), x, x_hat(1,:)); grid; legend("raw","filtered");
subplot(2,2,2); plot(x,y(2,:), x, x_hat(2,:)); grid; legend("raw","filtered");
subplot(2,2,3); plot(x,y(3,:), x, x_hat(3,:)); grid; legend("raw","filtered");
subplot(2,2,4); plot(x,y(4,:), x, x_hat(4,:)); grid; legend("raw","filtered");


                     
                     
 