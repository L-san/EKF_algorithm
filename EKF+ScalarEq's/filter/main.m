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
q = y(1:4,:);
w = y(5:7,:);
r = y(8:10,:);
V = y(11:13,:);
dt = 0.01;

sensor_data = [B_b; wb];
x_hat = zeros(6,length(x));
fx = zeros(6,length(x));
hx = zeros(6,length(x));
Jfx = zeros(6,length(x));
Jhx = zeros(6,length(x));

normP = zeros(1,length(x));
P(:,:,1) = eye(6);

%x_hat(:,1) = [wb(:,1)/norm(wb(:,1))*sin(norm(wb(:,1))/2*dt); wb(:,1)];
%x_hat(:,1) = [[0 0 0]'; wb(:,1)];
x_hat(:,1) = y(2:7,1);
sigma_b = var(B_b');
sigma_wb = var(wb');

%R = diag([sigma_b, sigma_wb]);
R = eye(6,6)*0;
Q = 0;
%t = floor(length(x));
t = 100;
for i = 2:t
    x_hat(:,i-1) = y(2:7,i-1);
    Ain2orb = in2orb(r(:,i-1),V(:,i-1));
    B_I = magneticField(r(:,i-1));
    Borb = Ain2orb*B_I;
    
    %
    JF = stateTransitionMatrix([x_hat(:,i-1); r(:,i-1); V(:,i-1)],dt);
    JH = observationMatrix(x_hat(:,i-1), Borb);
    
    f = getAttitudeVector([x_hat(:,i-1); r(:,i-1); V(:,i-1)], dt);
    h = getHiddenState([x_hat(:,i-1); r(:,i-1); V(:,i-1)], Borb);
    
    fx(:,i) = f;
    hx(:,i) = h;
    
    Jfx(:,i) = JF*x_hat(:,i-1);
    Jhx(:,i) = JH*x_hat(:,i-1);
    
    z = sensor_data(:,i);
    [x_hat(:,i), P(:,:,i)] = kalmanReal(R, Q, z, P(:,:,i-1), JF, JH, h, f);
    normP(i) = norm(P(:,:,i));
end



% err = [abs(( x_hat(1,1:t)-y(1,1:t))./y(1,1:t))*100;
%     abs(( x_hat(2,1:t)-y(2,1:t))./y(2,1:t))*100;
%     abs(( x_hat(3,1:t)-y(3,1:t))./y(3,1:t))*100;
%     abs(( x_hat(4,1:t)-y(4,1:t))./y(4,1:t))*100;
%     abs(( x_hat(5,1:t)-y(5,1:t))./y(5,1:t))*100;
%     abs(( x_hat(6,1:t)-y(6,1:t))./y(6,1:t))*100;
%     abs(( x_hat(7,1:t)-y(7,1:t))./y(7,1:t))*100];

% figure;
% subplot(2,4,1); plot(x(1:t),y(1,1:t), x(1:t), x_hat(1,1:t)); xlabel("Время, с"); ylabel("q0"); legend("данные","фильтр"); grid;
% subplot(2,4,2); plot(x(1:t),y(2,1:t), x(1:t), x_hat(2,1:t)); xlabel("Время, с"); ylabel("q1"); legend("данные","фильтр"); grid;
% subplot(2,4,3); plot(x(1:t),y(3,1:t), x(1:t), x_hat(3,1:t)); xlabel("Время, с"); ylabel("q2"); legend("данные","фильтр"); grid;
% subplot(2,4,4); plot(x(1:t),y(4,1:t), x(1:t), x_hat(4,1:t)); xlabel("Время, с"); ylabel("q3"); legend("данные","фильтр"); grid; 
% subplot(2,4,5); plot(x(1:t), y(5,1:t)*180/pi, x(1:t), x_hat(5,1:t)*180/pi, x(1:t), sensor_data(4,1:t)*180/pi); xlabel("Время, с"); ylabel("w_x,  градусы/c"); grid; legend("данные","фильтр","гироскоп");
% subplot(2,4,6); plot(x(1:t), y(6,1:t)*180/pi, x(1:t), x_hat(6,1:t)*180/pi, x(1:t), sensor_data(5,1:t)*180/pi); xlabel("Время, с"); ylabel("w_y,  градусы/c"); grid; legend("данные","фильтр","гироскоп");
% subplot(2,4,7); plot(x(1:t), y(7,1:t)*180/pi, x(1:t), x_hat(7,1:t)*180/pi, x(1:t), sensor_data(6,1:t)*180/pi); xlabel("Время, с"); ylabel("w_z,  градусы/c"); grid; legend("данные","фильтр","гироскоп");
% subplot(2,4,8); plot(x(1:t),vecnorm(x_hat(1:4,1:t)));legend("норма кватерниона");grid; 

% figure;
% subplot(2,3,1); plot(x(1:t),y(2,1:t), x(1:t), fx(1,1:t)); xlabel("Время, с"); ylabel("q1"); legend("данные","фильтр"); grid;
% subplot(2,3,2); plot(x(1:t),y(3,1:t), x(1:t), fx(2,1:t)); xlabel("Время, с"); ylabel("q2"); legend("данные","фильтр"); grid;
% subplot(2,3,3); plot(x(1:t),y(4,1:t), x(1:t), fx(3,1:t)); xlabel("Время, с"); ylabel("q3"); legend("данные","фильтр"); grid;
% subplot(2,3,4); plot(x(1:t),y(5,1:t), x(1:t), fx(4,1:t)); xlabel("Время, с"); ylabel("wx"); legend("данные","фильтр"); grid;
% subplot(2,3,5); plot(x(1:t),y(6,1:t), x(1:t), fx(5,1:t)); xlabel("Время, с"); ylabel("wy"); legend("данные","фильтр"); grid;
% subplot(2,3,6); plot(x(1:t),y(7,1:t), x(1:t), fx(6,1:t)); xlabel("Время, с"); ylabel("wz"); legend("данные","фильтр"); grid; 
% 
% figure;
% subplot(2,3,1); plot(x(1:t),sensor_data(1,1:t), x(1:t), hx(1,1:t)); xlabel("Время, с"); ylabel("q1"); legend("данные","фильтр"); grid;
% subplot(2,3,2); plot(x(1:t),sensor_data(2,1:t), x(1:t), hx(2,1:t)); xlabel("Время, с"); ylabel("q2"); legend("данные","фильтр"); grid;
% subplot(2,3,3); plot(x(1:t),sensor_data(3,1:t), x(1:t), hx(3,1:t)); xlabel("Время, с"); ylabel("q3"); legend("данные","фильтр"); grid;
% subplot(2,3,4); plot(x(1:t),sensor_data(4,1:t), x(1:t), hx(4,1:t)); xlabel("Время, с"); ylabel("wx"); legend("данные","фильтр"); grid;
% subplot(2,3,5); plot(x(1:t),sensor_data(5,1:t), x(1:t), hx(5,1:t)); xlabel("Время, с"); ylabel("wy"); legend("данные","фильтр"); grid;
% subplot(2,3,6); plot(x(1:t),sensor_data(6,1:t), x(1:t), hx(6,1:t)); xlabel("Время, с"); ylabel("wz"); legend("данные","фильтр"); grid; 

figure;
subplot(2,3,1); plot(x(1:t),y(2,1:t), x(1:t), Jfx(1,1:t)); xlabel("Время, с"); ylabel("q1"); legend("данные","фильтр"); grid;
subplot(2,3,2); plot(x(1:t),y(3,1:t), x(1:t), Jfx(2,1:t)); xlabel("Время, с"); ylabel("q2"); legend("данные","фильтр"); grid;
subplot(2,3,3); plot(x(1:t),y(4,1:t), x(1:t), Jfx(3,1:t)); xlabel("Время, с"); ylabel("q3"); legend("данные","фильтр"); grid;
subplot(2,3,4); plot(x(1:t),y(5,1:t), x(1:t), Jfx(4,1:t)); xlabel("Время, с"); ylabel("wx"); legend("данные","фильтр"); grid;
subplot(2,3,5); plot(x(1:t),y(6,1:t), x(1:t), Jfx(5,1:t)); xlabel("Время, с"); ylabel("wy"); legend("данные","фильтр"); grid;
subplot(2,3,6); plot(x(1:t),y(7,1:t), x(1:t), Jfx(6,1:t)); xlabel("Время, с"); ylabel("wz"); legend("данные","фильтр"); grid; 

figure;
subplot(2,3,1); plot(x(1:t),sensor_data(1,1:t), x(1:t), Jhx(1,1:t)); xlabel("Время, с"); ylabel("q1"); legend("данные","фильтр"); grid;
subplot(2,3,2); plot(x(1:t),sensor_data(2,1:t), x(1:t), Jhx(2,1:t)); xlabel("Время, с"); ylabel("q2"); legend("данные","фильтр"); grid;
subplot(2,3,3); plot(x(1:t),sensor_data(3,1:t), x(1:t), Jhx(3,1:t)); xlabel("Время, с"); ylabel("q3"); legend("данные","фильтр"); grid;
subplot(2,3,4); plot(x(1:t),sensor_data(4,1:t), x(1:t), Jhx(4,1:t)); xlabel("Время, с"); ylabel("wx"); legend("данные","фильтр"); grid;
subplot(2,3,5); plot(x(1:t),sensor_data(5,1:t), x(1:t), Jhx(5,1:t)); xlabel("Время, с"); ylabel("wy"); legend("данные","фильтр"); grid;
subplot(2,3,6); plot(x(1:t),sensor_data(6,1:t), x(1:t), Jhx(6,1:t)); xlabel("Время, с"); ylabel("wz"); legend("данные","фильтр"); grid; 


%figure;
% subplot(2,3,1); plot(x(1:t),y(2,1:t), x(1:t), x_hat(1,1:t)); xlabel("Время, с"); ylabel("q1"); legend("данные","фильтр"); grid;
% subplot(2,3,2); plot(x(1:t),y(3,1:t), x(1:t), x_hat(2,1:t)); xlabel("Время, с"); ylabel("q2"); legend("данные","фильтр"); grid;
% subplot(2,3,3); plot(x(1:t),y(4,1:t), x(1:t), x_hat(3,1:t)); xlabel("Время, с"); ylabel("q3"); legend("данные","фильтр"); grid; 
% subplot(2,3,4); plot(x(1:t), y(5,1:t)*180/pi, x(1:t), x_hat(4,1:t)*180/pi); xlabel("Время, с"); ylabel("w_x,  градусы/c"); grid; legend("данные","фильтр");
% subplot(2,3,5); plot(x(1:t), y(6,1:t)*180/pi, x(1:t), x_hat(5,1:t)*180/pi); xlabel("Время, с"); ylabel("w_y,  градусы/c"); grid; legend("данные","фильтр");
% subplot(2,3,6); plot(x(1:t), y(7,1:t)*180/pi, x(1:t), x_hat(6,1:t)*180/pi); xlabel("Время, с"); ylabel("w_z,  градусы/c"); grid; legend("данные","фильтр");
%subplot(2,4,8); plot(x(1:t),vecnorm(x_hat(1:4,1:t)));legend("норма кватерниона");grid; 

% 
% figure;
% subplot(2,4,1); plot(x(1:t), err(1,1:t)); xlabel("Время, с"); ylabel("q0"); grid;
% subplot(2,4,2); plot(x(1:t), err(2,1:t)); xlabel("Время, с"); ylabel("q1"); grid;
% subplot(2,4,3); plot(x(1:t), err(3,1:t)); xlabel("Время, с"); ylabel("q2"); grid;
% subplot(2,4,4); plot(x(1:t), err(4,1:t)); xlabel("Время, с"); ylabel("q3"); grid; 
% subplot(2,4,5); plot(x(1:t), err(5,1:t)*180/pi); xlabel("Время, с"); ylabel("w_x,  градусы/c"); grid; 
% subplot(2,4,6); plot(x(1:t), err(6,1:t)*180/pi); xlabel("Время, с"); ylabel("w_y,  градусы/c"); grid; 
% subplot(2,4,7); plot(x(1:t), err(7,1:t)*180/pi); xlabel("Время, с"); ylabel("w_z,  градусы/c"); grid;
% subplot(2,4,8); plot(x(1:t),vecnorm(x_hat(1:4,1:t)));legend("норма кватерниона");grid; 


                     
                     
 