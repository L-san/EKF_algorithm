clear variables; clc;
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
statvec = kepel_statvec([orbit_vec.sma;
                         orbit_vec.ecc;
                         orbit_vec.inc;
                         orbit_vec.raan;
                         orbit_vec.aop;
                         orbit_vec.ta]);

M = [1.32955081913368,0.0208565423416162,0.0114289793418493;0.0208565423416162,1.22838488999301,-0.00295709623585174;0.0114289793418493,-0.00295709623585174,0.612564210335257];
xshift = [1.90831820362797e-06,1.77462922157422e-06,-3.09982128810498e-06]';
%alphaT = [0.022032687411326; -0.080079545805854; -0.004240626483571];
alphaT = 0.02*pi/180;
A = angle2dcm(0.01*pi/180,0.02*pi/180,0.03*pi/180);
%integrator
h = 0.01;%step
N = 1/15;
T = 2*pi*sqrt(orbit_vec.sma^3/mu);
t = 0:h:N*T;
x = zeros(length(t));
y = zeros(13,length(t));
%%initial conditions
invel = [0, 0, 0.2*pi];
q0 = 1; q = 1/2*q0*invel;
y0 = [q0 q invel statvec]'; x0 = 0;
y(:,1) = y0; x(1) = x0;
%%integrating...
%%kalmaaaan F

[x,y] = ode4(@motionEquations,h,[0,N*T],[x0;y0]);
for i = 1:length(t)
    [a(i) b(i) c(i)] = quat2angle(y(1:4,i)','XYX');
    q_norm(i) = sqrt(y(1,i)^2+y(2,i)^2+y(3,i)^2+y(4,i)^2);
end
sigmaw = sqrt(0.05);
for j = 1:length(y)
    V = y(11:13,j);
    r = y(8:10,j);
    Ain2orb = in2orb(r,V);
    B_I = magneticField(r);
    Borb = Ain2orb*B_I;
    DCM = quat2DCM(y(1:4,j)');
    %B_b(:,j) = quat_mult(quat_mult(y(1:4,j),[0; Borb]),quat_conj(y(1:4,j)));
    B_b(:,j) = DCM*Borb;
end
mes = awgn([B_b; y(5:7,:)],151,'measured');
B_b = mes(1:3,:);
wb = mes(4:6,:);
nana = mean(abs(wb(1,2:end) - y(5,2:end))./abs(y(5,2:end))*100)
B_b = inv(M)*B_b+xshift;
wb = inv(A)*wb+alphaT;

% [pitch, roll, yaw] = quat2angle(y(1:4,:)');
% subplot(1,3,1); plot(x,pitch);
% subplot(1,3,2); plot(x,roll);
% subplot(1,3,3); plot(x,yaw);
%--------------------------------------------------------------------------

% figure;
% subplot(1,2,1); plot(x,a);xlabel("¬рем€, с"); ylabel("угол прецессии, с^-^1");grid;
% subplot(1,2,2); plot(x,b);xlabel("¬рем€, с"); ylabel("угол нутации, с^-^1");grid;
% 
% figure;
% subplot(1,2,1); plot(x,c);xlabel("Time,s"); ylabel("угол соб.вращени€, с^-^1");grid;
% subplot(1,2,2); plot(x, y(5,:)); xlabel("Time,s"); ylabel("w_x, рад/c"); grid;
% 
% figure;
% subplot(1,2,1); plot(x, y(6,:)); xlabel("Time,s"); ylabel("w_y, рад/c"); grid;
% subplot(1,2,2); plot(x, y(7,:)); xlabel("Time,s"); ylabel("w_z, рад/c"); grid;
% 
% figure; plot(x,q_norm); xlabel("Time,s"); ylabel("quaternion norm"); grid;
% 
% figure;
% plot3(y(8,:),y(9,:),y(10,:)); grid;
% 


% figure;
% subplot(1,3,1); plot(x,B_b(1,:)*1e6);xlabel("¬рем€, с"); ylabel("B_x, мк“л");grid;
% subplot(1,3,2); plot(x,B_b(2,:)*1e6);xlabel("¬рем€, с"); ylabel("B_y, мк“л");grid;
% subplot(1,3,3); plot(x,B_b(3,:)*1e6);xlabel("¬рем€, с"); ylabel("B_z, мк“л");grid;
% figure;
subplot(1,3,1); plot(x, wb(1,:)*180/pi, x, y(5,:)*180/pi); xlabel("¬рем€, с"); ylabel("w_x,  градусы/c"); grid;
subplot(1,3,2); plot(x, wb(2,:)*180/pi, x, y(6,:)*180/pi); xlabel("¬рем€, с"); ylabel("w_y,  градусы/c"); grid;
subplot(1,3,3); plot(x, wb(3,:)*180/pi, x, y(7,:)*180/pi); xlabel("¬рем€, с"); ylabel("w_z,  градусы/c"); grid;

save('data.mat','x','y','h','B_b','wb');