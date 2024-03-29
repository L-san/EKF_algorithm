function kl = main(filename)
%clear variables; clc;
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

SNR = 20*log10(1/0.05);
SNR_GNSS = 20*log10(orbit_vec.sma/10);

M = [1.32955081913368,0.0208565423416162,0.0114289793418493;0.0208565423416162,1.22838488999301,-0.00295709623585174;0.0114289793418493,-0.00295709623585174,0.612564210335257];
xshift = [1.90831820362797e-06,1.77462922157422e-06,-3.09982128810498e-06]';
%alphaT = [0.022032687411326; -0.080079545805854; -0.004240626483571];
alphaT = [0.02*pi/180; 0.02*pi/180; 0.02*pi/180];
A = angle2dcm(0.01*pi/180,0.02*pi/180,0.03*pi/180);
%integrator
h = 0.01;%step
N = 1/5;
T = 2*pi*sqrt(orbit_vec.sma^3/mu);
%T = 200;
t = 0:h:N*T;
x = zeros(1,length(t));
y = zeros(13,length(t));
B_b = zeros(3,length(t));
%%initial conditions
invel = [0, 0, 0.02*pi];
y0 = [[1 0 0 0] invel statvec]'; x0 = 0;
y(:,1) = y0; x(1) = x0;
%%integrating...
[x,y] = ode4(@motionEquations,h,[0,N*T],[x0;y0]);

for j = 1:length(y)
    V = y(11:13,j);
    r = y(8:10,j);
    Ain2orb = in2orb(r,V);
    B_I = magneticField(r);
    Borb = Ain2orb*B_I;
    DCM = quat2DCM(y(1:4,j)');
    %B_b(:,j) = quat_mult(quat_mult(y(1:4,j),[0; Borb]),quat_conj(y(1:4,j)));
    %B_b(:,j) = quat_mult(quat_mult(y(1:4,j),[0; Borb]),quat_conj(y(1:4,j)));
    B_b(:,j) = DCM*Borb;
end


MX = 0.05;
DX = 0.05^2;

B_b = inv(M)*(B_b)+xshift;
wb = inv(A)*(y(5:7,:))+alphaT;

dt = 0.01;
B_b1 = B_b;
wb1 = y(5:7,:);
% bx = normrnd(MX*mean(B_b(1,:)),DX*mean(B_b(1,:)),[1,length(t)]);
% by = normrnd(MX*mean(B_b(2,:)),DX*mean(B_b(2,:)),[1,length(t)]);
% bz = normrnd(MX*mean(B_b(3,:)),DX*mean(B_b(3,:)),[1,length(t)]);
% B_b = B_b+[bx;by;bz];% awgn(B_b,30,'measured');
% 
% wx = normrnd(MX*mean(wb(1,:)),DX*mean(wb(1,:)),[1,length(t)]);
% wy = normrnd(MX*mean(wb(2,:)),DX*mean(wb(2,:)),[1,length(t)]);
% wz = normrnd(MX*mean(wb(3,:)),DX*mean(wb(3,:)),[1,length(t)]);
% wb = wb+[wx; wy; wz]; %awgn(wb,30,'measured');
% %%
% r1 = y(8:10,:);
% r2 = y(11:13,:);
% yx = normrnd(MX*mean(y(8,:)),DX*mean(y(8,:)),[1,length(t)]);
% yy = normrnd(MX*mean(y(9,:)),DX*mean(y(9,:)),[1,length(t)]);
% yz = normrnd(MX*mean(y(10,:)),DX*mean(y(10,:)),[1,length(t)]);
% y(8:10,:) = y(8:10,:)+[yz; yy; yz]; %awgn(y(8:10,:),50,'measured');
% 
% vx = normrnd(MX*mean(y(11,:)),DX*mean(y(11,:)),[1,length(t)]);
% vy = normrnd(MX*mean(y(12,:)),DX*mean(y(12,:)),[1,length(t)]);
% vz = normrnd(MX*mean(y(13,:)),DX*mean(y(13,:)),[1,length(t)]);
% y(11:13,:) = y(11:13,:)+[vx; vy; vz]; %awgn(y(11:13,:),50,'measured');
y1 = y(8:10,:);
y2 = y(11:13,:);
B_b = awgn(B_b,SNR,'measured');
wb = awgn(wb,SNR,'measured');
y(8:10,:) = awgn(y(8:10,:),SNR_GNSS,'measured');
y(11:13,:) = awgn(y(11:13,:),SNR_GNSS,'measured');
%%
qw1 = wb1/norm(wb1).*sin(norm(wb1)/2*dt);
qw = wb/norm(wb).*sin(norm(wb)/2*dt);

varb = var((B_b-B_b1)');
varw = var((wb-wb1)');
varq = var((qw-qw1)');

meanb = mean((B_b-B_b1)');
meanw = mean((wb-wb1)');
meanq = mean((qw-qw1)');

% [pitch, roll, yaw] = quat2angle(y(1:4,:)');
% subplot(1,3,1); plot(x,pitch);
% subplot(1,3,2); plot(x,roll);
% subplot(1,3,3); plot(x,yaw);
%--------------------------------------------------------------------------

% figure;
% subplot(1,2,1); plot(x,a);xlabel("�����, �"); ylabel("���� ���������, �^-^1");grid;
% subplot(1,2,2); plot(x,b);xlabel("�����, �"); ylabel("���� �������, �^-^1");grid;
% 
% figure;
% subplot(1,2,1); plot(x,c);xlabel("Time,s"); ylabel("���� ���.��������, �^-^1");grid;
% subplot(1,2,2); plot(x, y(5,:)); xlabel("Time,s"); ylabel("w_x, ���/c"); grid;
% 
% figure;
% subplot(1,2,1); plot(x, y(6,:)); xlabel("Time,s"); ylabel("w_y, ���/c"); grid;
% subplot(1,2,2); plot(x, y(7,:)); xlabel("Time,s"); ylabel("w_z, ���/c"); grid;
% 
% figure; plot(x,q_norm); xlabel("Time,s"); ylabel("quaternion norm"); grid;
% 
% figure;
% plot3(y(8,:),y(9,:),y(10,:)); grid;
% 


figure;
subplot(4,3,1); plot(x,B_b(1,:)*1e6);xlabel("�����, �"); ylabel("B_x, ����");grid;
subplot(4,3,2); plot(x,B_b(2,:)*1e6);xlabel("�����, �"); ylabel("B_y, ����");grid;
subplot(4,3,3); plot(x,B_b(3,:)*1e6);xlabel("�����, �"); ylabel("B_z, ����");grid;
subplot(4,3,4); plot(x, wb(1,:)*180/pi, x, y(5,:)*180/pi); xlabel("�����, �"); ylabel("w_x,  �������/c"); grid;
subplot(4,3,5); plot(x, wb(2,:)*180/pi, x, y(6,:)*180/pi); xlabel("�����, �"); ylabel("w_y,  �������/c"); grid;
subplot(4,3,6); plot(x, wb(3,:)*180/pi, x, y(7,:)*180/pi); xlabel("�����, �"); ylabel("w_z,  �������/c"); grid;
subplot(4,3,7); plot(x, y(8,:), x, y1(1,:));
subplot(4,3,8); plot(x, y(9,:), x, y1(2,:));
subplot(4,3,9); plot(x, y(10,:), x, y1(3,:));
subplot(4,3,10); plot(x, y(11,:), x, y2(1,:));
subplot(4,3,11); plot(x, y(12,:), x,  y2(2,:));
subplot(4,3,12); plot(x, y(13,:), x, y2(3,:));
save(filename,'x','y','h','B_b','wb','varb','varw','varq','meanb','meanw','meanq');
kl = 0;
end