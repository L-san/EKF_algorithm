function [dzetaN,errN] = main(filename,Q, pl)
load(filename);
%load("data_distb2.mat");
global mu
        mu = 398600e+9; %gravitational constant
%%
%orbital elemets
%%
dt = 0.01;
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
%загоняем обратные матрицы сразу, чтобы не вычислять их дальше в коде
global M
       M = inv([1.32955081913368,0.0208565423416162,0.0114289793418493;0.0208565423416162,1.22838488999301,-0.00295709623585174;0.0114289793418493,-0.00295709623585174,0.612564210335257]);
global x0
       x0 = [1.90831820362797e-06,1.77462922157422e-06,-3.09982128810498e-06]';
global alphaT
       alphaT = [0.02*pi/180; 0.02*pi/180; 0.02*pi/180];
global A
       A = inv(angle2dcm(0.01*pi/180,0.02*pi/180,0.03*pi/180));
       %%
W = 500;
y =y(:,1+W:end);  
q = y(1:4,:);
w = y(5:7,:);
r = y(8:10,:);
V = y(11:13,:);

sensor_data = [B_b(:,1+W:end); wb(:,1+W:end)];
x_hat = ones(6,length(y))*NaN;

w1 = wb(:,1);
w1 = A\(w1-alphaT);
global q0
       q0 = cos(norm(w1)/2*dt);
x_hat(:,1) = [w1/norm(w1)*sin(norm(w1)/2*dt); w1];
%x_hat(:,1) = [[-0.5;-0.5;0.5]; w1];
R = diag([varb, varw]);
P = diag([varq, varw]);
%Q = eye(6,6)*1e-1;
hi_level = chi2inv(0.95,6);
t = 100001;%floor(length(y));
%t = 1000;
dzetaN = 0;
for i = 2:t
    w = y(2:7,i-1)-x_hat(:,i-1);
    dzeta = w'*inv(P)*w;

    if(dzeta<=hi_level)
        dzetaN = dzetaN+1;
    end

    Ain2orb = in2orb(r(:,i-1),V(:,i-1));
    B_I = magneticField(r(:,i-1));
    Borb = Ain2orb*B_I;
    f = getAttitudeVector([x_hat(:,i-1); r(:,i-1); V(:,i-1)]);
    h = getHiddenState([x_hat(:,i-1); r(:,i-1); V(:,i-1)], Borb);
   
    F = stateTransitionMatrix([x_hat(:,i-1); r(:,i-1); V(:,i-1)], dt);
    H = observationMatrix(x_hat(:,i-1), Borb);
    
    z = sensor_data(:,i);
    [x_hat(:,i), P] = kalmanReal(R, Q, z, P, F, H, h, f);
    x_hat(1:3,i) = x_hat(1:3,i)/norm([q0; x_hat(1:3,i)]);
%     if(max(x_hat(:,i))>1)
%         break
%     end
    %% covariance update   
%     wk = x_hat(:,i)-f;
%     vk = z-h;
%      R = 1/i*((i-1)*R+wk*wk');
%      Q = 1/i*((i-1)*Q+vk*vk');
end
dzetaN = dzetaN/t*100;
prv_angle_f = 2*acos(sqrt(1-sum(x_hat(1:3,1:t).^2)));
prv_angle_y = 2*acos(sqrt(1-sum(y(2:4,1:t).^2)));
errors = [mean(abs(prv_angle_y-prv_angle_f)); mean(abs(x_hat(4:6,1:t)-y(5:7,1:t))')'];

if(isnan(errors))
    errors = ones(6,1);
end

if(pl == 1)
    
    if(max(x_hat(:,end))>10)
        errN = ones(6,1);
    else
        [errN,err] = errcalc(x_hat(:,1:t),y(2:7,1:t));
    end
     
%     figure;
%     subplot(2,2,1); plot(x(1:t), prv_angle_y, x(1:t), prv_angle_f); xlabel("Время, с"); ylabel("Ф, градусы"); legend("данные","фильтр"); grid; 
%     subplot(2,2,2); plot(x(1:t), y(5,1:t)*180/pi, x(1:t), x_hat(4,1:t)*180/pi); xlabel("Время, с"); ylabel("w_x,  градусы/c"); grid; legend("данные","фильтр");
%     subplot(2,2,3); plot(x(1:t), y(6,1:t)*180/pi, x(1:t), x_hat(5,1:t)*180/pi); xlabel("Время, с"); ylabel("w_y,  градусы/c"); grid; legend("данные","фильтр");
%     subplot(2,2,4); plot(x(1:t), y(7,1:t)*180/pi, x(1:t), x_hat(6,1:t)*180/pi); xlabel("Время, с"); ylabel("w_z,  градусы/c"); grid; legend("данные","фильтр");
%     
%     
%     figure;
%     subplot(2,2,1); plot(x(1:t), err(1,1:t)); xlabel("Время, с"); ylabel("Ф, градусы"); grid;
%     subplot(2,2,2); plot(x(1:t), err(2,1:t)); xlabel("Время, с"); ylabel("w_x,  градусы/c"); grid;
%     subplot(2,2,3); plot(x(1:t), err(3,1:t)); xlabel("Время, с"); ylabel("w_y,  градусы/c"); grid;
%     subplot(2,2,4); plot(x(1:t), err(4,1:t)); xlabel("Время, с"); ylabel("w_z,  градусы/c"); grid;
figure;
subplot(2,3,1); plot(x(1:t),y(2,1:t), x(1:t), x_hat(1,1:t)); xlabel("Время, с"); ylabel("q1"); legend("данные","фильтр"); grid;
subplot(2,3,2); plot(x(1:t),y(3,1:t), x(1:t), x_hat(2,1:t)); xlabel("Время, с"); ylabel("q2"); legend("данные","фильтр"); grid;
subplot(2,3,3); plot(x(1:t),y(4,1:t), x(1:t), x_hat(3,1:t)); xlabel("Время, с"); ylabel("q3"); legend("данные","фильтр"); grid; 
subplot(2,3,4); plot(x(1:t), y(5,1:t)*180/pi, x(1:t), x_hat(4,1:t)*180/pi); xlabel("Время, с"); ylabel("w_x,  градусы/c"); grid; legend("данные","фильтр");
subplot(2,3,5); plot(x(1:t), y(6,1:t)*180/pi, x(1:t), x_hat(5,1:t)*180/pi); xlabel("Время, с"); ylabel("w_y,  градусы/c"); grid; legend("данные","фильтр");
subplot(2,3,6); plot(x(1:t), y(7,1:t)*180/pi, x(1:t), x_hat(6,1:t)*180/pi); xlabel("Время, с"); ylabel("w_z,  градусы/c"); grid; legend("данные","фильтр");

figure;
subplot(2,3,1); plot(x(1:t), err(1,1:t)); xlabel("Время, с"); ylabel("q1"); grid;
subplot(2,3,2); plot(x(1:t), err(2,1:t)); xlabel("Время, с"); ylabel("q2"); grid;
subplot(2,3,3); plot(x(1:t), err(3,1:t)); xlabel("Время, с"); ylabel("q3"); grid;
subplot(2,3,4); plot(x(1:t), err(4,1:t)); xlabel("Время, с"); ylabel("w_x,  градусы/c"); grid;
subplot(2,3,5); plot(x(1:t), err(5,1:t)); xlabel("Время, с"); ylabel("w_y,  градусы/c"); grid;
subplot(2,3,6); plot(x(1:t), err(6,1:t)); xlabel("Время, с"); ylabel("w_z,  градусы/c"); grid;
end
end
     
