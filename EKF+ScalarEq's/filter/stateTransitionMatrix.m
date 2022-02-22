function Jf = stateTransitionMatrix(attitude,dt)
global satellite orbit_vec
q0 = getQ0(attitude(1:3)); q1 = attitude(1); q2 = attitude(2); q3 = attitude(3);
wbx = attitude(4); wby = attitude(5); wbz = attitude(6);
I = satellite.I;

r = attitude(7:9);
V = attitude(10:12);
ro = Density(norm(r)-orbit_vec.Re);
Vx = V(1); Vy = V(2); Vz = V(3);
x = satellite.ra(1); y = satellite.ra(2); z = satellite.ra(3);
Vn = sqrt(Vx^2+Vy^2+Vz^2);
ka =  -0.5*ro* satellite.c0* satellite.Sm*Vn^2;
w0y = orbit_vec.w0;
w0 = orbit_vec.w0;
Ix = I(1,1); Iy = I(2,2); Iz = I(3,3);
kx=(Iy-Iz)/Ix;
ky=(Iz-Ix)/Iy;
kz=(Ix-Iy)/Iz;

dqdot1dq1=1+q0*q2*w0y;
dqdot1dq2=q1*q2*w0y-wbx/2;
dqdot1dq3=1/2*(q0^2*w0y+q1^2*w0y+3*q2^2*w0y+q3^2*w0y-wby);
dqdot1dq4=q2*q3*w0y-wbz/2;
dqdot1dwbx=-(q1/2);
dqdot1dwby=-(q2/2);
dqdot1dwbz=-(q3/2);

dqdot2dq1=1/2*(-4*q1*q2*w0y-6*q0*q3*w0y+wbx);
dqdot2dq2=1-2*q0*q2*w0y+q1*q3*w0y;
dqdot2dq3=-2*q0*q1*w0y+q2*q3*w0y-wbz/2;
dqdot2dq4=1/2*(-3*q0^2*w0y+q1^2*w0y+q2^2*w0y+3*q3^2*w0y+wby);
dqdot2dwbx=q0/2;
dqdot2dwby=q3/2;
dqdot2dwbz=-(q2/2);

dqdot3dq1=1/2*(-3*q0^2*w0y+3*q1^2*w0y-q2^2*w0y+3*q3^2*w0y+wby);
dqdot3dq2=1/2*(6*q0*q1*w0y+wbz);
dqdot3dq3=1-q0*q2*w0y;
dqdot3dq4=3*q0*q3*w0y-wbx/2;
dqdot3dwbx=-(q3/2);
dqdot3dwby=q0/2;
dqdot3dwbz=q1/2;

dqdot4dq1=1/2*(6*q0*q1*w0y-4*q2*q3*w0y+wbz);
dqdot4dq2=1/2*(3*q0^2*w0y-3*q1^2*w0y-q2^2*w0y-q3^2*w0y-wby);
dqdot4dq3=1/2*(-2*q1*q2*w0y-4*q0*q3*w0y+wbx);
dqdot4dq4=1-2*q0*q2*w0y-q1*q3*w0y;
dqdot4dwbx=q2/2;
dqdot4dwby=-(q1/2);
dqdot4dwbz=q0/2;

dwbxdq1=-6*dt*kx*(3*q0^2*q1+2*q0*q2*q3-q1*(q1^2+q2^2-q3^2))*w0^2+(2*dt*ka*(q2*Vx*y+q0*Vz*y+q3*Vx*z-q0*Vy*z-q1*(Vy*y+Vz*z)))/(Ix*Vn);
dwbxdq2=-6*dt*kx*(q0^3-2*q1*q2*q3+q0*(-3*q1^2-q2^2+q3^2))*w0^2-(2*dt*ka*(-q3*Vx*y+q0*Vy*y+q1*Vz*y+q2*Vx*z-q1*Vy*z+q0*Vz*z))/(Ix*Vn);
dwbxdq3=-6*dt*kx*(-2*q0*q1*q2+q0^2*q3-q1^2*q3-3*q2^2*q3+q3^3)*w0^2-(2*dt*ka*(-q0*Vx*y-q3*Vy*y+q2*Vz*y+q1*Vx*z+q2*Vy*z+q3*Vz*z))/(Ix*Vn);
dwbxdq4=-6*dt*kx*(q0^2*q2+2*q0*q1*q3-q2*(q1^2+q2^2-3*q3^2))*w0^2+(2*dt*ka*(q1*Vx*y+q2*Vy*y+q3*Vz*y+q0*Vx*z+q3*Vy*z-q2*Vz*z))/(Ix*Vn);
dwbxdwbx=1;
dwbxdwby=dt*kx*wbz;
dwbxdwbz=dt*kx*wby;

dwbydq1=6*dt*ky*(3*q0^2*q2-2*q0*q1*q3-q2*(q1^2+q2^2-q3^2))*w0^2-(2*dt*ka*(q0*Vx+q3*Vy-q2*Vz)*(x-z))/(Iy*Vn);
dwbydq2=-6*dt*ky*(2*q0*q1*q2+q0^2*q3-3*q1^2*q3-q2^2*q3+q3^3)*w0^2-(2*dt*ka*(q1*Vx+q2*Vy+q3*Vz)*(x-z))/(Iy*Vn);
dwbydq3=6*dt*ky*(q0^3+2*q1*q2*q3+q0*(-q1^2-3*q2^2+q3^2))*w0^2+(2*dt*ka*(q2*Vx-q1*Vy+q0*Vz)*(x-z))/(Iy*Vn);
dwbydq4=-6*dt*ky*(q0^2*q1-2*q0*q2*q3-q1*(q1^2+q2^2-3*q3^2))*w0^2+(2*dt*ka*(q3*Vx-q0*Vy-q1*Vz)*(x-z))/(Iy*Vn);
dwbydwbx=dt*ky*wbz;
dwbydwby=1;
dwbydwbz=dt*ky*wbx;

dwbzdq1=2*dt*(6*kz*(2*q0*q1*q2+(-q1^2+q2^2)*q3)*w0^2+(ka*(q0*Vy*x+q1*Vz*x-q0*Vx*y+q2*Vz*y-q3*(Vx*x+Vy*y)))/(Iz*Vn));
dwbzdq2=2*dt*(6*kz*(q0^2*q2-2*q0*q1*q3-q2*q3^2)*w0^2+(ka*(q2*Vx*x-q1*Vy*x+q0*Vz*x-q1*Vx*y-q2*Vy*y-q3*Vz*y))/(Iz*Vn));
dwbzdq3=2*dt*(6*kz*(q0^2*q1+2*q0*q2*q3-q1*q3^2)*w0^2+(ka*(q1*Vx*x+q2*Vy*x+q3*Vz*x+q2*Vx*y-q1*Vy*y+q0*Vz*y))/(Iz*Vn));
dwbzdq4=-(1/(Iz*Vn))*2*dt*(6*Iz*kz*(q0*(q1^2-q2^2)+2*q1*q2*q3)*Vn*w0^2+ka*(q0*Vx*x+q3*Vy*x-q2*Vz*x-q3*Vx*y+q0*Vy*y+q1*Vz*y));
dwbzdwbx=dt*kz*wby;
dwbzdwby=dt*kz*wbx;
dwbzdwbz=1;

% dqdq = [dqdot1dq1 dqdot1dq2 dqdot1dq3 dqdot1dq4;
%         dqdot2dq1 dqdot2dq2 dqdot2dq3 dqdot2dq4;
%         dqdot3dq1 dqdot3dq2 dqdot3dq3 dqdot3dq4
%         dqdot4dq1 dqdot4dq2 dqdot4dq3 dqdot4dq4];
% 
% dqdw = [dqdot1dwbx dqdot1dwby dqdot1dwbz;
%         dqdot2dwbx dqdot2dwby dqdot2dwbz;
%         dqdot3dwbx dqdot3dwby dqdot3dwbz;
%         dqdot4dwbx dqdot4dwby dqdot4dwbz];
% 
% dwdq = [dwbxdq1 dwbxdq2 dwbxdq3 dwbxdq4;
%         dwbydq1 dwbydq2 dwbydq3 dwbydq4;
%         dwbzdq1 dwbzdq2 dwbzdq3 dwbzdq4];
% 
% dwdw = [dwbxdwbx dwbxdwby dwbxdwbz;
%         dwbydwbx dwbydwby dwbydwbz;
%         dwbzdwbx dwbzdwby dwbzdwbz];

dqdq = [dqdot2dq2 dqdot2dq3 dqdot2dq4;
        dqdot3dq2 dqdot3dq3 dqdot3dq4
        dqdot4dq2 dqdot4dq3 dqdot4dq4];

dqdw = [dqdot2dwbx dqdot2dwby dqdot2dwbz;
        dqdot3dwbx dqdot3dwby dqdot3dwbz;
        dqdot4dwbx dqdot4dwby dqdot4dwbz];

dwdq = [dwbxdq2 dwbxdq3 dwbxdq4;
        dwbydq2 dwbydq3 dwbydq4;
        dwbzdq2 dwbzdq3 dwbzdq4];

dwdw = [dwbxdwbx dwbxdwby dwbxdwbz;
        dwbydwbx dwbydwby dwbydwbz;
        dwbzdwbx dwbzdwby dwbzdwbz];

Jf = [dqdq dqdw;
     dwdq dwdw];
end