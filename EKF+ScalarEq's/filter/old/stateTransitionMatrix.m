function Jf = stateTransitionMatrix(attitude,dt)
global satellite orbit_vec q0
q1 = attitude(1); q2 = attitude(2); q3 = attitude(3);
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

dq1dq1=-2*q0*q2*w0y+q1*q3*w0y;
dq1dq2=-2*q0*q1*w0y+q2*q3*w0y-wbz/2;
dq1dq3=1/2*(-3*q0^2*w0y+q1^2*w0y+q2^2*w0y+3*q3^2*w0y+wby);
dq1dw1=q0/2;
dq1dw2=q3/2;
dq1dw3=-(q2/2);

dq2dq1=1/2*(6*q0*q1*w0y+wbz);
dq2dq2=-q0*q2*w0y;
dq2dq3=3*q0*q3*w0y-wbx/2;
dq2dw1=-(q3/2);
dq2dw2=q0/2;
dq2dw3=q1/2;

dq3dq1=1/2*(3*q0^2*w0y-3*q1^2*w0y-q2^2*w0y-q3^2*w0y-wby);
dq3dq2=1/2*(-2*q1*q2*w0y-4*q0*q3*w0y+wbx);
dq3dq3=-((2*q0*q2+q1*q3)*w0y);
dq3dw1=q2/2;
dq3dw2=-(q1/2);
dq3dw3=q0/2;

dwxdq1=-6*kx*(q0^3-2*q1*q2*q3+q0*(-3*q1^2-q2^2+q3^2))*w0^2-(2*ka*(-q3*Vx*y+q0*Vy*y+q1*Vz*y+q2*Vx*z-q1*Vy*z+q0*Vz*z))/(Ix*Vn);
dwxdq2=-6*kx*(-2*q0*q1*q2+q0^2*q3-q1^2*q3-3*q2^2*q3+q3^3)*w0^2-(2*ka*(-q0*Vx*y-q3*Vy*y+q2*Vz*y+q1*Vx*z+q2*Vy*z+q3*Vz*z))/(Ix*Vn);
dwxdq3=-6*kx*(q0^2*q2+2*q0*q1*q3-q2*(q1^2+q2^2-3*q3^2))*w0^2+(2*ka*(q1*Vx*y+q2*Vy*y+q3*Vz*y+q0*Vx*z+q3*Vy*z-q2*Vz*z))/(Ix*Vn);
dwxdw1=0;
dwxdw2=kx*wbz;
dwxdw3=kx*wby;

dwydq1=-6*ky*(2*q0*q1*q2+q0^2*q3-3*q1^2*q3-q2^2*q3+q3^3)*w0^2-(2*ka*(q1*Vx+q2*Vy+q3*Vz)*(x-z))/(Iy*Vn);
dwydq2=6*ky*(q0^3+2*q1*q2*q3+q0*(-q1^2-3*q2^2+q3^2))*w0^2+(2*ka*(q2*Vx-q1*Vy+q0*Vz)*(x-z))/(Iy*Vn);
dwydq3=-6*ky*(q0^2*q1-2*q0*q2*q3-q1*(q1^2+q2^2-3*q3^2))*w0^2+(2*ka*(q3*Vx-q0*Vy-q1*Vz)*(x-z))/(Iy*Vn);
dwydw1=ky*wbz;
dwydw2=0;
dwydw3=ky*wbx;

dwzdq1=12*kz*(q0^2*q2-2*q0*q1*q3-q2*q3^2)*w0^2+(2*ka*(q2*Vx*x-q1*Vy*x+q0*Vz*x-q1*Vx*y-q2*Vy*y-q3*Vz*y))/(Iz*Vn);
dwzdq2=12*kz*(q0^2*q1+2*q0*q2*q3-q1*q3^2)*w0^2+(2*ka*(q1*Vx*x+q2*Vy*x+q3*Vz*x+q2*Vx*y-q1*Vy*y+q0*Vz*y))/(Iz*Vn);
dwzdq3=-12*kz*(q0*(q1^2-q2^2)+2*q1*q2*q3)*w0^2-(2*ka*(q0*Vx*x+q3*Vy*x-q2*Vz*x-q3*Vx*y+q0*Vy*y+q1*Vz*y))/(Iz*Vn);
dwzdw1=kz*wby;
dwzdw2=kz*wbx;
dwzdw3=0;

dqdq = [dq1dq1 dq1dq2 dq1dq3;
	dq2dq1 dq2dq2 dq2dq3;
	dq3dq1 dq3dq2 dq3dq3];

dqdw = [dq1dw1 dq1dw2 dq1dw3;
	dq2dw1 dq2dw2 dq2dw3;
	dq3dw1 dq3dw2 dq3dw3];

dwdq = [dwxdq1 dwxdq2 dwxdq3;
	dwydq1 dwydq2 dwydq3;
	dwzdq1 dwzdq2 dwzdq3];

dwdw = [dwxdw1 dwxdw2 dwxdw3;
	dwydw1 dwydw2 dwydw3;
	dwzdw1 dwzdw2 dwzdw3];


Jf_hat = [dqdq dqdw;
          dwdq dwdw];
  
Jf = eye(6,6)+Jf_hat*dt;
%Jf = Jf_hat*dt;
end