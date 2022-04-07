function F = stateTransitionMatrix(attitude,dt)
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

dq0dq0=q0*q2*w0y;
dq0dq1=1/2*(-2*q0*q3*w0y+2*(q1*q2+q0*q3)*w0y-wbx);
dq0dq2=1/2*(2*q1^2*w0y+2*q2^2*w0y+2*q3^2*w0y+(q0^2-q1^2+q2^2-q3^2)*w0y-wby);
dq0dq3=1/2*(2*q0*q1*w0y+2*(-q0*q1+q2*q3)*w0y-wbz);
dq0dw1=-q1/2;
dq0dw2=-q2/2;
dq0dw3=-q3/2;

dq1dq0=1/2*(-2*q1*q2*w0y-4*q0*q3*w0y-2*(q1*q2+q0*q3)*w0y+wbx);
dq1dq1=1/2*(-4*q0*q2*w0y+2*q1*q3*w0y);
dq1dq2=1/2*(-2*q0*q1*w0y+2*(-q0*q1+q2*q3)*w0y-wbz);
dq1dq3=1/2*(-2*q0^2*w0y+2*q2^2*w0y+2*q3^2*w0y-(q0^2-q1^2+q2^2-q3^2)*w0y+wby);
dq1dw1=q0/2;
dq1dw2=q3/2;
dq1dw3=-q2/2;

dq2dq0=1/2*(-2*q0^2*w0y+2*q1^2*w0y+2*q3^2*w0y-(q0^2-q1^2+q2^2-q3^2)*w0y+wby);
dq2dq1=1/2*(4*q0*q1*w0y+2*q2*q3*w0y-2*(-q0*q1+q2*q3)*w0y+wbz);
dq2dq2=-q0*q2*w0y;
dq2dq3=1/2*(-2*q1*q2*w0y+4*q0*q3*w0y+2*(q1*q2+q0*q3)*w0y-wbx);
dq2dw1=-q3/2;
dq2dw2=q0/2;
dq2dw3=q1/2;

dq3dq0=1/2*(4*q0*q1*w0y-2*q2*q3*w0y-2*(-q0*q1+q2*q3)*w0y+wbz);
dq3dq1=1/2*(2*q0^2*w0y-2*q1^2*w0y-2*q2^2*w0y+(q0^2-q1^2+q2^2-q3^2)*w0y-wby);
dq3dq2=1/2*(-2*q0*q3*w0y-2*(q1*q2+q0*q3)*w0y+wbx);
dq3dq3=1/2*(-4*q0*q2*w0y-2*q1*q3*w0y);
dq3dw1=q2/2;
dq3dw2=-q1/2;
dq3dw3=q0/2;

dwxdq0=-12*kx*q0*(q0*q1+q2*q3)*w0^2-6*kx*q1*(q0^2-q1^2-q2^2+q3^2)*w0^2+(ka*(((2*q2*Vx-2*q1*Vy+2*q0*Vz)*y)/Vn-((-2*q3*Vx+2*q0*Vy+2*q1*Vz)*z)/Vn))/Ix;
dwxdq1=12*kx*q1*(q0*q1+q2*q3)*w0^2-6*kx*q0*(q0^2-q1^2-q2^2+q3^2)*w0^2+(ka*(((2*q3*Vx-2*q0*Vy-2*q1*Vz)*y)/Vn-((2*q2*Vx-2*q1*Vy+2*q0*Vz)*z)/Vn))/Ix;
dwxdq2=12*kx*q2*(q0*q1+q2*q3)*w0^2-6*kx*q3*(q0^2-q1^2-q2^2+q3^2)*w0^2+(ka*(((2*q0*Vx+2*q3*Vy-2*q2*Vz)*y)/Vn-((2*q1*Vx+2*q2*Vy+2*q3*Vz)*z)/Vn))/Ix;
dwxdq3=-12*kx*q3*(q0*q1+q2*q3)*w0^2-6*kx*q2*(q0^2-q1^2-q2^2+q3^2)*w0^2+(ka*(((2*q1*Vx+2*q2*Vy+2*q3*Vz)*y)/Vn-((-2*q0*Vx-2*q3*Vy+2*q2*Vz)*z)/Vn))/Ix;
dwxdw1=0;
dwxdw2=kx*wbz;
dwxdw3=kx*wby;

dwydq0=-12*ky*q0*(-q0*q2+q1*q3)*w0^2+6*ky*q2*(q0^2-q1^2-q2^2+q3^2)*w0^2+(ka*(-(((2*q0*Vx+2*q3*Vy-2*q2*Vz)*x)/Vn)+((2*q0*Vx+2*q3*Vy-2*q2*Vz)*z)/Vn))/Iy;
dwydq1=12*ky*q1*(-q0*q2+q1*q3)*w0^2-6*ky*q3*(q0^2-q1^2-q2^2+q3^2)*w0^2+(ka*(-(((2*q1*Vx+2*q2*Vy+2*q3*Vz)*x)/Vn)+((2*q1*Vx+2*q2*Vy+2*q3*Vz)*z)/Vn))/Iy;
dwydq2=12*ky*q2*(-q0*q2+q1*q3)*w0^2+6*ky*q0*(q0^2-q1^2-q2^2+q3^2)*w0^2+(ka*(-(((-2*q2*Vx+2*q1*Vy-2*q0*Vz)*x)/Vn)+((-2*q2*Vx+2*q1*Vy-2*q0*Vz)*z)/Vn))/Iy;
dwydq3=-12*ky*q3*(-q0*q2+q1*q3)*w0^2-6*ky*q1*(q0^2-q1^2-q2^2+q3^2)*w0^2+(ka*(-(((-2*q3*Vx+2*q0*Vy+2*q1*Vz)*x)/Vn)+((-2*q3*Vx+2*q0*Vy+2*q1*Vz)*z)/Vn))/Iy;
dwydw1=ky*wbz;
dwydw2=0;
dwydw3=ky*wbx;

dwzdq0=-12*kz*q1*(-q0*q2+q1*q3)*w0^2+12*kz*q2*(q0*q1+q2*q3)*w0^2+(ka*(((-2*q3*Vx+2*q0*Vy+2*q1*Vz)*x)/Vn-((2*q0*Vx+2*q3*Vy-2*q2*Vz)*y)/Vn))/Iz;
dwzdq1=-12*kz*q0*(-q0*q2+q1*q3)*w0^2-12*kz*q3*(q0*q1+q2*q3)*w0^2+(ka*(((2*q2*Vx-2*q1*Vy+2*q0*Vz)*x)/Vn-((2*q1*Vx+2*q2*Vy+2*q3*Vz)*y)/Vn))/Iz;
dwzdq2=-12*kz*q3*(-q0*q2+q1*q3)*w0^2+12*kz*q0*(q0*q1+q2*q3)*w0^2+(ka*(((2*q1*Vx+2*q2*Vy+2*q3*Vz)*x)/Vn-((-2*q2*Vx+2*q1*Vy-2*q0*Vz)*y)/Vn))/Iz;
dwzdq3=-12*kz*q2*(-q0*q2+q1*q3)*w0^2-12*kz*q1*(q0*q1+q2*q3)*w0^2+(ka*(((-2*q0*Vx-2*q3*Vy+2*q2*Vz)*x)/Vn-((-2*q3*Vx+2*q0*Vy+2*q1*Vz)*y)/Vn))/Iz;
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
  
F = eye(6,6)+Jf_hat*dt;
%Jf = Jf_hat*dt;
end