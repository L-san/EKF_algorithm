function Jf = jacobiF(attitude,dt)
global satellite orbit_vec
q1 = attitude(1); q2 = attitude(2); q3 = attitude(3); q4 = attitude(4);
wbx = attitude(5); wby = attitude(6); wbz = attitude(7);
I = satellite.I;

r = attitude(8:10);
V = attitude(11:13);
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

dqdot1dq1=1+q1*q3*w0y-2*q2*q4*w0y;
dqdot1dq2=q2*q3*w0y-2*q1*q4*w0y-wbx/2;
dqdot1dq3=1/2*(q1^2*w0y+q2^2*w0y+3*q3^2*w0y+q4^2*w0y-wby);
dqdot1dq4=-2*q1*q2*w0y+q3*q4*w0y-wbz/2;
dqdot1dwbx=-(q2/2);
dqdot1dwby=-(q3/2);
dqdot1dwbz=-(q4/2);

dqdot2dq1=1/2*(6*q1*q4*w0y+wbx);
dqdot2dq2=1-q2*q4*w0y;
dqdot2dq3=1/2*(-2*q3*q4*w0y+wbz);
dqdot2dq4=1/2*(3*q1^2*w0y-q2^2*w0y-q3^2*w0y-3*q4^2*w0y-wby);
dqdot2dwbx=q1/2;
dqdot2dwby=-(q4/2);
dqdot2dwbz=q3/2;

dqdot3dq1=1/2*(-3*q1^2*w0y-q2^2*w0y-q3^2*w0y+3*q4^2*w0y+wby);
dqdot3dq2=-q1*q2*w0y-wbz/2;
dqdot3dq3=1-q1*q3*w0y;
dqdot3dq4=1/2*(6*q1*q4*w0y+wbx);
dqdot3dwbx=q4/2;
dqdot3dwby=q1/2;
dqdot3dwbz=-(q2/2);

dqdot4dq1=q1*q2*w0y-2*q3*q4*w0y+wbz/2;
dqdot4dq2=1/2*(q1^2*w0y+3*q2^2*w0y+q3^2*w0y+q4^2*w0y+wby);
dqdot4dq3=q2*q3*w0y-2*q1*q4*w0y-wbx/2;
dqdot4dq4=1-2*q1*q3*w0y+q2*q4*w0y;
dqdot4dwbx=-(q3/2);
dqdot4dwby=q2/2;
dqdot4dwbz=q1/2;

dwbxdq1=2*dt*(-6*(Iy-Iz)*q1*(q1*q2+q3*q4)*w0^2+3*(Iy-Iz)*q2*(-q1^2+q2^2+q3^2-q4^2)*w0^2+(ka*(q4*Vx-q2*Vy+q1*Vz)*y)/Vn+(ka*(q4*Vx-q1*Vy-q2*Vz)*z)/Vn);
dwbxdq2=dt*(12*(Iy-Iz)*q2*(q1*q2+q3*q4)*w0^2-6*(Iy-Iz)*q1*(q1^2-q2^2-q3^2+q4^2)*w0^2-(2*ka*(q1*Vy*y+q2*Vz*y-q2*Vy*z+q1*Vz*z+q3*Vx*(-y+z)))/Vn);
dwbxdq3=dt*(12*(Iy-Iz)*q3*(q1*q2+q3*q4)*w0^2-6*(Iy-Iz)*q4*(q1^2-q2^2-q3^2+q4^2)*w0^2+(2*ka*(q2*Vx*(y-z)-q3*(Vz*y+Vy*z)+q4*(Vy*y-Vz*z)))/Vn);
dwbxdq4=dt*(-12*(Iy-Iz)*q4*(q1*q2+q3*q4)*w0^2+6*(Iy-Iz)*q3*(-q1^2+q2^2+q3^2-q4^2)*w0^2+(2*ka*(q3*Vy*y+q4*Vz*y+q4*Vy*z-q3*Vz*z+q1*Vx*(y+z)))/Vn);
dwbxdwbx=1;
dwbxdwby=dt*kx*wbz;
dwbxdwbz=dt*kx*wby;

dwbydq1=(1/Vn)*2*dt*(3*Iz*(-3*q1^2*q3+q2^2*q3+q3^3+2*q1*q2*q4-q3*q4^2)*Vn*w0^2+3*Ix*(3*q1^2*q3-2*q1*q2*q4-q3*(q2^2+q3^2-q4^2))*Vn*w0^2+ka*(q1*Vx+q4*Vy-q3*Vz)*(x-z));
dwbydq2=6*dt*Ix*(2*q1*q2*q3+q1^2*q4-3*q2^2*q4-q3^2*q4+q4^3)*w0^2-6*dt*Iz*(2*q1*q2*q3+q1^2*q4-3*q2^2*q4-q3^2*q4+q4^3)*w0^2-(2*dt*ka*(q2*Vx+q3*Vy+q4*Vz)*(x-z))/Vn;
dwbydq3=2*dt*(-3*Ix*(q1^3+2*q2*q3*q4+q1*(-q2^2-3*q3^2+q4^2))*w0^2+3*Iz*(q1^3+2*q2*q3*q4+q1*(-q2^2-3*q3^2+q4^2))*w0^2+(ka*(q3*Vx-q2*Vy+q1*Vz)*(x-z))/Vn);
dwbydq4=2*dt*(6*(Ix-Iz)*q4*(-q1*q3+q2*q4)*w0^2+3*(Ix-Iz)*q2*(q1^2-q2^2-q3^2+q4^2)*w0^2+(ka*(q4*Vx-q1*Vy-q2*Vz)*x)/Vn+(ka*(-q4*Vx+q1*Vy+q2*Vz)*z)/Vn);
dwbydwbx=dt*ky*wbz;
dwbydwby=1;
dwbydwbz=dt*ky*wbx;

dwbzdq1=(1/Vn)*2*dt*(6*Ix*(2*q1*q2*q3+(-q2^2+q3^2)*q4)*Vn*w0^2-6*Iy*(2*q1*q2*q3+(-q2^2+q3^2)*q4)*Vn*w0^2+ka*(q1*Vy*x+q2*Vz*x-q1*Vx*y+q3*Vz*y-q4*(Vx*x+Vy*y)));
dwbzdq2=(1/Vn)*2*dt*(6*(Ix-Iy)*q1*(q1*q3-q2*q4)*Vn*w0^2-6*(Ix-Iy)*q4*(q1*q2+q3*q4)*Vn*w0^2+ka*(q3*Vx-q2*Vy+q1*Vz)*x-ka*(q2*Vx+q3*Vy+q4*Vz)*y);
dwbzdq3=(1/Vn)*2*dt*(6*Ix*(q1^2*q2+2*q1*q3*q4-q2*q4^2)*Vn*w0^2-6*Iy*(q1^2*q2+2*q1*q3*q4-q2*q4^2)*Vn*w0^2+ka*(q2*Vx*x+q3*Vy*x+q4*Vz*x+q3*Vx*y-q2*Vy*y+q1*Vz*y));
dwbzdq4=-(1/Vn)*2*dt*(-6*Iy*(q1*q2^2-q1*q3^2+2*q2*q3*q4)*Vn*w0^2+6*Ix*(q1*(q2^2-q3^2)+2*q2*q3*q4)*Vn*w0^2+ka*(q1*Vx*x+q4*Vy*x-q3*Vz*x-q4*Vx*y+q1*Vy*y+q2*Vz*y));
dwbzdwbx=dt*kz*wby;
dwbzdwby=dt*kz*wbx;
dwbzdwbz=1;

dqdq = [dqdot1dq1 dqdot1dq2 dqdot1dq3 dqdot1dq4;
        dqdot2dq1 dqdot2dq2 dqdot2dq3 dqdot2dq4;
        dqdot3dq1 dqdot3dq2 dqdot3dq3 dqdot3dq4
        dqdot4dq1 dqdot4dq2 dqdot4dq3 dqdot4dq4];

dqdw = [dqdot1dwbx dqdot1dwby dqdot1dwbz;
        dqdot2dwbx dqdot2dwby dqdot2dwbz;
        dqdot3dwbx dqdot3dwby dqdot3dwbz;
        dqdot4dwbx dqdot4dwby dqdot4dwbz];

dwdq = [dwbxdq1 dwbxdq2 dwbxdq3 dwbxdq4;
        dwbydq1 dwbydq2 dwbydq3 dwbydq4;
        dwbzdq1 dwbzdq2 dwbzdq3 dwbzdq4];

dwdw = [dwbxdwbx dwbxdwby dwbxdwbz;
        dwbydwbx dwbydwby dwbydwbz;
        dwbzdwbx dwbzdwby dwbzdwbz];

Jf = [dqdq dqdw;
     dwdq dwdw];
end