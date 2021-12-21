function F = stateTransitionMatrix(attitude,dt)
global satellite orbit_vec 
q0 = attitude(1); q1 = attitude(2); q2 = attitude(3); q3 = attitude(4);
% wx = attitude(5); wy = attitude(6); wz = attitude(7);
wbx = attitude(5); wby = attitude(6); wbz = attitude(7);
I = satellite.I;

r = attitude(8:10);
V = attitude(11:13);
ro = Density(norm(r)-orbit_vec.Re);
Vx = V(1); Vy = V(2); Vz = V(3);
x = satellite.ra(1); y = satellite.ra(2); z = satellite.ra(3);
Vn = sqrt(Vx^2+Vy^2+Vz^2);
ka =  -0.5*ro* satellite.c0* satellite.Sm*Vn^2;
w0 = orbit_vec.w0;
Ix = I(1,1); Iy = I(2,2); Iz = I(3,3);
kx=(Iy-Iz)/Ix;
ky=(Iz-Ix)/Iy;
kz=(Ix-Iy)/Iz;
%
% kx=dt*(Iy-Iz)/Ix;
% ky=dt*(Iz-Ix)/Iy;
% kz=dt*(Ix-Iy)/Iz;
% 
% w11=kx/2*(-q1*w0*((q0^2-q1^2+q2^2-q3^2)*w0+wy)+q0*w0*((-2*q0*q1+2*q2*q3)*w0+wz));
% w12=kx/2*(-q0*w0*((q0^2-q1^2+q2^2-q3^2)*w0+wy)-q1*w0*((-2*q0*q1+2*q2*q3)*w0+wz));
% w13=kx/2*(q3*w0*((q0^2-q1^2+q2^2-q3^2)*w0+wy)+q2*w0*((-2*q0*q1+2*q2*q3)*w0+wz));
% w14=kx/2*(q2*w0*((q0^2-q1^2+q2^2-q3^2)*w0+wy)-q3*w0*((-2*q0*q1+2*q2*q3)*w0+wz));
% w15=1;
% w16=kx/2*((-2*q0*q1+2*q2*q3)*w0+wz);
% w17=kx/2*((q0^2-q1^2+q2^2-q3^2)*w0+wy);
% 
% w21=ky*(q1*w0*((2*q1*q2-2*q0*q3)*w0+wx)-q3*w0*((-2*q0*q1+1/2*q2*q3)*w0+1/2*wz));
% w22=ky*(-q0*w0*((2*q1*q2-2*q0*q3)*w0+wx)+q2*w0*((-2*q0*q1+1/2*q2*q3)*w0+2/3*wz));
% w23=ky*(q3*w0*((2*q1*q2-q0*q3)*w0+wx)+q1*w0*((-2*q0*q1+1/2*q2*q3)*w0+2/3*wz))*q2;
% w24=ky*(q2*w0*((1*q1*q2-2*q0*q3)*w0+2*wx)-q0*w0*((-2*q0*q1+1/2*q2*q3)*w0+1/2*wz));
% w25=ky*((-2*q0*q1-q2*q3)*w0+1/2*wz);
% w26=1;
% w27=ky*((2/3*q1*q2-q0*q3)*w0+1/2*wx);
% 
% w31=kz/2*(q0*w0*((2*q1*q2-2*q0*q3)*w0+wx)-q3*w0*((q0^2-q1^2+q2^2-q3^2)*w0+wy));
% w32=kz/2*(-q1*w0*((2*q1*q2-2*q0*q3)*w0+wx)+q2*w0*((q0^2-q1^2+q2^2-q3^2)*w0+wy));
% w33=kz/2*(q2*w0*((2*q1*q2-2*q0*q3)*w0+wx)+q1*w0*((q0^2-q1^2+q2^2-q3^2)*w0+wy));
% w34=kz/2*(-q3*w0*((2*q1*q2-2*q0*q3)*w0+wx)-q0*w0*((q0^2-q1^2+q2^2-q3^2)*w0+wy));
% w35=kz/2*((q0^2-q1^2+q2^2-q3^2)*w0+wy);
% w36=kz/2*((2*q1*q2-2*q0*q3)*w0+wx);
% w37=1;

% 
% F1 = [1, -(1/4)*wx*dt, -(1/4)*wy*dt, -(1/4)*wz*dt, -(1/4)*q1*dt, -(1/4)*q2*dt,-(1/4)*q3*dt;
%     (1/4)*wx*dt, 1, (1/4)*wz*dt, -(1/4)*wy*dt, (1/4)*q0*dt, -(1/4)*q3*dt, (1/4)*q2*dt;
%     (1/4)*wy*dt, -(1/4)*wz*dt, 1, (1/4)*wx*dt, (1/4)*q3*dt, (1/4)*q0*dt, -(1/4)*q1*dt;
%     (1/4)*wz*dt, (1/4)*wy*dt, -(1/4)*wx*dt, 1 -(1/4)*q2*dt, (1/4)*q1*dt, (1/4)*q0*dt;
%      0, 0, 0, 0, 1, (1/2)*(Iy-Iz)*wz/Ix*dt, (1/2)*(Iy-Iz)*wy/Ix*dt;
%      0, 0, 0, 0, (1/2)*(Iz-Ix)*wz/Iy*dt, 1, (1/2)*(Iz-Ix)*wx/Iy*dt;
%      0, 0, 0, 0, (1/2)*(Ix-Iy)*wy/Iz*dt, (1/2)*(Ix-Iy)*wx/Iz*dt, 1];

% F1 = [1, -(1/4)*wx*dt, -(1/4)*wy*dt, -(1/4)*wz*dt, -(1/4)*q1*dt, -(1/4)*q2*dt,-(1/4)*q3*dt;
%     (1/4)*wx*dt, 1, (1/4)*wz*dt, -(1/4)*wy*dt, (1/4)*q0*dt, -(1/4)*q3*dt, (1/4)*q2*dt;
%     (1/4)*wy*dt, -(1/4)*wz*dt, 1, (1/4)*wx*dt, (1/4)*q3*dt, (1/4)*q0*dt, -(1/4)*q1*dt;
%     (1/4)*wz*dt, (1/4)*wy*dt, -(1/4)*wx*dt, 1 -(1/4)*q2*dt, (1/4)*q1*dt, (1/4)*q0*dt;
%      w11, w12, w13, w14, w15, w16, w17;
%      w21, w22, w23, w24, w25, w26, w27;
%      w31, w32, w33, w34, w35, w36, w37];


Ma11 = ka*(((q2*Vx-q1*Vy+q0*Vz)*y)-((-q3*Vx+q0*Vy+q1*Vz)*z))/Vn;
Ma12 = ka*(((q3*Vx-q0*Vy-q1*Vz)*y)-((q2*Vx-q1*Vy+q0*Vz)*z))/Vn;
Ma13 = ka*(((q0*Vx+q3*Vy-q2*Vz)*y)-((q1*Vx+q2*Vy+q3*Vz)*z))/Vn;
Ma14 = ka*(((q1*Vx+q2*Vy+q3*Vz)*y)-((-q0*Vx-q3*Vy+q2*Vz)*z))/Vn;
Ma15 = 0;
Ma16 = 0;
Ma17 = 0;

Ma21 = ka*(-(((q0*Vx-q3*Vy-q2*Vz)*x))+((q0*Vx-q3*Vy-q2*Vz)*z))/Vn;
Ma22 = ka*(-(((q1*Vx+q2*Vy+q3*Vz)*x))+((q1*Vx+q2*Vy+q3*Vz)*z))/Vn;
Ma23 = ka*(-(((-q2*Vx+q1*Vy-q0*Vz)*x))+((-q2*Vx+q1*Vy-q0*Vz)*z))/Vn;
Ma24 = ka*(-(((-q3*Vx-q0*Vy+q1*Vz)*x))+((-q3*Vx-q0*Vy+q1*Vz)*z))/Vn;
Ma25 = 0;
Ma26 = 0;
Ma27 = 0;

Ma31 = ka*(((-q3*Vx+q0*Vy+q1*Vz)*x)-((q0*Vx-q3*Vy-q2*Vz)*y))/Vn;
Ma32 = ka*(((q2*Vx-q1*Vy+q0*Vz)*x)-((q1*Vx+q2*Vy+q3*Vz)*y))/Vn;
Ma33 = ka*(((q1*Vx+q2*Vy+q3*Vz)*x)-((-q2*Vx+q1*Vy-q0*Vz)*y))/Vn;
Ma34 = ka*(((-q0*Vx-q3*Vy+q2*Vz)*x)-((-q3*Vx-q0*Vy+q1*Vz)*y))/Vn;
Ma35 = 0;
Ma36 = 0;
Ma37 = 0;

Ma = [Ma11/Ix Ma12/Ix Ma13/Ix Ma14/Ix Ma15 Ma16 Ma17;
     Ma21/Iy Ma22/Iy Ma23/Iy Ma24/Iy Ma25 Ma26 Ma27;
     Ma31/Iz Ma32/Iz Ma33/Iz Ma34/Iz Ma35 Ma36 Ma37];
 
Mg11 = 6*(-Iy+Iz)*q0*(q0*q1+q2*q3)*w0^2+3*(-Iy+Iz)*q1*(q0^2-q1^2-q2^2+q3^2)*w0^2;
Mg12 = -6*(-Iy+Iz)*q1*(q0*q1+q2*q3)*w0^2+3*(-Iy+Iz)*q0*(q0^2-q1^2-q2^2+q3^2)*w0^2;
Mg13 = -6*(-Iy+Iz)*q2*(q0*q1+q2*q3)*w0^2+3*(-Iy+Iz)*q3*(q0^2-q1^2-q2^2+q3^2)*w0^2;
Mg14 = 6*(-Iy+Iz)*q3*(q0*q1+q2*q3)*w0^2+3*(-Iy+Iz)*q2*(q0^2-q1^2-q2^2+q3^2)*w0^2;
Mg15 = 0;
Mg16 = 0;
Mg17 = 0;

Mg21 = 6*(Ix-Iz)*q0*(-q0*q2+q1*q3)*w0^2-3*(Ix-Iz)*q2*(q0^2-q1^2-q2^2+q3^2)*w0^2;
Mg22 = -6*(Ix-Iz)*q1*(-q0*q2+q1*q3)*w0^2+3*(Ix-Iz)*q3*(q0^2-q1^2-q2^2+q3^2)*w0^2;
Mg23 = -6*(Ix-Iz)*q2*(-q0*q2+q1*q3)*w0^2-3*(Ix-Iz)*q0*(q0^2-q1^2-q2^2+q3^2)*w0^2;
Mg24 = 6*(Ix-Iz)*q3*(-q0*q2+q1*q3)*w0^2+3*(Ix-Iz)*q1*(q0^2-q1^2-q2^2+q3^2)*w0^2;
Mg25 = 0;
Mg26 = 0;
Mg27 = 0;

Mg31 = 6*(-Ix+Iy)*q1*(-q0*q2+q1*q3)*w0^2-6*(-Ix+Iy)*q2*(q0*q1+q2*q3)*w0^2;
Mg32 = 6*(-Ix+Iy)*q0*(-q0*q2+q1*q3)*w0^2+6*(-Ix+Iy)*q3*(q0*q1+q2*q3)*w0^2;
Mg33 = 6*(-Ix+Iy)*q3*(-q0*q2+q1*q3)*w0^2-6*(-Ix+Iy)*q0*(q0*q1+q2*q3)*w0^2;
Mg34 = 6*(-Ix+Iy)*q2*(-q0*q2+q1*q3)*w0^2+6*(-Ix+Iy)*q1*(q0*q1+q2*q3)*w0^2;
Mg35 = 0;
Mg36 = 0;
Mg37 = 0;

Mg = [Mg11/Ix Mg12/Ix Mg13/Ix Mg14/Ix Mg15 Mg16 Mg17;
     Mg21/Iy Mg22/Iy Mg23/Iy Mg24/Iy Mg25 Mg26 Mg27;
     Mg31/Iz Mg32/Iz Mg33/Iz Mg34/Iz Mg35 Mg36 Mg37];

M = (Mg+Ma)*dt;

A11=1+1/4*(2*q0*q2*w0-4*q1*q3*w0);
A12=1/4*(-2*q0*q3*w0+(2*q1*q2-2*q0*q3)*w0-wbx);
A13=1/4*(2*q1^2*w0+2*q2^2*w0+2*q3^2*w0+(q0^2-q1^2+q2^2-q3^2)*w0-wby);
A14=1/2*(-2*q0*q1*w0+(-2*q0*q1+2*q2*q3)*w0-wbz);
A15=-q1/4;
A16=-q2/4;
A17=-q3/4;

A21=1/4*(2*q1*q2*w0+4*q0*q3*w0-(2*q1*q2-2*q0*q3)*w0+wbx);
A22=(1-1/2*q1*q3*w0);
A23=1/4*(-2*q0*q1*w0-(-2*q0*q1+2*q2*q3)*w0+wbz);
A24=1/4*(2*q0^2*w0-2*q2^2*w0-2*q3^2*w0+(q0^2-q1^2+q2^2-q3^2)*w0-wby);
A25=q0/4;
A26=-q3/4;
A27=q2/4;

A31=1/4*(-2*q0^2*w0-2*q1^2*w0+2*q3^2*w0-(q0^2-q1^2+q2^2-q3^2)*w0+wby);
A32=1/4*(-2*q2*q3*w0+(-2*q0*q1+2*q2*q3)*w0-wbz);
A33=1-1/2*q0*q2*w0;
A34=1/4*(2*q1*q2*w0+4*q0*q3*w0-(2*q1*q2-2*q0*q3)*w0+wbx);
A35=q3/4;
A36=q0/4;
A37=-q1/4;

A41=1/4*(-2*q2*q3*w0-(-2*q0*q1+2*q2*q3)*w0+wbz);
A42=1/4*(2*q0^2*w0+2*q1^2*w0+2*q2^2*w0-(q0^2-q1^2+q2^2-q3^2)*w0+wby);
A43=1/4*(-2*q0*q3*w0+(2*q1*q2-2*q0*q3)*w0-wbx);
A44=1+1/4*(-4*q0*q2*w0+2*q1*q3*w0);
A45=-q2/4;
A46=q1/4;
A47=q0/4;

B11=M(1,1);
B12=M(1,2);
B13=M(1,3);
B14=M(1,4);
B15=1;
B16=dt*kx*wbz/2;
B17=dt*kx*wby/2;

B21=M(2,1);
B22=M(2,2);
B23=M(2,3);
B24=M(2,4);
B25=dt*ky*wbz/2;
B26=1;
B27=dt*ky*wbx/2;

B31=M(3,1);
B32=M(3,2);
B33=M(3,3);
B34=M(3,4);
B35=dt*kz*wby/2;
B36=dt*kz*wbx/2;
B37=1;

% F2 = F1(5:7,:)+M*dt;
% F = [F1(1:4,:);F2];

F = [A11 A12 A13 A14 A15 A16 A17;
    A21 A22 A23 A24 A25 A26 A27;
    A31 A32 A33 A34 A35 A36 A37;
    A41 A42 A43 A44 A45 A46 A47;
    B11 B12 B13 B14 B15 B16 B17;
    B21 B22 B23 B24 B25 B26 B27;
    B31 B32 B33 B34 B35 B36 B37];
 
end