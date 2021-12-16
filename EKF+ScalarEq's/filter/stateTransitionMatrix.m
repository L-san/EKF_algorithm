function F = stateTransitionMatrix(attitude,dt)
global satellite orbit_vec 
q0 = attitude(1); q1 = attitude(2); q2 = attitude(3); q3 = attitude(4);
wx = attitude(5); wy = attitude(6); wz = attitude(7);
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

% 
F1 = [1, -(1/4)*wx*dt, -(1/4)*wy*dt, -(1/4)*wz*dt, -(1/4)*q1*dt, -(1/4)*q2*dt,-(1/4)*q3*dt;
    (1/4)*wx*dt, 1, (1/4)*wz*dt, -(1/4)*wy*dt, (1/4)*q0*dt, -(1/4)*q3*dt, (1/4)*q2*dt;
    (1/4)*wy*dt, -(1/4)*wz*dt, 1, (1/4)*wx*dt, (1/4)*q3*dt, (1/4)*q0*dt, -(1/4)*q1*dt;
    (1/4)*wz*dt, (1/4)*wy*dt, -(1/4)*wx*dt, 1 -(1/4)*q2*dt, (1/4)*q1*dt, (1/4)*q0*dt;
     0, 0, 0, 0, 1, (1/2)*(Iy-Iz)*wz/Ix*dt, (1/2)*(Iy-Iz)*wy/Ix*dt;
     0, 0, 0, 0, (1/2)*(Iz-Ix)*wz/Iy*dt, 1, (1/2)*(Iz-Ix)*wx/Iy*dt;
     0, 0, 0, 0, (1/2)*(Ix-Iy)*wy/Iz*dt, (1/2)*(Ix-Iy)*wx/Iz*dt, 1];

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

M = Mg+Ma;
F2 = F1(5:7,:)+M;
F = [F1(1:4,:);F2];
 
end