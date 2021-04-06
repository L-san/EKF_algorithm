function F = stateTransitionMatrix(attitude, I,dt)
Ix = I(1,1); Iy = I(2,2); Iz = I(3,3);
q0 = attitude(1); q1 = attitude(2); q2 = attitude(3); q3 = attitude(4);
wx = attitude(6); wy = attitude(7); wz = attitude(8);
% 
F = [1, -(1/4)*wx*dt, -(1/4)*wy*dt, -(1/4)*wz*dt, 0, -(1/4)*q1*dt, -(1/4)*q2*dt,-(1/4)*q3*dt;
    (1/4)*wx*dt, 1, (1/4)*wz*dt, -(1/4)*wy*dt, 0, (1/4)*q0*dt, -(1/4)*q3*dt, (1/4)*q2*dt;
    (1/4)*wy*dt, -(1/4)*wz*dt, 1, (1/4)*wx*dt, 0, (1/4)*q3*dt, (1/4)*q0*dt, -(1/4)*q1*dt;
    (1/4)*wz*dt, (1/4)*wy*dt, -(1/4)*wx*dt, 1, 0, -(1/4)*q2*dt, (1/4)*q1*dt, (1/4)*q0*dt;
     0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 1, (1/2)*(Iy-Iz)*wz/Ix*dt, (1/2)*(Iy-Iz)*wy/Ix*dt;
     0, 0, 0, 0, 0, (1/2)*(Iz-Ix)*wz/Iy*dt, 1, (1/2)*(Iz-Ix)*wx/Iy*dt;
     0, 0, 0, 0, 0, (1/2)*(Ix-Iy)*wy/Iz*dt, (1/2)*(Ix-Iy)*wx/Iz*dt, 1];
 
end