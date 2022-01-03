function z = getHiddenState(attitude, Borb)
B_Ax = Borb(1); B_Ay = Borb(2); B_Az = Borb(3);
q1 = attitude(1); q2 = attitude(2); q3 = attitude(3); q4 = attitude(4);
global orbit_vec
w0 = orbit_vec.w0;
Wx = attitude(5); Wy = attitude(6); Wz = attitude(7);

DCM = quat2dcm(attitude(1:4)');
B_b = DCM*Borb;
% 
% wox = (2*q2*q3 + 2*q1*q4)*w0;
% woy = (q1^2 - q2^2 + q3^2 - q4^2)*w0;
% woz = (-2*q1*q2 + 2*q3*q4)*w0;
% 
% w=[Wx+wox;Wy+woy;Wz+woz];
%w=[Wx;Wy;Wz];
w = attitude(5:7);
z = [0; B_b; w];
end