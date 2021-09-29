function z = getHiddenState(attitude, Borb)
B_Ax = Borb(1); B_Ay = Borb(2); B_Az = Borb(3);
q1 = attitude(1); q2 = attitude(2); q3 = attitude(3); q4 = attitude(4);
global orbit_vec
w0 = orbit_vec.w0;
Wx = attitude(5); Wy = attitude(6); Wz = attitude(7);
Bb1 = 0;
Bb2 = 2*(B_Az*q1*q3+B_Ay*q2*q3-B_Ay*q1*q4+B_Az*q2*q4)+B_Ax*(q1^2+q2^2-q3^2-q4^2);
Bb3 = 2*(-B_Az*q1*q2+B_Ax*q2*q3+B_Ax*q1*q4+B_Az*q3*q4)+B_Ay*(q1^2-q2^2+q3^2-q4^2);
Bb4 = 2*(B_Ay*q1*q2-B_Ax*q1*q3+B_Ax*q2*q4+B_Ay*q3*q4)+B_Az*(q1^2-q2^2-q3^2+q4^2);
B_b = [Bb1; Bb2; Bb3; Bb4];
%B_b = quat_mult(quat_mult(attitude(1:4),[0; Borb]),quat_conj(attitude(1:4)));
% 
% wox = (2*q2*q3 - 2*q1*q4)*w0;
% woy = (q1^2 - q2^2 + q3^2 - q4^2)*w0;
% woz = (-2*q1*q2 + 2*q3*q4)*w0;
% 
% w=[Wx+2*wox;Wy+2*woy;Wz+2*woz];
w=[Wx;Wy;Wz];
z = [B_b; w];
end