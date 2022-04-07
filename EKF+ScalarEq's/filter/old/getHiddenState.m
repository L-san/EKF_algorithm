function z = getHiddenState(attitude, Borb)
global q0
DCM = quat2DCM([q0; attitude(1:3)]);
B_b = DCM*Borb;
w = attitude(4:6);
z = [B_b; w];
end