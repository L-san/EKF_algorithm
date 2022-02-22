function z = getHiddenState(attitude, Borb)
q0 = getQ0(attitude(1:3));
DCM = quat2DCM([q0; attitude(1:3)]);
B_b = DCM*Borb;
w = attitude(4:6);
z = [B_b; w];
end