function z = getHiddenState(attitude, Borb)

DCM = quat2DCM(attitude(1:4)');
B_b = DCM*Borb;

w = attitude(5:7);
z = [0; B_b; w];
end