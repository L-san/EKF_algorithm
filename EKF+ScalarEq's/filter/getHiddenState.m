function z = getHiddenState(attitude, Borb)
global q0 M A x0 alphaT
DCM = quat2DCM([q0; attitude(1:3)]);
B_b = DCM*Borb;
w = attitude(4:6);
m = M*(B_b)+x0;
g = A*(w)+alphaT;
z = [m; g];
end