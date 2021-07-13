function [M] = calculateTorques(attitude)
global orbit_vec satellite

Aorb2b = quat2dcm(attitude(1:4)');
%%
%aerodynamics
ro = Density(orbit_vec.h);
V = attitude(11:13);%circular velocity
e_V = Aorb2b*V/norm(V);
Ma = -0.5*ro* satellite.c0* satellite.Sm* norm(V)^2*cross(satellite.ra,e_V);
%%
%gravity
e_B = Aorb2b*[0;0;1];
Mg = 3*orbit_vec.w0^2*cross(e_B,satellite.I*e_B);
%%all moments
%M = Ma+Mg;
M = Mg;
end