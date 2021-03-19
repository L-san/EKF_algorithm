function [M] = calculateTorques(attitude)
global orbit_vec satellite

Aorb2b = quat2dcm(attitude(1:4)');
%%
%aerodynamics
ro = Density(orbit_vec.h);
statvec = kepel_statvec([orbit_vec.sma;
                         orbit_vec.ecc;
                         orbit_vec.inc;
                         orbit_vec.raan;
                         orbit_vec.aop;
                         orbit_vec.ta]);%convert to carthesian
V = norm(statvec(4:6));%circular velocity
e_V = Aorb2b*statvec(4:6)'/V;
Ma = -0.5*ro* satellite.c0* satellite.Sm* V^2*cross(satellite.ra,e_V);
%%
%gravity
e_B = Aorb2b*[0;0;1];
Mg = 3*orbit_vec.w0^2*cross(e_B,satellite.I*e_B);
%%all moments
M = Ma+Mg;

end