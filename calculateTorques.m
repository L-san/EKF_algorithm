function [M] = calculateTorques(attitude)
global orbit_vec sattelite mu
q0 = attitude(1);
q1 = attitude(2);
q2 = attitude(3);
q3 = attitude(4);

R = orbit_vec(1);
h = orbit_vec(8);

c0 = sattelite(1);
Sm = sattelite(2);
dx = sattelite(3);
Ix = sattelite(4);
Iy = sattelite(5);
Iz = sattelite(6);

ro = Density(h);
statvec = kepel_statvec(orbit_vec(1:6));%convert to carthesian
V = norm(statvec(4:6));%circular velocity
%%
%aerodynamics
Max = 0;
May = 0.5*c0*Sm*(abs(q0^2 + q1^2 - q2^2 - q3^2) +...
                3*(2*abs(q0*q2 +q1*q3) + 2*abs(q1*q2 - q0*q3)))*...
                ro*V^2*dx*(2*(q0*q2 + q1*q3));
Maz = -0.5*c0*Sm*(abs(q0^2 + q1^2 - q2^2 - q3^2) +...
                  3*(2*abs(q0*q2 +q1*q3) + 2*abs(q1*q2 - q0*q3)))*...
                  ro*V^2*dx*(2*(q1*q2 - q0*q3)); 
Ma = [Max; May; Maz];
%%
%gravity
Mgx = 3*mu/R^3*(Iz-Iy)*(2*(q0*q1+q2*q3))*(q0^2 + q3^2 - q2^2 - q1^2);
Mgy = 3*mu/R^3*(Ix-Iz)*(2*(q1*q3-q0*q2))*(q0^2 + q3^2 - q2^2 - q1^2);
Mgz = 3*mu/R^3*(Iy-Ix)*(2*(q1*q3-q0*q2))*(2*(q0*q1+q2*q3));
Mg = [Mgx; Mgy; Mgz];
%%all moments
M = Ma+Mg;
%M = [0 0 0];
end