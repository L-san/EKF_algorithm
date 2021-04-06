function B = magneticField(R)
mu_0 = 1.257e-6;
mu_e = 7.94e22;
k_e = [0; 0; -1];
B = (mu_0*mu_e/(4*pi*norm(R)^3))*(3*R*dot(k_e,R)/norm(R)^2-k_e);
end