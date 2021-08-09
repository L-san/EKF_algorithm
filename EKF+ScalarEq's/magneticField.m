function B = magneticField(R)
[m,n] = size(R);
mu_0 = 1.257e-6;
mu_e = 7.94e22;
k_e = [0 0 -1];
K_e = repmat(k_e',1,n); 
B = mu_0*mu_e./(4*pi*dot(R,R).^(3/2)).*((3.*R.*dot(K_e,R)/dot(R,R).^(3/2)-K_e));
end