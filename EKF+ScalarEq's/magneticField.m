function B = magneticField(R)
[m,n] = size(R);
mu_0 = 1.257e-6;
mu_e = 7.94e22;
%k_e = [0 0 -1];
%K_e = repmat(k_e',1,n); 
%B = mu_0*mu_e./(4*pi*dot(R,R).^(3/2)).*((3.*R.*dot(K_e,R)/dot(R,R).^(3/2)-K_e));
x = R(1); y = R(2); z = R(3);
r = sqrt(x^2+y^2+z^2);
km = mu_0*mu_e./(4*pi*r^5);
Bx = 3*km*x*z;
By = 3*km*y*z;
Bz = km*(2*z^2-x^2-y^2);
B = [Bx; By; Bz];
end