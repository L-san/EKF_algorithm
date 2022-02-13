%%receiving data-array sized m-by-n and filtering the simplest way
%%returns filtered array
%%initArr should consist data, where m - number of sensor(channel) in array
%%and n is for number of each element of sensor's array
%%x - previous filtered
%%P - previous cov. matrix
%%z - now unfiltered
%%R_coeff and Qcoeff are coeff-s for filtering. Guess them.
%%note: R_coeff is for measurement uncertainty - increase for more
%%smoothing and vice versa(plot -> direct line)
%------------------------------------
%%plans: add variance-dependence here
function [x_hat_f, P] = kalmanReal(Rcoeff, Qcoeff, z, P, JF, JH, h, f)
m = 7; n = 7;
Im = eye(m);
In = eye(n);
R = In*Rcoeff;%measurement covariance
Q = Im*Qcoeff;%process covariance
%prediction
x_hat = f;
P = JF*P*JF'+Q;
%correction
S = JH*P*JH'+R;
%answers can be close to singular matrixes 
if det(S) ==0  
    K = zeros(m,n);
elseif isnan(det(S))
    K = zeros(m,n);
elseif det(S)==Inf
    K = zeros(m,n);
elseif det(S)==-Inf
    K = zeros(m,n);
else 
    K = P*JH'*inv(S)';
end

x_hat_f = x_hat+ K*(z-h);
%P = (eye(m)-K*JH)*P;
P = (P-P*K*JH);
end