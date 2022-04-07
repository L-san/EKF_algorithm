function [x_hat_f, P] = kalmanReal(R, Q, z, P, F, H, h, f)
m = 6; n = 6;
%prediction
x_hat = f;
P = F*P*F'+Q;
%correction
S = H*P*H'+R;
Sx = inv(S);
%Sx = (S'*S)\S';
K = P*H'*Sx;
x_hat_f = x_hat+ K*(z-h);
P = (P-P*K*H);
end