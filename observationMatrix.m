function H = observationMatrix(attitude, A_I, f)
q0 = attitude(1); q1 = attitude(2); q2 = attitude(3); q3 = attitude(4);
wx = attitude(6); wy = attitude(7); wz = attitude(8);
ax = A_I(2); ay = A_I(3); az = A_I(4);

M21 = 2*ax*q0+2*az*q2-2*ay*q3;
M22 = 2*ax*q1+2*ay*q2+2*az*q3;
M23 = 2*az*q0+2*ay*q1-2*ax*q2;
M24 = -2*ay*q0+2*az*q1-2*ax*q3;

M31 = 2*ay*q0-2*az*q1+2*ax*q3;
M32 = -2*az*q0-2*ay*q1+2*ax*q2;
M33 = 2*ax*q1+2*ay*q2+2*az*q3;
M34 = 2*ax*q0+2*az*q2-2*ay*q3;

M41 = 2*az*q0+2*ay*q1-2*ax*q2;
M42 = 2*ay*q0-2*az*q1+2*ax*q3;
M43 = -2*ax*q0-2*az*q2+2*ay*q3;
M44 = 2*ax*q1+2*ay*q2+2*az*q3;


M = [0 0 0 0;
     M21 M22 M23 M24;
     M31 M32 M33 M34;
     M41 M42 M43 M44];

N = zeros(4,4);

K21 = az*(2*q3*wx+2*q0*wy)+ay*(2*q2*wx-2*q0*wz)-2*ax*(q2*wy+q3*wz);
K22 = az*(-2*q2*wx+2*q1*wy)+ay*(2*q3*wx-2*q1*wz)-ax*(q3*wy-q2*wz);
K23 = az*(-2*q1*wx-2*q2*wy)-2*ax*(q0*wy-q1*wz)+ay*(2*q0*wx+2*q2*wz);
K24 = az*(2*q0*wx-2*q3*wy)-2*ax*(q1*wy+q0*wz)+ay*(2*q1*wx+2*q3*wz);
L21 = 0;
L22 = az*(-2*q1*q2+2*q0*q3)+ay*(2*q0*q2+2*q1*q3);
L23 = -2*ax*(q0*q2+q1*q3)+az*(q0^2+q1^2-q2^2-q3^2);
L24 = -2*ax*(-q1*q2+q0*q3)+ay*(-q0^2-q1^2+q2^2+q3^2);


K31 = az*(-2*q0*wx+2*q3*wy)+ax*(2*q1*wy+2*q0*wz)-2*ay*(q1*wx+q3*wz);
K32 = az*(2*q1*wx+2*q2*wy)+ax*(2*q0*wy-2*q1*wz)-2*ay*(q0*wx+q2*wz);
K33 = az*(-2*q2*wx+2*q1*wy)-2*ay*(-q3*wx+q1*wz)+ax*(-2*q3*wy+2*q2*wz);
K34 = az*(2*q3*wx+2*q0*wy)-2*ay*(-q2*wx+q0*wz)+ax*(-2*q2*wy-2*q3*wz);
L31 = 0;
L32 = -2*ay*(q0*q1-q2*q3)+az*(-q0^2+q1^2-q2^2+q3^2);
L33 = az*(2*q1*q2+2*q0*q3)+ax*(2*q0*q1-2*q2*q3);
L34 = -2*ay*(q1*q2+q0*q3)+ax*(q0^2-q1^2+q2^2-q3^2);

K41 = -2*az*(q1*wx+q2*wy)+ax*(-2*q0*wy+2*q1*wz)+ay*(2*q0*wx+2*q2*wz);
K42 = -2*az*(q0*wx-q3*wy)+ax*(2*q1*wy+2*q0*wz)+ay*(-2*q1*wx-2*q3*wz);
K43 = -2*az*(q3*wx+q0*wy)+ay*(-2*q2*wx+2*q0*wz)+ax*(2*q2*wy+2*q3*wz);
K44 = -2*az*(q2*wx-q1*wy)+ay*(2*q3*wx-2*q1*wz)+ax*(-2*q3*wy+2*q2*wz);
L41 = 0;
L42 = -2*az*(q0*q1+q2*q3)+ay*(q0^2-q1^2-q2^2+q3^2);
L43 = -2*az*(q0*q2-q1*q3)+ax*(-q0^2+q1^2+q2^2-q3^2);
L44 = ay*(2*q0*q2-2*q1*q3)+ax*(2*q0*q1+2*q2*q3);

K = [0 0 0 0;
     K21 K22 K23 K24;
     K31 K32 K33 K34;
     K41 K42 K43 K44];
 
L = [0 0 0 0;
    L21 L22 L23 L24;
    L31 L32 L33 L34;
    L41 L42 L43 L44];

F = diag(f);
 if det(F)==0
     invF = zeros(8);
 else
     invF = diag(F);
 end

H1 = 0.5*[M N];
H2 = [zeros(4,4) eye(4,4)];
H2 = [K L]*(eye(8)-diag(attitude)*invF)*0.3335;
H = [H1;H2];
end

