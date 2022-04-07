function Jh = observationMatrix(attitude, Borb)
global orbit_vec q0 M A
M11 = M(1,1); M12 = M(1,2); M13 = M(1,3);
M21 = M(2,1); M22 = M(2,2); M23 = M(2,3);
M31 = M(3,1); M32 = M(3,2); M33 = M(3,3);
A11 = A(1,1); A12 = A(1,2); A13 = A(1,3);
A21 = A(2,1); A22 = A(2,2); A23 = A(2,3);
A31 = A(3,1); A32 = A(3,2); A33 = A(3,3);

BAx = Borb(1); BAy = Borb(2); BAz = Borb(3);
q1 = attitude(1); q2 = attitude(2); q3 = attitude(3);
w0y = orbit_vec.w0;

dbxdq1=M12*(2*BAz*q0-2*BAy*q1+2*BAx*q2)+M13*(-2*BAy*q0-2*BAz*q1+2*BAx*q3)+M11*(2*BAx*q1+2*BAy*q2+2*BAz*q3);
dbxdq2=M11*(-2*BAz*q0+2*BAy*q1-2*BAx*q2)+M13*(2*BAx*q0-2*BAz*q2+2*BAy*q3)+M12*(2*BAx*q1+2*BAy*q2+2*BAz*q3);
dbxdq3=M11*(2*BAy*q0+2*BAz*q1-2*BAx*q3)+M12*(-2*BAx*q0+2*BAz*q2-2*BAy*q3)+M13*(2*BAx*q1+2*BAy*q2+2*BAz*q3);

dbydq1=M22*(2*BAz*q0-2*BAy*q1+2*BAx*q2)+M23*(-2*BAy*q0-2*BAz*q1+2*BAx*q3)+M21*(2*BAx*q1+2*BAy*q2+2*BAz*q3);
dbydq2=M21*(-2*BAz*q0+2*BAy*q1-2*BAx*q2)+M23*(2*BAx*q0-2*BAz*q2+2*BAy*q3)+M22*(2*BAx*q1+2*BAy*q2+2*BAz*q3);
dbydq3=M21*(2*BAy*q0+2*BAz*q1-2*BAx*q3)+M22*(-2*BAx*q0+2*BAz*q2-2*BAy*q3)+M23*(2*BAx*q1+2*BAy*q2+2*BAz*q3);

dbzdq1=M32*(2*BAz*q0-2*BAy*q1+2*BAx*q2)+M33*(-2*BAy*q0-2*BAz*q1+2*BAx*q3)+M31*(2*BAx*q1+2*BAy*q2+2*BAz*q3);
dbzdq2=M31*(-2*BAz*q0+2*BAy*q1-2*BAx*q2)+M33*(2*BAx*q0-2*BAz*q2+2*BAy*q3)+M32*(2*BAx*q1+2*BAy*q2+2*BAz*q3);
dbzdq3=M31*(2*BAy*q0+2*BAz*q1-2*BAx*q3)+M32*(-2*BAx*q0+2*BAz*q2-2*BAy*q3)+M33*(2*BAx*q1+2*BAy*q2+2*BAz*q3);

dwxdwx=A11;
dwxdwy=A12;
dwxdwz=A13;

dwydwx=A21;
dwydwy=A22;
dwydwz=A23;

dwzdwx=A31;
dwzdwy=A32;
dwzdwz=A33;

dbdq = [dbxdq1 dbxdq2 dbxdq3;
        dbydq1 dbydq2 dbydq3;
        dbzdq1 dbzdq2 dbzdq3];

dbdw = zeros(3,3);

dwdq = zeros(3,3);

dwdw = [dwxdwx dwxdwy dwxdwz;
        dwydwx dwydwy dwydwz;
        dwzdwx dwzdwy dwzdwz];

Jh = [dbdq dbdw;
      dwdq dwdw]; 
end

