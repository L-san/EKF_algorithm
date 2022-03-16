function Jh = observationMatrix(attitude, Borb)
global orbit_vec q0
BAx = Borb(1); BAy = Borb(2); BAz = Borb(3);
q1 = attitude(1); q2 = attitude(2); q3 = attitude(3);
w0y = orbit_vec.w0;

dbxdq1=2*(BAx*q1+BAy*q2+BAz*q3);
dbxdq2=-2*(BAz*q0-BAy*q1+BAx*q2);
dbxdq3=2*(BAy*q0+BAz*q1-BAx*q3);

dbydq1=2*(BAz*q0-BAy*q1+BAx*q2);
dbydq2=2*(BAx*q1+BAy*q2+BAz*q3);
dbydq3=-2*(BAx*q0-BAz*q2+BAy*q3);

dbzdq1=-2*(BAy*q0+BAz*q1-BAx*q3);
dbzdq2=2*(BAx*q0-BAz*q2+BAy*q3);
dbzdq3=2*(BAx*q1+BAy*q2+BAz*q3);

dwxdq1=-2*q2*w0y;
dwxdq2=-2*q1*w0y;
dwxdq3=-2*q0*w0y;

dwydq1=2*q1*w0y;
dwydq2=-2*q2*w0y;
dwydq3=2*q3*w0y;

dwzdq1=2*q0*w0y;
dwzdq2=-2*q3*w0y;
dwzdq3=-2*q2*w0y;

dbdq = [dbxdq1 dbxdq2 dbxdq3;
	dbydq1 dbydq2 dbydq3;
	dbzdq1 dbzdq2 dbzdq3];

dbdw = zeros(3,3);

dwdq = [dwxdq1 dwxdq2 dwxdq3;
	dwydq1 dwydq2 dwydq3;
	dwzdq1 dwzdq2 dwzdq3];

dwdw = ones(3,3);

Jh = [dbdq dbdw;
      dwdq dwdw];
      
end

