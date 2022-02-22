function Jh = observationMatrix(attitude, Borb)
global orbit_vec
BAx = Borb(1); BAy = Borb(2); BAz = Borb(3);
q0 = getQ0(attitude(1:3)); q1 = attitude(1); q2 = attitude(2); q3 = attitude(3);
w0y = orbit_vec.w0;

dbxdq0=2*(BAx*q0-BAz*q2+BAy*q3);
dbxdq1=2*(BAx*q1+BAy*q2+BAz*q3);
dbxdq2=-2*(BAz*q0-BAy*q1+BAx*q2);
dbxdq3=2*(BAy*q0+BAz*q1-BAx*q3);

dbydq0=2*(BAy*q0+BAz*q1-BAx*q3);
dbydq1=2*(BAz*q0-BAy*q1+BAx*q2);
dbydq2=2*(BAx*q1+BAy*q2+BAz*q3);
dbydq3=-2*(BAx*q0-BAz*q2+BAy*q3);

dbzdq0=2*(BAz*q0-BAy*q1+BAx*q3);
dbzdq1=-2*(BAy*q0+BAz*q1-BAx*q2);
dbzdq2=2*(BAx*q1-BAz*q2+BAy*q3);
dbzdq3=2*(BAx*q0+BAy*q2+BAz*q3);

dwxdq0=-2*q3*w0y;
dwxdq1=-2*q2*w0y;
dwxdq2=-2*q2*w0y;
dwxdq3=-2*q0*w0y;
dwxdwx=1;
dwxdwy=0;
dwxdwz=0;

dwydq0=-2*q0*w0y;
dwydq1=2*q1*w0y;
dwydq2=-2*q2*w0y;
dwydq3=2*q3*w0y;
dwydwx=0;
dwydwy=1;
dwydwz=0;

dwzdq0=2*q1*w0y;
dwzdq1=2*q0*w0y;
dwzdq2=-2*q3*w0y;
dwzdq3=-2*q2*w0y;
dwzdwx=0;
dwzdwy=0;
dwzdwz=1;

% dbdq = [dbxdq0 dbxdq1 dbxdq2 dbxdq3;
% 	dbydq0 dbydq1 dbydq2 dbydq3;
% 	dbzdq0 dbzdq1 dbzdq2 dbzdq3];
% 
% dbdw = [0 0 0;
% 	0 0 0;
% 	0 0 0];
% 
% dwdq = [dwxdq0 dwxdq1 dwxdq2 dwxdq3;
% 	dwydq0 dwydq1 dwydq2 dwydq3
% 	dwzdq0 dwzdq1 dwzdq2 dwzdq3];
% 
% dwdw = [dwxdwx dwxdwy dwxdwz;
%         dwydwx dwydwy dwydwz;
%         dwzdwx dwzdwy dwzdwz];

dbdq = [dbxdq1 dbxdq2 dbxdq3;
        dbydq1 dbydq2 dbydq3;
        dbzdq1 dbzdq2 dbzdq3];

dbdw = [0 0 0;
	0 0 0;
	0 0 0];

dwdq = [dwxdq1 dwxdq2 dwxdq3;
        dwydq1 dwydq2 dwydq3
        dwzdq1 dwzdq2 dwzdq3];

dwdw = [dwxdwx dwxdwy dwxdwz;
        dwydwx dwydwy dwydwz;
        dwzdwx dwzdwy dwzdwz];
    
Jh = [dbdq dbdw;
      dwdq dwdw];
end

