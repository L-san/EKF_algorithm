function dcm = quat2DCM( q )
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    dcm11 = q0^2 + q1^2 - q2^2 - q3^2;
    dcm12 = 2*(q1*q2 + q0*q3);
    dcm13 = 2*(q1*q3 - q0*q2);
    dcm21 = 2*(q1*q2 - q0*q3);
    dcm22 = q0^2 - q1^2 + q2^2 - q3^2;
    dcm23 = 2*(q2*q3 + q0*q1);
    dcm31 = 2*(q1*q3 + q0*q2);
    dcm32 = 2*(q2*q3 - q0*q1);
    dcm33 = q0^2 - q1^2 - q2^2 + q3^2;
    dcm = [dcm11 dcm12 dcm13;
           dcm21 dcm22 dcm23;
           dcm31 dcm32 dcm33];
end