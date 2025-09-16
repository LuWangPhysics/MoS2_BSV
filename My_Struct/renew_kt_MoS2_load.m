function MF=renew_kt_MoS2_load(MF,C,kx_iter,ky_iter,kz_iter)
        %the efield is rotated vs x axis by angle_m
        t_range=MF.t_start:MF.t_end;
        Ax=cos(MF.angle_m).*MF.At(t_range);
        Ay=sin(MF.angle_m).*MF.At(t_range);
        %the momentum shift in kx ky space shift to dx, dy monkhorst space

        Ktx = MF.kx(kx_iter)+C.e.*(MF.B_inv(1,1).*Ax+MF.B_inv(1,2).*Ay)./C.hbar;
        Kty = MF.ky(ky_iter)+C.e.*(MF.B_inv(2,1).*Ax+MF.B_inv(2,2).*Ay)./C.hbar;
        %the periodic boundary of the calculation arry not the data array
        T_x=abs(MF.kx_d(1));
        T_y=abs(MF.ky_d(1));

        flag_sign=fix(Ktx./T_x);
        kx_eff=rem(Ktx,T_x)+(-1).^flag_sign.*rem(flag_sign,2).*T_x;

        flag_sign=fix(Kty./T_y);
        ky_eff=rem(Kty,T_y)+(-1).^flag_sign.*rem(flag_sign,2).*T_y;

       %construct shift vector
       % R_x12=transpose(interpn(obj.kx_d,obj.ky_d,obj.R_x12_d,kx_eff,ky_eff,'spline'));


       %time dependent d_eff use absolute value
       MF.etax12=interp2(MF.kx_d,MF.ky_d,MF.etax_12_d,kx_eff,ky_eff,'spline');
       MF.etay12=interp2(MF.kx_d,MF.ky_d,MF.etay_12_d,kx_eff,ky_eff,'spline');
       MF.v_x11=interp2(MF.kx_d,MF.ky_d,MF.v_x11_d,kx_eff,ky_eff,'spline');
        MF.v_x12=interp2(MF.kx_d,MF.ky_d,MF.v_x12_d,kx_eff,ky_eff,'spline');
       MF.v_y12=interp2(MF.kx_d,MF.ky_d,MF.v_y12_d,kx_eff,ky_eff,'spline');
       MF.v_x22=interp2(MF.kx_d,MF.ky_d,MF.v_x22_d,kx_eff,ky_eff,'spline');
       MF.v_y11=interp2(MF.kx_d,MF.ky_d,MF.v_y11_d,kx_eff,ky_eff,'spline');
       MF.v_y22=interp2(MF.kx_d,MF.ky_d,MF.v_y22_d,kx_eff,ky_eff,'spline');
       MF.Eps_0=interp2(MF.kx_d,MF.ky_d,MF.Eps_d,kx_eff,ky_eff,'spline');

      theta_r=(cos(MF.angle_m).*real(MF.etax12)+sin(MF.angle_m).*real(MF.etay12));
     theta_i=(cos(MF.angle_m).*imag(MF.etax12)+sin(MF.angle_m).*imag(MF.etay12));
     angle_fix=myatan_0Topi(theta_r,theta_i);

      MF.v_x12=MF.v_x12.*exp(-1i.*angle_fix);
       MF.v_y12=MF.v_y12.*exp(-1i.*angle_fix);
       MF.etax12=MF.etax12.*exp(-1i.*angle_fix);
       MF.etay12=MF.etay12.*exp(-1i.*angle_fix);
       MF.Omega_dipole=2*C.q.*MF.Et(t_range).*real(cos(MF.angle_m).*MF.etax12+sin(MF.angle_m).*MF.etay12);

 

  


       MF.Eps_total=  sqrt(MF.Eps_0.^2 +abs(MF.Omega_dipole).^2 );


        theta_sq=0;
        MF.v_sq=(MF.v_x22(1)-MF.v_x11(1))*cos(theta_sq)+sin(theta_sq)*(MF.v_y22(1)-MF.v_y11(1));

end