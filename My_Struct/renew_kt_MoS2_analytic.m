function MF=renew_kt_MoS2_analytic(MF,C,kx_iter,ky_iter,kz_iter)
        %the efield is rotated vs x axis by angle_m

        Ax=cos(MF.angle_m).*MF.At(MF.t_range);
        Ay=sin(MF.angle_m).*MF.At(MF.t_range);

        Kx = MF.kx(kx_iter)*MF.B(1,1)+MF.ky(ky_iter).*MF.B(1,2);
        Ky = MF.kx(kx_iter)*MF.B(2,1)+MF.ky(ky_iter).*MF.B(2,2);

        ddx=MF.dx_mk/1000;
        ddy=MF.dy_mk/1000;
        %the periodic boundary of the calculation arry not the data array
        for t_iter=1:length(MF.t_range)
                  %prepare the derivative of the H_k
                  Ktx = Kx+abs(C.e).*Ax(t_iter)./C.hbar;
                  Kty = Ky+abs(C.e).*Ay(t_iter)./C.hbar;

                   Hkt=MF.H_k(Ktx,Kty);
                   Hkt_xp=MF.H_k(Ktx+ddx,Kty);
                   Hkt_xn=MF.H_k(Ktx-ddx,Kty);
                   Hkt_yp=MF.H_k(Ktx,Kty+ddy);
                   Hkt_yn=MF.H_k(Ktx,Kty-ddy);
  
                      [eigen_vec,eigen_value] = eig(Hkt);
                      MF.Eps_0(t_iter)=eigen_value(2,2)-eigen_value(1,1);
                     eigen_vec=MF.gauge_fix(eigen_vec);       
                      dHx=(Hkt_xp-Hkt_xn)/(2*ddx*C.hbar);
                      dHy=(Hkt_yp-Hkt_yn)/(2*ddy*C.hbar);
                      MF.v_x11(t_iter)=real(MF.bra_ket_sum(eigen_vec(:,1),eigen_vec(:,1),dHx));
                      MF.v_x22(t_iter)=real(MF.bra_ket_sum(eigen_vec(:,2),eigen_vec(:,2),dHx));
                      MF.v_x12(t_iter)=MF.bra_ket_sum(eigen_vec(:,1),eigen_vec(:,2),dHx);
                      MF.v_y11(t_iter)=real(MF.bra_ket_sum(eigen_vec(:,1),eigen_vec(:,1),dHy));
                      MF.v_y22(t_iter)=real(MF.bra_ket_sum(eigen_vec(:,2),eigen_vec(:,2),dHy));
                      MF.v_y12(t_iter)=MF.bra_ket_sum(eigen_vec(:,1),eigen_vec(:,2),dHy);

                      MF.etax12(t_iter)=MF.v_x12(t_iter)./(-1i*  MF.Eps_0(t_iter)/C.hbar);
                      MF.etay12(t_iter)=MF.v_y12(t_iter)./(-1i*  MF.Eps_0(t_iter)/C.hbar);
                          


        end

    
      theta_r=(cos(MF.angle_m).*real(MF.etax12)+sin(MF.angle_m).*real(MF.etay12));
     theta_i=(cos(MF.angle_m).*imag(MF.etax12)+sin(MF.angle_m).*imag(MF.etay12));
     angle_fix=myatan_0Topi(theta_r,theta_i);

      MF.v_x12=MF.v_x12.*exp(-1i.*angle_fix);
       MF.v_y12=MF.v_y12.*exp(-1i.*angle_fix);
       MF.etax12=MF.etax12.*exp(-1i.*angle_fix);
       MF.etay12=MF.etay12.*exp(-1i.*angle_fix);
       MF.Omega_dipole=2*C.q.*MF.Et(MF.t_range).*real(cos(MF.angle_m).*MF.etax12+sin(MF.angle_m).*MF.etay12);


   
       MF.Eps_total=  sqrt(MF.Eps_0.^2 +abs(MF.Omega_dipole).^2 );
        %---------------------------------------------
       %the velocity part for BSV is not time dependent!
       %---------------------------------------------

       theta_sq=MF.angle_m;
        MF.v_sq=(MF.v_x22(1)-MF.v_x11(1))*cos(theta_sq)+sin(theta_sq)*(MF.v_y22(1)-MF.v_y11(1));


end