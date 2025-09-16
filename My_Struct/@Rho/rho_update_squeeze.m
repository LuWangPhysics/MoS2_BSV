          function obj=rho_update_squeeze(obj,MF,C)
           %  t_range=MF.t_start:MF.t_end;
             % V1=sqrt((MF.Eps_total(t_range)+MF.Eps_0(t_range))./(2*MF.Eps_total(t_range)));
             % V2=-MF.Omega_dipole(t_range)./sqrt(2*MF.Eps_total(t_range).*(MF.Eps_total(t_range)+MF.Eps_0(t_range)));
                  V1=sqrt((MF.Eps_total+MF.Eps_0)./(2*MF.Eps_total));
             V2=-MF.Omega_dipole./sqrt(2*MF.Eps_total.*(MF.Eps_total+MF.Eps_0));
             
             %---------------------------------------
             %define the f_1 functin for integration
             %---------------------------------------

             f1=-1i.*MF.Omega_dipole./C.hbar./2;
             %---------------------------------------
             %perform the integration
             %---------------------------------------

             %make sure the right dimension is multiplied
             int_eps_total=cumtrapz(MF.Eps_total);
             f1_HI=f1.*exp(-1i*(MF.dt*int_eps_total./C.hbar));
             first_int=obj.magic_operation_squeeze(C,f1_HI,MF);
            
             int_temp=2.*real(cumtrapz(conj(f1_HI).*first_int.*MF.dt));
             rho_k_s=V1.^2.*int_temp;
             obj.rho_2_s=obj.rho_2_s+rho_k_s./obj.Norm_k;
             Ct=exp(1i*(MF.dt*int_eps_total./C.hbar)).*first_int+V1.*V2;
             %---------------------------------------------
             %calculate the current 
             %---------------------------------------------
             v_arr={MF.v_x11,MF.v_x22,MF.v_x12;MF.v_y11,MF.v_y22,MF.v_y12};
             for d_iter=1:2
             j0_intra=-C.e.*(v_arr{d_iter,1}.*V1.^2+abs(V2).^2.*v_arr{d_iter,2});
             j0_inter=-C.e.*(2.*real(V1.*V2.*conj(v_arr{d_iter,3})));
             j1_intra=2*C.e.*(v_arr{d_iter,2}-v_arr{d_iter,1}).*real(V1.*conj(V2).*Ct);
             j1_inter=-2.*C.e.*real(conj(v_arr{d_iter,3}).*V1.^2.*Ct-conj(v_arr{d_iter,3}).*V2.^2.*conj(Ct));
             j2_intra=C.e.*(v_arr{d_iter,2}-v_arr{d_iter,1}).*(V1.^2-abs(V2).^2).*int_temp;
             j2_inter=-4.*C.e.*real(V2.*V1.*conj(v_arr{d_iter,3})).*int_temp;
   
             j_intra=(j0_intra+j1_intra+j2_intra)*obj.J_norm;
             j_inter=(j0_inter+j1_inter+j2_inter)*obj.J_norm;
             if d_iter==1
             obj.j_s_k=j_intra+j_inter;
             end
             obj.J_xy_no_t_s(4+d_iter,:)= obj.J_xy_no_t_s(4+d_iter,:)+  j_intra+j_inter;
             end 



          end