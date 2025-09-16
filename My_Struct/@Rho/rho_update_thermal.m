          function obj=rho_update_thermal(obj,MF,C,jo)
             % t_range=MF.t_start:MF.t_end;
             % V1=sqrt((MF.Eps_total(t_range)+MF.Eps_0(t_range))./(2*MF.Eps_total(t_range)));
             % V2=-MF.Omega_dipole(t_range)./sqrt(2*MF.Eps_total(t_range).*(MF.Eps_total(t_range)+MF.Eps_0(t_range)));
               V1=sqrt((MF.Eps_total+MF.Eps_0)./(2*MF.Eps_total));
             V2=-MF.Omega_dipole./sqrt(2*MF.Eps_total.*(MF.Eps_total+MF.Eps_0));
             %---------------------------------------
             %define the f_1 functin for integration
             %---------------------------------------

             f1=-1i.*MF.Omega_dipole./C.hbar./2-2*V2.*[diff(V1),V1(end)-V1(end-1)]./MF.dt;
             %---------------------------------------
             %perform the integration
             %---------------------------------------

             %make sure the right dimension is multiplied
             int_eps_total=cumtrapz(MF.Eps_total);
             f1_HI=f1.*exp(-1i*(MF.dt*int_eps_total./C.hbar));
             first_int=obj.magic_operation(f1_HI,MF);


             int_temp=2.*real(cumtrapz(conj(f1_HI).*first_int.*MF.dt));
             %for a given k
             obj.rho_k=V1.^2.*int_temp;
             obj.rho_2=obj.rho_2+obj.rho_k./obj.Norm_k;
             % 
             % %---------------------------------------------
             % %calculate the current 
             % %---------------------------------------------
             % %since E_band is always symmetry, j0 +k -k always cancels out
             % %and does not contirubte to overall HHG
             Ct=exp(1i*(MF.dt*int_eps_total./C.hbar)).*first_int+V1.*V2;
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
             obj.j_t_k=j_intra+j_inter;
              end
             obj.J_xy_no_t_s(2+d_iter,:)=   obj.J_xy_no_t_s(2+d_iter,:)+j_intra+j_inter;
             end
             obj.j_thermal_intra=obj.j_thermal_intra+j_intra;
             obj.j_thermal_inter=obj.j_thermal_inter+j_inter;
       
             %---------------------------------------------
             %no heat bath
             %---------------------------------------------          

              int_temp=abs(cumtrapz(f1_HI).*MF.dt).^2;
             obj.rho_k_no=V1.^2.*int_temp;
             obj.rho_2_no_heat=obj.rho_2_no_heat+obj.rho_k_no./obj.Norm_k;
         
             Ct_no=exp(1i*(MF.dt*int_eps_total./C.hbar)).*cumtrapz(f1_HI).*MF.dt+V1.*V2;
           
             for d_iter=1:2
              j0_intra=-C.e.*(v_arr{d_iter,1}.*V1.^2+abs(V2).^2.*v_arr{d_iter,2});
              j0_inter=-C.e.*(2.*real(V1.*V2.*conj(v_arr{d_iter,3})));
               j1_no_intra=0*2*C.e.*(v_arr{d_iter,2}-v_arr{d_iter,1}).*real(V1.*conj(V2).*Ct_no);
               j1_no_inter=-2.*C.e.*real(conj(v_arr{d_iter,3}).*V1.^2.*Ct_no-conj(v_arr{d_iter,3}).*V2.^2.*conj(Ct_no));
             j2_no_intra=C.e.*(v_arr{d_iter,2}-v_arr{d_iter,1}).*(V1.^2-abs(V2).^2).*int_temp;
              j2_no_inter=-4.*C.e.*real(V2.*V1.*conj(v_arr{d_iter,3})).*int_temp;
             
              % 
              j_no_intra=(j0_intra+j1_no_intra+j2_no_intra)*obj.J_norm;
              j_no_inter=(j0_inter+j1_no_inter+j2_no_inter)*obj.J_norm;


	   obj.j_no_bath_intra=obj.j_no_bath_intra+j_no_intra;
             obj.j_no_bath_inter=obj.j_no_bath_inter+j_no_inter;   
             if d_iter==1
                   obj.j_no_k=j_no_intra+j_no_inter;
             end
              obj.J_xy_no_t_s(d_iter,:)= obj.J_xy_no_t_s(d_iter,:)+  j_no_intra+j_no_inter;
             end
             %------------------------------------------------
             %analytical resuslts for two-level not two band!!!
             %------------------------------------------------
              % T2=C.hbar/(2*jo*pi*C.Kb*C.T);
              % f_two_level=MF.Omega_dipole./C.hbar/2;
              % S=cumtrapz(MF.Eps_0).*MF.dt./C.hbar;

             %obj.rho_nm=-1i.*exp(1i*(MF.dt*int_eps_total./C.hbar)).*cumtrapz(f1_HI).*MF.dt./obj.Norm_k;
             % obj.rho_nm=-1i.*exp(-MF.t./T2+1i.*S).*cumtrapz(f_two_level.*exp(MF.t./T2-1i.*S)).*MF.dt./obj.Norm_k;
             % obj.rho_nn=-2.*cumtrapz(f_two_level.*imag(obj.rho_nm).*MF.dt);
             % obj.rho_2=obj.rho_2+obj.rho_nn(t_range);
          end