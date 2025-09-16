function obj=renew_kt_ZnO(obj,C,kx_iter,ky_iter,kz_iter)
     %my E field is x-polarized
    Kt = obj.kx(kx_iter)+C.e.*obj.At./C.hbar;                                      %the k after shifted by A                          % kt = K-q A(t)/hbar 

    Ec=obj.eg;
    Ev=0;
    obj.coe_a={Kt.*obj.ax,obj.ky(ky_iter).*obj.ay,obj.kz(kz_iter).*obj.az};
    for i_d=1:length(obj.coe_c)
        for i_iter=1:length(obj.coe_c{i_d})
        Ec=Ec+obj.coe_c{i_d}(i_iter).*cos((i_iter-1).*obj.coe_a{i_d});
        Ev=Ev+obj.coe_v{i_d}(i_iter).*cos((i_iter-1).*obj.coe_a{i_d});
        end
    end
          
        obj.Eps_0=Ec-Ev;

         
        obj.d=obj.eg.*obj.d0./obj.Eps_0;
        obj.Omega_dipole=2*C.q.*obj.Et.*obj.d;                               
        obj.Eps_total=  sqrt(obj.Eps_0.^2 +obj.Omega_dipole.^2 );
        %obj.Eps_0_dot=obj.vx.*C.q.*obj.Et;
        obj.Omega_dipole_dot=diff([obj.Omega_dipole,obj.Omega_dipole(end)])./obj.dt;
        
        %current related parameters, not time dependent
        obj.v_x11=0;
        obj.v_x22=0;
     
        Ec=obj.eg;
        Ev=0;

            %derivateve for kx, v is not time dependent
           
            for i_iter=1:length(obj.coe_c{1})
            obj.v_x22=obj.v_x22-obj.coe_c{1}(i_iter).*(i_iter-1)*obj.ax.*sin((i_iter-1).*obj.kx(kx_iter).*obj.ax)./C.hbar;
            obj.v_x11=obj.v_x11-obj.coe_v{1}(i_iter).*(i_iter-1)*obj.ax.*sin((i_iter-1).*obj.kx(kx_iter).*obj.ax)./C.hbar;
             end

          obj.v_x12=1i*obj.eg.*obj.d0/C.hbar;
           obj.v_diff_x=obj.v_x22-obj.v_x11;

end