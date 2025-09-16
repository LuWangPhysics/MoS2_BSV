function plot_J_total_pv(obj,E,R,T, n_order)
            t=E.t(E.t_start:(E.t_start+length(R.rho_2)-1));
            D_x={R.J_xy_no_t_s(1,:),R.J_xy_no_t_s(3,:),...
                      R.J_xy_no_t_s(5,:)};
            D_y={R.J_xy_no_t_s(2,:),R.J_xy_no_t_s(4,:),...
                      R.J_xy_no_t_s(6,:)};

            name_arr={'no heat bath', 'thermal ', 'squeeze'}; 
             %save as non pv, thermal pv, sq pv
            I_value=zeros(2*length(D_x),length(n_order));
           figure(Name='Jp,Jv')
            for i_iter=1:length(D_x)
   %               
                      jp=cos(E.angle_m).*D_x{i_iter}+sin(E.angle_m).*D_y{i_iter};
                      jv=-sin(E.angle_m).*D_x{i_iter}+cos(E.angle_m).*D_y{i_iter};

                      dt=t(2)-t(1);
                      f=linspace(-1/(2*dt),1/(2*dt),length(t));
                      df=1/dt/length(t);

  
                      subplot(3,2,2*i_iter-1)
                      jp_f=fftshift(fft(ifftshift(jp)))./(length(t)*df);
                      plot(f./(E.om/(2*pi)),abs(jp_f))
                      set(gca, 'YScale', 'log')
                      title([name_arr{i_iter} 'p'])
                      xlabel('Harmonics orders')
                      xlim([0,20])
           
                      subplot(3,2,2*i_iter)
                      jv_f=fftshift(fft(ifftshift(jv)))./(length(t)*df);
                      plot(f./(E.om/(2*pi)),abs(jv_f))
                      set(gca, 'YScale', 'log')
                      title([name_arr{i_iter} 'v'])
                      xlabel('Harmonics orders')
                      xlim([0,20])

                     f_location=fix(length(f)/2)+fix(n_order.*E.om/(2*pi*df));
                 
                     for f_iter=1:length(n_order)
                          I_value(2*i_iter-1,f_iter)=jp_f(f_location(f_iter));
                          I_value(2*i_iter,f_iter)=jv_f(f_location(f_iter));
                     end

                     
	    end


                         if obj.save_flag==1
                           h = gcf; % gcf returns the handle to the current figure
                           savefig(h, [obj.save_str   'J_spec_pv.fig']);
                            save([obj.save_str  '_E_harmonic_orders_value_pv.mat'],"I_value"); 
                           save([obj.save_str  '_E_harmonic_orders.mat'],"n_order"); 
                       end
        end
