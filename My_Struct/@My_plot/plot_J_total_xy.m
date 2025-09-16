        function plot_J_total_xy(obj,E,R,T, n_order)
            t=E.t(E.t_start:(E.t_start+length(R.rho_2)-1));
            % if (sum(sum(abs(R.j_no_bath_intra)))~=0)
            D_arr={R.J_xy_no_t_s(1,:),R.J_xy_no_t_s(2,:),R.J_xy_no_t_s(3,:),R.J_xy_no_t_s(4,:)...
                      R.J_xy_no_t_s(5,:),R.J_xy_no_t_s(6,:)};
            name_arr={'no heat bathx','no heat bathy',' thermal x',' thermal y', 'squeezex','squeeze y'}; 
            % else
            % D_arr={    R.J_xy_no_t_s(5,:),R.J_xy_no_t_s(6,:)};
            % name_arr={'squeeze'}; 
            % end
            I_value=zeros(length(D_arr),length(n_order));
             figure(Name='Jx,Jy')
            for i_iter=1:length(D_arr)
   %
                      subplot(3,2,i_iter)
                      dt=t(2)-t(1);
                      f=linspace(-1/(2*dt),1/(2*dt),length(t));
                      df=1/dt/length(t);
                      jf=fftshift(fft(ifftshift(D_arr{i_iter})))./(length(t)*df);
                      plot(f./(E.om/(2*pi)),abs(jf))
                      set(gca, 'YScale', 'log')
                      title(name_arr{i_iter})
                      xlabel('Harmonics orders')
                      xlim([0,20])
           
                     f_location=fix(length(f)/2)+fix(n_order.*E.om/(2*pi*df));
                 
                     for f_iter=1:length(n_order)
                          I_value(i_iter,f_iter)=jf(f_location(f_iter));
                     end
                     
	    end

                        if obj.save_flag==1
                           h = gcf; % gcf returns the handle to the current figure
                           savefig(h, [obj.save_str 'J_spec_xy.fig']);
                            save([obj.save_str  '_E_harmonic_orders_value_xy.mat'],"I_value"); 
                            save([obj.save_str  '_E_harmonic_orders.mat'],"n_order"); 
                        end

        end
