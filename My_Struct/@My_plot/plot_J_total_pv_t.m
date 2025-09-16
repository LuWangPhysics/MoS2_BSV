function plot_J_total_pv_t(obj,E,R,T, n_order)
            t=E.t(E.t_start:(E.t_start+length(R.rho_2)-1));
            D_x={R.J_xy_no_t_s(1,:),R.J_xy_no_t_s(3,:),...
                      R.J_xy_no_t_s(5,:)};
            D_y={R.J_xy_no_t_s(2,:),R.J_xy_no_t_s(4,:),...
                      R.J_xy_no_t_s(6,:)};

            name_arr={'no heat bath', 'thermal ', 'squeeze'}; 
             %save as non pv, thermal pv, sq pv
           figure(Name='Jp,Jv')
            for i_iter=1:length(D_x)
   %               
                      jp=cos(E.angle_m).*D_x{i_iter}+sin(E.angle_m).*D_y{i_iter};
                      jv=-sin(E.angle_m).*D_x{i_iter}+cos(E.angle_m).*D_y{i_iter};



  
                      subplot(3,2,2*i_iter-1)
                      plot(t.*1e15,jp)
                      title([name_arr{i_iter} 'p'])
                      xlabel('time (fs)')

           
                    subplot(3,2,2*i_iter)
                   plot(t.*1e15,jv)
                      title([name_arr{i_iter} 'v'])
                      xlabel('time (fs)')


	    end


                         if obj.save_flag==1
                           h = gcf; % gcf returns the handle to the current figure
                           savefig(h, [obj.save_str   'J_t_pv.fig']);

                       end
        end
