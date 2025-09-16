function plot_J_total_pv_t_given_order(obj,E,R,T, n_order)
            t=E.t(E.t_start:(E.t_start+length(R.rho_2)-1));
            D_x={R.J_xy_no_t_s(1,:),R.J_xy_no_t_s(3,:),...
                      R.J_xy_no_t_s(5,:)};
            D_y={R.J_xy_no_t_s(2,:),R.J_xy_no_t_s(4,:),...
                      R.J_xy_no_t_s(6,:)};

            name_arr={'no heat bath', 'thermal ', 'squeeze'}; 
             %save as non pv, thermal pv, sq pv
             dt=t(2)-t(1);
             f=linspace(-1/(2*dt),1/(2*dt),length(t));
             df=1/dt/length(t);
             N_order=2;
             f_win=-fix(0.5*E.om/(2*pi*df)):1:fix(0.5*E.om/(2*pi*df));

 for i_iter=1:length(D_x)
           figure
                 for i_order=1:N_order
                      jp=cos(E.angle_m).*D_x{i_iter}+sin(E.angle_m).*D_y{i_iter};
                      jv=-sin(E.angle_m).*D_x{i_iter}+cos(E.angle_m).*D_y{i_iter};
                       f_location_o=f_win+fix(length(f)/2)+fix((2*i_order+1).*E.om/(2*pi*df));
                       f_location_e=f_win+fix(length(f)/2)+fix((2*i_order).*E.om/(2*pi*df));
  

                      
                      subplot(2,N_order,i_order)
                      temp=(fftshift(fft(ifftshift(jp)))./(length(t)*df));
                      jp_f=zeros(size(temp));
                      jp_f(f_location_o)=temp(f_location_o);
                      jp_t=fftshift(ifft(ifftshift(jp_f))).*length(t)*df;
                      plot(t.*1e15,real(jp_t))
                      title([name_arr{i_iter} 'p' num2str(i_order*2+1)])
                      xlabel('t (fs)')
                      ylabel('E(t)')
                     
             
           
                      subplot(2,N_order,N_order+i_order)
                      temp=(fftshift(fft(ifftshift(jv)))./(length(t)*df));

                      jv_f=zeros(size(temp));
                      jv_f(f_location_e)=temp(f_location_e);
                      jv_t=fftshift(ifft(ifftshift(jv_f))).*length(t)*df;
                      plot(t.*1e15,real(jv_t))
                      title([name_arr{i_iter} 'v' num2str(i_order*2)])
                      xlabel('t (fs)')
                      ylabel('E(t)')
                     
                     
                 end
                     if obj.save_flag==1
                           h = gcf; % gcf returns the handle to the current figure
                           savefig(h, [obj.save_str   'J_spec_pv_t_given_orderabs' name_arr{i_iter}  '.fig']);
                       end
 end



end
