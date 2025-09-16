classdef My_plot
    properties
              save_str;
              save_flag;
    end
    methods
              function    obj=save(obj, save_str)
                        obj.save_str=save_str;
                        if isempty(obj.save_str)
                                 obj.save_flag=0;
                        else
                                  obj.save_flag=1;
                        end
              end
        function plot_E_field(obj,E)


            figure('Name','laser electric field','NumberTitle','off');
            plot(E.t./1e-15,E.Et)
            xlabel('Time (fs)')
            ylabel('E field strengh (V/m)')
        end



        function plot_rho2(obj,E,R,T)
           t=E.t(E.t_start:(E.t_start+length(R.rho_2)-1));
            figure
            plot(t./1e-15,[R.rho_2_no_heat;R.rho_2;R.rho_2_s])
            set(gca, 'YScale', 'log')
            xlim([-4*E.tau*1e15,4*E.tau*1e15])
            legend('no heat bath','heat bath','squeeze')
            xlabel('Time (fs)')
            ylabel('Ionization ')
            title(['Temperature=' num2str(T) 'K'])

                  if obj.save_flag==1
                           h = gcf; % gcf returns the handle to the current figure
                           savefig(h, [obj.save_str,'rho.fig']);
                  end
        end

        function plot_J_inter_intra(obj,E,R,T)
            t=E.t(E.t_start:(E.t_start+length(R.rho_2)-1));
           if (sum(sum(abs(R.j_no_bath_intra)))~=0)
              D_arr={R.j_thermal_intra,R.j_thermal_inter,R.j_no_bath_intra,R.j_no_bath_inter,...
                      R.j_s_intra,R.j_s_inter};
            name_arr={'heatbath-intra',' heatbath-inter','no-heat-intra',' no-heat-inter', 'squeeze-intra','squeeze-inter'}; 
            else
            D_arr={R.j_s_intra,R.j_s_inter};
            name_arr={ 'squeeze-intra','squeeze-inter'}; 
            end
     
            for i_iter=1:length(D_arr)
                       if mod(i_iter,2)~=0
                                 figure
                                 p_iter=1;
                       end
                      subplot(2,2,(p_iter-1)*2+1)
                      plot(t./1e-15,D_arr{i_iter})
                      %set(gca, 'YScale', 'log')
                      xlim([-4*E.tau*1e15,4*E.tau*1e15])
                      xlabel('Time (fs)')
                      ylabel('Current J')
                      title([name_arr{i_iter}])
                      
                      subplot(2,2,p_iter*2)
                      dt=t(2)-t(1);
                      f=linspace(-1/(2*dt),1/(2*dt),length(t));
                      df=1/dt/length(t);
                      jf=fftshift(fft(ifftshift(D_arr{i_iter})))./(length(t)*df);
                      plot(f./(E.om/(2*pi)),abs(jf))
                      set(gca, 'YScale', 'log')
                      xlabel('Harmonics orders')
                      xlim([0,20])
                      p_iter=p_iter+1;

                  if obj.save_flag==1&&(mod(i_iter,2)==0)
                           h = gcf; % gcf returns the handle to the current figure
                           savefig(h, [obj.save_str  name_arr{i_iter} '.fig']);
                  end
	    end

                              
        end




                 

        function J_t_int_FT_compare(obj,Jt_FT,Jt_int,MF,name_str,C,j0)
             name_select=[name_str,'_T=',num2str(C.T),'_Jo=' num2str(j0)];
            N_half=fix(length(MF.t)/2);
            plot_range =(N_half):(N_half+MF.t_win_effect);    
            figure('Name',name_select,'NumberTitle','off');
            subplot(1,2,1)
            plot(MF.t(plot_range)./1e-15,[real(Jt_FT(plot_range));imag(Jt_FT(plot_range))],'LineWidth',1,'Marker','square')
            xlim([0,5])
            xlabel('fs')
            ylim([-1e-3,1])
            title('numerical FFT')

            subplot(1,2,2)
            plot(MF.t(plot_range)./1e-15,[real(Jt_int(plot_range));imag(Jt_int(plot_range))],'LineWidth',1,'Marker','square')
            xlim([0,5])
               ylim([-1e-3,1])
                xlabel('fs')
            title('numerical integration')
    

        end
    end 
end