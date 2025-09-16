function t_windown_Fourier(MF,R,name_str)
t= MF.t(MF.t_start:(MF.t_start+length(R.rho_2)-1));
% t_win=0.5/MF.om;
% t_win_shift=t_win/2;

t_win=0.2/MF.om;
t_win_shift=0.2*t_win;
t_range=4*2*pi/MF.om;
N_t=fix((2*t_range)/t_win_shift);

%the Fourier resolution is the same as the entire long t array
dt=t(2)-t(1);
df=1/dt/length(t);
f_limit=length(t).*df/2;
f_arr=linspace(-f_limit,f_limit,length(t));

N_half=fix(length(t)/2)+1;
N_f=length(f_arr(N_half:end));

t_plot=zeros(1,N_t);


Jt={R.j_no_bath_intra+R.j_no_bath_inter,R.j_thermal_intra+R.j_thermal_inter,...
          R.j_s_intra+R.j_s_inter};
J_name={'no_bath','thermal','sqeeze'};
for i_iter=1:length(Jt)
            J=Jt{i_iter};
            scan_mat=zeros(N_f,N_t);
            for scan_iter=1:N_t
  
                t_location=t(abs(J)==max(abs(J)))-t_range+scan_iter.*t_win_shift;
                t_plot(scan_iter)=t_location;
                N_win=fix(t_win/dt);
                win_p=(-fix(N_win/2):1:fix(N_win/2))+fix((t_location-t(1))/dt);
                filter_win=zeros(size(t));
                %-------------------------
                %blackman window  
                %-------------------------
               % filter_temp=0.42-0.5*cos(2*pi*(t-t_location-t_win/2)/(t_win))+0.08*cos(4*pi*(t-t_location-t_win/2)/(t_win));
                % filter_win(win_p)=filter_temp(win_p);
               %---------------------------------------
               %or a super gaussian window
               %----------------------------------------
                filter_temp=exp(-((t-t_location).^2./(t_win/2).^2).^1);
                filter_win=filter_temp;
                %----------------------------------------
                %hamming window
                %----------------------------------------
                % filter_temp=0.54-0.46*cos(2*pi*(t-t_location-t_win/2)/(t_win));
                % filter_win(win_p)=filter_temp(win_p);
            
            
                temp=fftshift(fft(ifftshift(J.*filter_win)))./(length(t)*df);
                %normalized the energy to be the same as t domain
                scan_mat(:,scan_iter)=temp(N_half:end);
       
        
            end
                figure
                h=surf(t_plot./1e-15,f_arr(N_half:end)./(MF.om/(2*pi)),(abs(scan_mat)));
     
                grid off
                box on
                xlabel('Time (fs)')
                view(2);%shading interp; 
                colormap('jet')
                %set(gca,'ColorScale','log')
       
                ylabel('Harmonic Order')
                xlim([-2*pi*1.5e15/MF.om+1e15*t(abs(J)==max(abs(J))),1.5e15*t(abs(J)==max(abs(J)))+1e15*2*pi/MF.om])
                ylim([0,60])
                set(h,'LineStyle','none')
                title(['Linear color scale '  J_name{i_iter} ])
                hold on
    
 
                yyaxis right
                plot(t./1e-15,J,'w')
                hold on
                %input norm
                Et=max(J).*MF.Et./max(MF.Et);
                plot(MF.t/1e-15,Et,'w--')
                savefig(gcf,[name_str  J_name{i_iter}  '_tf.fig'])
end



end