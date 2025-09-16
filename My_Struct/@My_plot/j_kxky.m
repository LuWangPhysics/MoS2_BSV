function j_kxky(obj,MF,j_data1,j_data2,j_data3)

                  if (sum(sum(abs(j_data1)))~=0)
                  str_name={'noEnv','thermal','squeeze'};
                  data={j_data1,j_data2,j_data3};
                  else
                    str_name={'squeeze'};
                  data={j_data3};
                  end
                  i_d=length(str_name);
                  Nt=length(j_data1(1,1,:));
                  %------------------------------------------------------------------
                  %find the maximum slice and minimum slice and plot the
                  %difference
                   %------------------------------------------------------------------
                 %  max_store=zeros(3,Nt);
                 % max_p=zeros(1,3);
                 % min_p=zeros(1,3);
                 %  for i_iter=1:i_d
                 % 
                 %         for t_iter=1:Nt
                 %                     max_store(i_iter,t_iter)=max(max(abs(data{i_iter}(:,:,t_iter))));
                 %         end
                 %      max_p(i_iter)=find(max_store(i_iter,:)==max(max_store(i_iter,:)));
                 %      min_p(i_iter)=find(max_store(i_iter,:)==min(max_store(i_iter,:)));
                 %  end
                 % 
                 % 
                 % 
                 %  figure      
                 %  for iter=1:i_d
                 %            subplot(2,i_d,iter)
                 %            imagesc(MF.kx./(pi/MF.a),MF.ky./(pi/MF.a),data{iter}(:,:,max_p(iter))'-data{iter}(:,:,min_p(iter))')
                 %             xlabel('kx/(pi/a)')
                 %            ylabel('ky/(pi/a)')
                 %            set(gca,'YDir','normal')
                 %            view(2)
                 %            colormap(jet)
                 %            title([str_name{iter} 'max-min contrast'])
                 %            subplot(2,i_d,iter+3)
                 %            imagesc(MF.kx./(pi/MF.a),MF.ky./(pi/MF.a),data{iter}(:,:,max_p(iter))')
                 %            xlabel('kx/(pi/a)')
                 %            ylabel('ky/(pi/a)')
                 %            set(gca,'YDir','normal')
                 %            view(2)
                 %            colormap(jet)
                 %             title([str_name{iter} 'max'])
                 %  end
 
                 %---------------------------------------
                 %spectrum in kx,ky distribution
                 %---------------------------------------
                      t=MF.t(MF.t_start:MF.t_end);
                      dt=t(2)-t(1);
                      f=linspace(-1/(2*dt),1/(2*dt),length(t));
                      df=1/dt/length(t);
                      order_select=[3,4,5,6];
                      f_location=fix(length(f)/2)+fix(order_select.*MF.om/(2*pi*df));
                      surf_temp=zeros(length(MF.kx),length(MF.ky),length(f_location)); 
                        
for d_iter=1:i_d
                     figure  
                        for x_iter=1:length(MF.kx)
                                  for y_iter=1:length(MF.ky)

                                          jf=fftshift(fft(ifftshift(data{d_iter}(x_iter,y_iter,:))))./(length(t)*df);                                        
                                          surf_temp(x_iter,y_iter,:)=(abs(jf(f_location)));
                                  end 
                        end
                        for p_iter=1:length(f_location)
                                  subplot(1,length(f_location),p_iter)
                                  imagesc(MF.kx./(pi/MF.a),MF.ky./(pi/MF.a),transpose(surf_temp(:,:,p_iter)))
                               % imagesc(transpose(surf_temp(:,:,p_iter)))
                                 xlabel('kx/(pi/a)')
                                 ylabel('ky/(pi/a)')
                                 set(gca,'ColorScale','log')
                                 set(gca,'YDir','normal')
                                 colormap(jet)
                                 title([str_name{d_iter} 'harmonic' num2str(order_select(p_iter))])
                        end
                        if obj.save_flag==1
                           h = gcf; % gcf returns the handle to the current figure
                           savefig(h, [obj.save_str str_name{d_iter} 'kxky.fig']);
                        end
end
                      
                      
         
        end