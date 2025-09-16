 function output=manipulate_long_thermal(obj,f1_HI,MF)

                %calculate matrix multiplication for each sub sections
                    output=zeros(size(f1_HI));
                  % output(1)=0*f1_HI(1)*1/2;
                        for iter_tau=MF.t_start+1:MF.t_end-1
                          temp1=obj.J_t(obj.N_t:(obj.N_t+iter_tau-MF.t_start-1));
                          temp2=flip([1/2,temp1]);
                          output(iter_tau-MF.t_start+1)=sum(f1_HI(1:length(temp2)).*temp2).*MF.dt;
                        end


          end