 function output=manipulate_long_squeeze(obj,C,f1_HI,MF)

                %calculate matrix multiplication for each sub sections
                    output=zeros(size(f1_HI));
                  % output(1)=0*f1_HI(1)*1/2;
                        for iter_tau=MF.t_start+1:MF.t_end-1
                          temp1=obj.numerical_squeeze(C,MF,MF.t(iter_tau),MF.t(MF.t_start:iter_tau));
                          temp2=flip([1/2,temp1]);
                          output(iter_tau-MF.t_start+1)=sum(f1_HI(1:length(temp2)).*temp2).*MF.dt;
                        end


          end