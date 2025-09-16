       function output=manipulate_short_squeeze(obj,C,f1_HI,MF)

             temp=exp(obj.J_exp_squeeze.*MF.v_sq.^2);
             %here the log(0)*0 will give NAN, thus remove it manully
             temp(isnan(temp))=1;
             
             temp(1:(obj.N_rho+1):end)=0.5;
             output =(f1_HI*temp).*MF.dt;
          end