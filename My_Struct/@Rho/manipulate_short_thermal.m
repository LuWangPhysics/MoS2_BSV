          function output=manipulate_short_thermal(obj,f1_HI,MF)
             output =(f1_HI*(obj.J_E_t_tau)).*MF.dt;
          end