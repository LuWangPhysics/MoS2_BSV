function [exp_Jt_ft,exp_Jt_int]=numerical_general_heat_bath(C,E,om_cutoff,jo,name_select)

bet = 1/(C.Kb*C.T);                                                        % 1/(K*bT) 
N_half=fix(length(E.w)/2);                                
% ========================== define different bath type ===================%
 if strcmp(name_select,'T2')
     T2=5e-15;%2*pi/(E.om*2);%bet*C.hbar/(2*jo*pi);
     exp_Jt_int=exp(-E.t./T2);%.*exp(1i*jo*pi.*exp(-om_cutoff.*E.t));
     exp_Jt_int(1:N_half)=flip(exp_Jt_int((N_half+2):end));
    %exp_Jt_int(N_half+1)=1;
     exp_Jt_ft=exp_Jt_int;
 else
        switch name_select
            case 'Dybe'
                w_J =jo* E.w./(E.w.^2+om_cutoff^2);                                  %Debye spectral density referred to J(omega)
            case 'Shift_Dybe'
                w_J =jo* E.w./((E.w.^2-E.om^2)+om_cutoff^2);   
            case 'OM'
                w_J= (jo*E.w./om_cutoff^2).*exp(-abs(E.w)./om_cutoff);                      %Ohm spectral density
            case 'Gauss'
                w_J=(jo.*E.w./om_cutoff^2).*exp(-E.w.^2./om_cutoff^2);                      %Gaussian spectral density
            case 'Shift_Gauss'
        %-----------------
        %shift gausian
        %-----------------
                omw = om_cutoff/2;
                env=exp(-(E.w-om_cutoff).^2/omw^2);
                env(1:N_half)=flip(env((N_half+2):end));
                w_J=(jo.* E.w/omw^2).*env;  
            case 'under_damp'
        
                w_J=jo.*om_cutoff.^2.*E.w./(om_cutoff^2.*E.w.^2+(E.w.^2-E.om^2).^2);
        
        end        
         
        
        % ========================== integration calculation from J(w) to J(t) ====
        
        w_range=1:length(E.w);%1:length(E.w);
        %-------------------------------------
        %Fourier method for the integral
        %-------------------------------------
        if C.T==0
            jdeb1 = w_J.*sign(E.w);
        else
            coth_fix=coth(C.hbar.*bet*E.w/2);
            coth_fix(isinf(coth_fix))=coth_fix(fix(length(E.w)/2+2));
            jdeb1 = coth_fix.*w_J;                            
        end
        
        I0d = -trapz(jdeb1(w_range))*E.dw;                                 % time independent term in exponent which is part of C_{\beta}; (1/2) --> w runs over all frequencies; om only over positive frequencies. 
        Icd = fftshift(real(ifft(ifftshift(jdeb1(w_range)))))*2*pi/E.dt;   % noise exponent Debye real; time dependent part of correlation function in the exponant calculated by IFT (e^{iwt}) 
        Isd = fftshift(imag(ifft(ifftshift(w_J(w_range)))))*2*pi/E.dt;     % noise exponent Debye imag

       exp_Jt_ft=exp(I0d+Icd+1i.*Isd);
       
          % Isd_ana=jo*4*om_cutoff.*E.t.*(1-0*4*om_cutoff^2.*E.t.^2);
          % exp_Jt_ft=exp(0*I0d+0*Icd+1i.*Isd_ana);
        %-------------------------------------
        %numerial mesh method for the integral
        %--------------------------------------
        if E.type_N>length(E.t)
            if C.T==0
                C_wt=(-(1-cos(E.w'*E.t)).*sign(E.w')+1i.*sin(E.w'*E.t)).*w_J'; 
            else 
                C_wt=(-(1-cos(E.w'*E.t)).*coth_fix'+1i.*sin(E.w'*E.t)).*w_J';
            end
            
            C_t=sum(C_wt(w_range,:),1).*E.dw;
            exp_Jt_int=exp(C_t);
        else
           % C_t=0;
           % N_num=300;
           % N_sec=fix(length(E.w)/N_num);
           % w_start=(0:N_num:N_sec*N_num)+1;
           % w_end=[w_start(2:end)-1,length(E.w)];
           % for sec_iter=1:N_sec
           %      w_iter=w_start(sec_iter):w_end(sec_iter);
           %          if C.T==0
           %              C_wt=(-(1-cos(E.w(w_iter)'*E.t))+1i.*sin(E.w(w_iter)'*E.t)).*w_J(w_iter)'; 
           %          else 
           %              C_wt=(-(1-cos(E.w(w_iter)'*E.t)).*coth_fix(w_iter)'+1i.*sin(E.w(w_iter)'*E.t)).*w_J(w_iter)';
           %          end
           %        C_t=C_t+E.dw*sum(C_wt,1);
           % 
           % end
            exp_Jt_int=0.*exp_Jt_ft;
        end
end

end