function analytical_heat_bath_T(C,MF,om_cutoff,jo)
bet = 1/(C.Kb*C.T);  
temp1=C.hbar*om_cutoff*bet;
temp_l=1+2*exp(-MF.w.*C.hbar.*bet);
%Debye
J_t_h_Debye=jo*pi*exp(-om_cutoff.*abs(MF.t)).*(2/temp1+1i*sign(MF.t))-2*pi*jo/temp1;
J_t_l_Debye=jo*pi*exp(-om_cutoff.*abs(MF.t)).*1i.*sign(MF.t);
%Ome
J_t_h_OM=jo.*(4/temp1+4i*MF.t*om_cutoff./(1+om_cutoff^2.*MF.t.^2))./(1+om_cutoff^2.*MF.t.^2)...
    -4*jo/temp1;
J_t_l_OM=4i*MF.t*om_cutoff./(1+om_cutoff^2.*MF.t.^2).^2;
%Gaussian
J_t_h_Gauss=sqrt(pi).*jo.*(-2/temp1+(2/temp1+1i*om_cutoff.*MF.t./2).*exp(-MF.t.^2.*om_cutoff^2./4));
J_t_l_Gauss=1i*jo*sqrt(pi)*om_cutoff.*MF.t./2;
%shift Gaussian

J_t={J_t_h_Debye,J_t_l_Debye,J_t_h_OM,J_t_l_OM,J_t_h_Gauss,J_t_l_Gauss};
name_str={'Debye High T','Debye Low T','OM High T','OM Low T','Gauss High T','Gauss Low T'};
figure('Name','Heat bath analytical','NumberTitle','off')
for i_iter=1:length(J_t)
subplot(3,2,i_iter)
plot(MF.t./1e-15,[real(exp(J_t{i_iter}));imag(exp(J_t{i_iter}))],'LineWidth',1)
legend('Re','Im')
xlabel('time (fs)')
title( name_str{i_iter})
xlim([0,5])
end
end