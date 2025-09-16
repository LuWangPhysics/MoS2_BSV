function save_data(name_select,om_cutoff,MF,rho_scan,T_arr,jo_arr)
save_str=['My_Output/' name_select 'cutoff_freq' num2str(om_cutoff/(2*pi*1e12)) 'THz_E0' num2str(MF.E0/1e9) 'V_pernm'];
save([save_str+".mat"],"T_arr","jo_arr","rho_scan","-v7.3");

if(length(jo_arr)>1)
h=figure
fig=surf(T_arr,jo_arr,log10(rho_scan'))
set(fig,'LineStyle','none')
colormap(jet)
xlabel('Temperrature (K)')
ylabel('coupling coefficient jo')
zlabel('rho')
savefig(h,[save_str '.fig'])
end
end