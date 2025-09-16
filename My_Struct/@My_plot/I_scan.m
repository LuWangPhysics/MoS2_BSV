function I_scan(obj,E0_arr,C,I_store)
figure
plot(E0_arr.^2*0.5*C.eps0*C.c/1e16,I_store)
xlabel('intensity (TW/cm^2)')
ylabel('|J(w)|^2')
set(gca, 'YScale', 'log')

figure
plot((E0_arr./1e9),I_store)
xlabel('E (V/nm)')
ylabel('|J(w)|^2')
set(gca, 'YScale', 'log')

                  if obj.save_flag==1
                           h = gcf; % gcf returns the handle to the current figure
                           savefig(h, [obj.save_str  'I_scan.fig']);
                  end
end