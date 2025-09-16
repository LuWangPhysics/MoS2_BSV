function sum_t=bra_ket_sum(obj,phi_l,phi_r,H)
sum_t=0;
dim=length(phi_l);
for x_iter=1:dim
          for y_iter=1:dim
          sum_t=sum_t+conj(phi_l(x_iter))*H(x_iter,y_iter)*phi_r(y_iter); 
          end 
end
end