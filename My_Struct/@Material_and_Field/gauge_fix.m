function eigen_vec=gauge_fix(obj,eigen_vec)                     
   [N_ind,~]=size(eigen_vec);
   for i_iter=1:N_ind
                phase_fix=angle(eigen_vec(obj.gauge_choice(i_iter),i_iter));
                eigen_vec(:,i_iter)=exp(-1i.*phase_fix).*eigen_vec(:,i_iter);
   end

end