function v=myatan_0Topi(theta_r,theta_i)

v=zeros(size(theta_r));
%range [-pi,pi]
for t_iter=1:length(v)
          x=theta_r(t_iter);
          y=theta_i(t_iter);
          if  (y>=0 && x>=0)
              v(t_iter)=mod(atan(y/x),pi/2);
          end
          if (y>=0 && x<0)
              v(t_iter)=pi+atan(y/x);
          end
           if (y<0 && x>=0)
              v(t_iter)=atan(y/x);
          end
          if (y<0 && x<=0)
              v(t_iter)=-pi+atan(y/x);
          end


end

end