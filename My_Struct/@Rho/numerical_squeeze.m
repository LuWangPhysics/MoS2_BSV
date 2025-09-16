function re_exp=numerical_squeeze(obj,C,MF,t1_num,t2_arr)
w_0=MF.om*obj.squeeze_om_factor;
E0=obj.squeeze_E*3e7;
tau=MF.tau;
Ew0=tau*E0/sqrt(2);
beam_size=pi*(4.5e-4)^2; %  [m^2]
U=beam_size*C.c*C.eps0*E0^2*tau*sqrt(pi/8);

%-------------------------
%for multi frequencies
%-------------------------
%dt=t(2)-t(1);
% df=1/(dt*length(t));
% f_limit=1/(2*dt);
% f=linspace(-f_limit,f_limit,length(t));

bw_limit=2*2/tau;%f_limit*2*pi;
w_arr=w_0+linspace(-bw_limit,bw_limit,length(MF.t));
dw=w_arr(2)-w_arr(1);

if(w_arr(1)<0)
          error('squeezed vacumm frequency overload')
end
E_w_arr=Ew0*exp(-tau^2.*(w_arr-w_0).^2./4);

N_bin=151;
N_bin_size=fix(length(w_arr)/N_bin);
Energy_bin=zeros(1,N_bin);
wc_bin=zeros(1,N_bin);
for n_iter=1:N_bin
          w_iter=(n_iter-1)*N_bin_size+1;
          Energy_bin(n_iter)=dw*beam_size*(0.5*C.c*C.eps0)*sum(E_w_arr(w_iter:(w_iter+N_bin_size)).^2);
          wc_bin(n_iter)=w_arr(w_iter+fix(N_bin_size/2));
end
dw_bin=wc_bin(2)-wc_bin(1);

cosh_2r=1+2.*Energy_bin./(C.hbar.*wc_bin);
sinh_2r=sqrt(cosh_2r.^2-1);

theta_s=0.5*pi.*ones(size(cosh_2r));



J_im=zeros(size(t2_arr));
J_re=zeros(size(t2_arr));
for iter_i=1:N_bin
    w_c=wc_bin(iter_i);
    %E_s_delta=sqrt(cosh_2r(iter_i)-sinh_2r(iter_i)*cos(2.*chi-theta_s_arr(iter_i)));

    J_w_c=C.e^2/(16*pi^2*C.c^3*C.hbar*w_c*C.eps0);

    %----------------
    %imaginary part
    %----------------
    J_im=J_im+dw_bin.*J_w_c.*(sin(w_c.*(t1_num-t2_arr))+sin(w_c.*t2_arr)-sin(w_c.*t1_num));

    %----------------
    %real part
    %----------------
    re_1=1-cos(w_c.*(t1_num-t2_arr));
    re_2=(sin(w_c.*t1_num)-sin(w_c.*t2_arr))./(cos(w_c.*t1_num)- cos(w_c.*t2_arr));
    re_2(isnan(re_2))=0;
    re_3=cosh_2r(iter_i)+cos(2.*atan(re_2)-theta_s(iter_i)).*sinh_2r(iter_i);

    J_re=J_re+dw_bin.*J_w_c.*re_1.*re_3;

end

%J_t=exp(MF.v_diff_x^2*(1i.*J_im-J_re));
re_exp=-J_re+1i.*J_im;


end
