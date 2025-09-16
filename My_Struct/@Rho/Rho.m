classdef Rho
properties 

    rho_2;
    rho_2_s;

    J_t;
    int_t1_arr;
    hbar;
    q;
    rho_2_no_heat;
    J_E_t_tau;

    J_exp_squeeze;


    N_rho;
    N_t;
    Norm_k;
    rho_nm;
    rho_k;
    rho_k_no;
    rho_nn;
    j_thermal_intra;
    j_thermal_inter;
    j_no_bath_intra;
    j_no_bath_inter; 
    j_s_k;
    j_t_k;
    j_no_k;
     j_s_intra;
    j_s_inter;
   J_xy_no_t_s;
    J_norm;
    %all functions declear
   magic_operation;
   magic_operation_squeeze;
   squeeze_om_factor;
   squeeze_E;

end
  methods
            function obj=init(obj,J_t,MF,C,kx_select,ky_select,kz_select,s_w,s_E)
              %normalization occurs when the 
              %1---   dkx dky /V_bz
              %2---   delt(B)dbx dby/V_bz
              %the integration of 1*dkxdky an 1*delt(B)dbx dby are equal
              %which is the area of the BZ. For a orthogonal lattix, kx*ky
              %is directly the area. For skewed lattice there is a del(B)
              %but delt(B) and V_bz cancel out and the N_kx,N_ky is just
              %right
              obj.Norm_k=length(kx_select)*length(ky_select)*length(kz_select);

              obj.N_rho= MF.t_end-MF.t_start+1;
              if mod(length(MF.t),2)==0
                  obj.N_t=fix(length(MF.t)/2)+1; 
              else
                  obj.N_t=fix(length(MF.t)/2)+2;
              end

              obj.rho_2=zeros(1,obj.N_rho);
              obj.rho_2_s=zeros(1,obj.N_rho);
              obj.rho_2_no_heat=zeros(1,obj.N_rho);
              obj.J_xy_no_t_s=zeros(6,obj.N_rho);
              obj.rho_nm=0;
              obj.rho_nn=0;
              obj.j_thermal_intra=0;
              obj.j_thermal_inter=0;
              obj.j_s_intra=0;
              obj.j_s_inter=0;
              obj.j_no_bath_intra=0;
              obj.j_no_bath_inter=0;
              %here the 2 is the spin degeneracy!
              obj.J_norm=MF.J_norm*2;
              obj.squeeze_om_factor=s_w;
              obj.squeeze_E=s_E;
       
              obj.J_t=J_t;
              %here for squeezed pulse the exp(J_exp*v^2)=J_t the
              %log(0)=-inf is to make sure the J(t) not needed arrary parts
              %goes to zero
              obj.J_exp_squeeze=zeros(obj.N_rho,obj.N_rho)+log(0);
              obj.J_E_t_tau=diag(ones(1,obj.N_rho).*0.5);
              obj.hbar=C.hbar;
              obj.q=C.q;

          end
         


         
    end
end