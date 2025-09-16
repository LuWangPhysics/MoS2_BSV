%Lulu's best code ever!
clear all;
restoredefaultpath;
addpath('My_Struct');
addpath('My_Functions');
addpath('My_Output');
addpath('My_Data');

	    % ========================== define constants==============================
	    C=Const;
	    C=C.init();
	    % ========================== material and driving laser parameters ========
	    MF=Material_and_Field;
	    M_name='MoS2';
	    E0=5e8;


	    angle_arr=[0*pi/3];
	    order_select=[2,3,4,5,6];


	    MF=MF.init_field(C,E0);
	    MF=MF.init_material(C,M_name,angle_arr);


	    % ========================== set-up heat bath general T numerical =========
	    C.T=300;                                  %heat bath temperature K
	    jo=2;         				              %heat bath coupling strength
	    om_cutoff =200*(2*pi*1e12);               %heat bath cutoff frequency                                 
	    s_iter=1;  				                  %BSV center frequency
	    s_E_arr=1;   			 	              %BSV energy      


        %analytical_heat_bath_T(C,MF,om_cutoff,jo);
        % ========================== set-up heat bath general T numerical =========
        %choose name among 'Dybe','OM','Gauss' 'Shift_Gauss','under_damp'
        name_select='OM';
        [Jt_FT,Jt_int]=numerical_general_heat_bath(C,MF,om_cutoff,jo,name_select);
        save_str=[name_select 'd250_BSV_r_together' M_name 'SV1e7_Sw' num2str(s_iter) 'om_SE'  num2str(s_E_arr) '_E' num2str(E0/1e9) ...
                  'Vpernm_angle' num2str((180*angle_arr/pi))  ];
        P=My_plot;
        P=P.save(save_str);
 

        %--------------------
        %all K domains  
        %--------------------
         kx_select=1:length(MF.kx);
         ky_select= 1:length(MF.ky);
  


         R=Rho;
         R=R.init(Jt_FT,MF,C,MF.ky,MF.ky,MF.kz,s_iter,s_E_arr);
         R=R.Jt_construct(C,MF);    

         N_t=length(MF.t_start:MF.t_end);%17;
         j_kxky_no=zeros(length(MF.kx),length(MF.ky),N_t);       
         j_kxky_t=zeros(length(MF.kx),length(MF.ky),N_t);       
         j_kxky_s=zeros(length(MF.kx),length(MF.ky),N_t);       



for ky_iter=ky_select
        for kx_iter=kx_select
        %----------------------------------------
        %renew the time dependent material
        %----------------------------------------
        MF=MF.renew_kt(MF,C,kx_iter,ky_iter,1);

        %----------------------------------------
        %calculate the ionization for a given k
        %----------------------------------------
        R=R.rho_update_thermal(MF,C,jo);
        R=R.rho_update_squeeze(MF,C);
        %only store the x direction
        j_kxky_no(kx_iter,ky_iter,:)=R.j_no_k;
        j_kxky_t(kx_iter,ky_iter,:)=R.j_t_k;
        j_kxky_s(kx_iter,ky_iter,:)=R.j_s_k;
 
        end
end




P.plot_rho2(MF,R,C.T)
P.plot_J_total_xy(MF,R,C.T,order_select);
P.plot_J_total_pv(MF,R,C.T,order_select);
P.j_kxky(MF,j_kxky_no,j_kxky_t,j_kxky_s);




