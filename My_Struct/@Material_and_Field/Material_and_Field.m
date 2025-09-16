classdef Material_and_Field
    properties
        %--------------------
        %field parameters
        %--------------------
        dt;
        t;
        om;
        dw;
        w;
        tau;
        E0;
        lambda_si;
        t_win_effect;
        t_start;
        t_end;
        t_range;
        Et;
        At;
        type_N;
        %--------------------
        %material
        %--------------------
        kx;ky;kz;
        dkx;
        coe_c;coe_v;coe_a;
        eg;

        ax;ay;az;




        %---------
        %MoS2
        %---------
        a;
        H_k;
        %with d ending means the imported data
        B;
        B_inv;
        Omega_dipole;
        Eps_0;
        Eps_total;
       etax12;
       etay12;
       gauge_choice;
        v_x11;
        v_x12;
        v_x22;

        v_y11;
        v_y12;
        v_y22;

        v_sq;

        dx_mk;
        dy_mk;
        J_norm;
         angle_m;
         %load data
         kx_d;
         ky_d;
         etax_12_d;
         etay_12_d;
         v_x11_d;
         v_x22_d;
          v_x12_d;
         v_y11_d;
         v_y22_d;
         v_y12_d;
         Eps_d;
 %deal with functiono handel
     renew_kt;
    end

    methods 


              function obj=init_material(obj,C,M_name,angle_select)
                 obj.angle_m=angle_select; 
              switch M_name
                case 'ZnO'
                    obj.eg = 0.129.*C.Energy_si; 
                    obj.d0 = 3.64*C.r_si;  
                    % min bandgap for ZnO
                    delx_c = [0.089,-0.0814,-0.0024,-0.0048,-0.003].*C.Energy_si;                                                   % bandwidth along direction x, ZnO
                    delx_v=[-0.093,0.0705,0.02,-0.0012,0.0029].*C.Energy_si;
                    obj.ax = 5.32.*C.r_si; 
                    
                    dely_c = [0.1147,-0.1147].*C.Energy_si;                                                   % bandwidth along direction x, ZnO
                    dely_v=[-0.0307,0.0307]*C.Energy_si;
                    obj.ay = 6.14.*C.r_si; 
                    
                    delz_c = [0.0435,-0.0435].*C.Energy_si;                                                   % bandwidth along direction x, ZnO
                    delz_v=[-0.0059,0.0059]*C.Energy_si;
                    obj.az=9.83.*C.r_si; 
                    N=230;
                    obj.kx=linspace(-pi/obj.ax,pi/obj.ax,N+1);
                    %obj.ky=linspace(-pi/obj.ay,pi/obj.ay,201)';

                    obj.ky=linspace(0,pi/obj.ay,fix(N/4))';
                    %obj.kz=[0,0.5*pi/obj.az,pi/obj.az];
                    obj.kz=linspace(0,pi/obj.az,N/2);
                    
                    obj.coe_c={delx_c,dely_c,delz_c};
                    obj.coe_v={delx_v,dely_v,delz_v};

                     obj.renew_kt=@ renew_kt_ZnO;
                     %the 4 is to count for the full BZ because 
                    %ky and kz and only take possitive part.
                     obj.J_norm=4*(obj.kx(2)-obj.kx(1))*(obj.ky(2)-obj.ky(1))*(obj.kz(2)-obj.kz(1))/((2*pi)^3);
                                                                  
                case 'MoS2'

              obj.a=3.19e-10;
             eps1=1.046*abs(C.e);
             eps2=2.104*abs(C.e);
            t0=-0.184*abs(C.e);
            t1=0.401*abs(C.e);
            t2=0.507*abs(C.e);
            t11=0.218*abs(C.e);
            t12=0.338*abs(C.e);
            t22=0.057*abs(C.e);

            h0=@(kx,ky) 2.*t0.*(cos(obj.a.*kx)+2.*cos(0.5*kx.*obj.a).*cos(0.5*sqrt(3).*ky.*obj.a))+eps1;
            h1=@(kx,ky) -2*sqrt(3)*t2.*sin(0.5.*obj.a*kx).*sin(0.5*sqrt(3).*ky.*obj.a)...
                +2i*t1.*(sin(kx.*obj.a)+sin(0.5*obj.a.*kx).*cos(0.5*sqrt(3).*ky.*obj.a));
            h2=@(kx,ky) 2*t2.*(cos(obj.a.*kx)-cos(0.5.*obj.a.*kx).*cos(0.5*sqrt(3).*ky.*obj.a))...
               +2*sqrt(3)*1i*t1.*cos(0.5.*obj.a.*kx).*sin(0.5*sqrt(3).*ky.*obj.a);
            h11=@(kx,ky) 2*t11.*cos(obj.a.*kx)+(t11+3*t22).*cos(0.5.*obj.a.*kx).*cos(0.5*sqrt(3).*ky.*obj.a)+eps2;
            h22=@(kx,ky) 2*t22.*cos(obj.a.*kx)+(t11*3+t22).*cos(0.5*obj.a.*kx).*cos(0.5*sqrt(3).*ky.*obj.a)+eps2;
            h12=@(kx,ky) sqrt(3).*(t22-t11).*sin(0.5.*obj.a.*kx).*sin(0.5*sqrt(3).*ky.*obj.a)...
                +4i*t12.*sin(0.5.*obj.a.*kx).*(cos(0.5.*kx.*obj.a)...
               -cos(0.5*sqrt(3).*ky.*obj.a));

               obj.H_k=@(kx,ky) [h0(kx,ky),h1(kx,ky),h2(kx,ky);
                conj(h1(kx,ky)),h11(kx,ky),h12(kx,ky);
                  conj(h2(kx,ky)),conj(h12(kx,ky)),h22(kx,ky)];
               obj.gauge_choice=[3,2,1];
            % 

                    %-----------------------------  
                    %graphene
                     %-----------------------------  
                      %  obj.a=obj.a/sqrt(3);
                      %   h1=@(kx,ky)  (C.hbar.*1e6/obj.a).*(exp(1i*ky*obj.a/2)*2.*cos(sqrt(3)*kx*obj.a/2)+exp(-1i*ky*obj.a));
                      %   h0=@(kx,ky)  (C.hbar.*1e6/obj.a).*0.1;
                      %  obj.H_k=@(kx,ky) [-h0(kx,ky),h1(kx,ky);
                      %                       conj(h1(kx,ky))  ,h0(kx,ky)];
                      % obj.gauge_choice=[1,1];

                    %-----------------------------        
                    %define the mesh for either Monkhorst or K space
                    %-----------------------------  
                    %    ------%Monkhorst
                      N=251;
                     d_win=2*pi/(sqrt(3)*obj.a);
                      obj.kx=linspace(-d_win,d_win,N+1);
                       obj.ky=obj.kx;
             
                     %name_head='Monkhorstg321d300';
                      obj.B= [sqrt(3)/2,-sqrt(3)/2;1/2,1/2];  
                      obj.B_inv=inv(obj.B);
                    obj.J_norm=det(obj.B)*(obj.kx(2)-obj.kx(1))*(obj.ky(2)-obj.ky(1))/((2*pi)^2);
                       %-----------------------------  
                       %    ------%K space
                        % N=500;
                        % d_winx=1.5*(4*pi/obj.a/3);
                        % d_winy=4*pi/(sqrt(3)*obj.a);
                        % obj.kx=linspace(-d_winx,d_winx,N+1);
                        % obj.ky=linspace(-d_winy,d_winy,N+1);
                        % obj.B= [1,0;0,1];  
                        % obj.B_inv=inv(obj.B);
                        % 
                        % obj.J_norm=(1/4)*(obj.kx(2)-obj.kx(1))*(obj.ky(2)-obj.ky(1))/((2*pi)^2);
                  %   %-----------------------------  
                  %   %Monkhorst analytic
                  %    %-----------------------------  

                       obj.kx=obj.kx(1:end-1);
                       obj.ky=obj.ky(1:end-1);
                       obj.kz=0;
                       obj.dx_mk=obj.kx(2)-obj.kx(1);
                       obj.dy_mk=obj.ky(2)-obj.ky(1);
                       obj.v_x11=zeros(size(obj.t_range));
                       obj.v_x12=zeros(size(obj.t_range));
                       obj.v_x22=zeros(size(obj.t_range));
                       obj.v_y11=zeros(size(obj.t_range));
                       obj.v_y12=zeros(size(obj.t_range));
                       obj.v_y22=zeros(size(obj.t_range));

                       obj.Omega_dipole=zeros(size(obj.t_range));
                       obj.Eps_0=zeros(size(obj.t_range));
               
                       obj.renew_kt=@ renew_kt_MoS2_analytic;

                                        
                  %   %-----------------------------  
                  %   %Monkhorst load data
                  %    %-----------------------------  
                    % 
                    %                    loop_d={'kx_d','ky_d','etax_12_d','etay_12_d', 'v_x11_d','v_x22_d','v_x12_d',...
                    %      'v_y11_d','v_y22_d','v_y12_d' };
                    % 
                    % loop_text={'kx','ky',['etax12' ],'etay12', ['vx11'],['vx22' ],['vx12' ],...
                    %   ['vy11'  ],['vy22' ],['vy12']};
                    % 
                    %    name_head='Monkhorstg321300';
                    %  for i_iter=1:length(loop_d)
                    % 
                    %              obj.(loop_d{i_iter})=readmatrix([name_head loop_text{i_iter} '.txt']);
                    %  end
                    %  obj.kx=linspace(obj.kx_d(1),obj.kx_d(end),N+1);
                    %   obj.ky=obj.kx;
                    %   obj.kx=obj.kx(1:end-1);
                    %   obj.ky=obj.ky(1:end-1);
                    %   obj.kz=0;
                    %   obj.dx_mk=obj.kx(2)-obj.kx(1);
                    %   obj.dy_mk=obj.ky(2)-obj.ky(1);
                    % 
                    % 
                    % 
                    %   E1=readmatrix([name_head 'E1.txt']);
                    %   E2=readmatrix([name_head 'E2.txt']);
                    %  obj.Eps_d= E2-E1;
                    % 
                    % obj.renew_kt=@ renew_kt_MoS2_load;




          
         end


end



    end 
end