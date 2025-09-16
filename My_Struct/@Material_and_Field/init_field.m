function obj=init_field(obj,C,E0)
      
        obj.E0 = E0;                                                    % peak electric field strength     

        obj.lambda_si=3.2e-6;%3200e-9;   
        obj.om = C.c*(2*pi/obj.lambda_si);                                 % circular frequency; 0.057 [at.u.] approx 800nm

        %20fs for 32um   

        obj.tau=40e-15; %20e-15;                                                    %pulse duration 

       

        %for fine mesh

              % N_mesh=13;
              % t_lim=13*obj.tau;
               % t_lim=9*obj.tau;
               % N_mesh=3;
           

        % %%for loose mesh
         t_lim=5.8*obj.tau;
         N_mesh=1.1;
        %the value decide to get which calculation methods
        obj.type_N=40000;
  
        obj.t = linspace(-t_lim,t_lim,N_mesh*8e3+1);                       % time vector                                          
        obj.dt = obj.t(2)-obj.t(1);    
        
        %this t_win has to be smller than 0.5*t_lim
        t_win=2.8*obj.tau;
        obj.t_win_effect=fix(t_win/obj.dt);                          %effetive time windwo for ionizaiton 
        
        obj.dw=2*pi/obj.dt/length(obj.t);
        w_limit=length(obj.t).*obj.dw/2;
        obj.w=linspace(-w_limit,w_limit,length(obj.t));
        
        obj.Et = obj.E0*exp(-(obj.t/obj.tau).^2).*cos(obj.om*obj.t);       % field function of time
        obj.At = -cumtrapz(obj.Et)*obj.dt;                                   % vector potential 


        obj.t_start=fix(length(obj.t)/2)-obj.t_win_effect;
        obj.t_end=fix(length(obj.t)/2)+obj.t_win_effect;
        obj.t_range=obj.t_start:obj.t_end;
        


 end 

