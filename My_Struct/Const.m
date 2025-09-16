classdef Const
    properties
        
        c=299792458;
        e=1.602176634e-19;
        q=-1.602176634e-19;
        hbar=1.05457181e-34;
        Kb=1.380649e-23; 
        E0_si=5.14*1e11;
        r_si=0.528e-10;
        m=9.1093837e-31;
        eps0=8.8541878188e-12;
        Energy_si;
        omega_si;
        time_si;
        T;

    end

    methods
        function obj=init(obj)
              obj.Energy_si=27.211386245981*obj.e;
              obj.omega_si=obj.Energy_si/obj.hbar;
              obj.time_si=1/obj.omega_si;

        end
  end
end