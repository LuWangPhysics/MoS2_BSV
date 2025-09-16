function obj=Jt_construct(obj,C,MF)   


            if length(MF.t)<MF.type_N
                              
                                 % -------------------------------------------------------

                                  for iter_tau=MF.t_start+1:MF.t_end
                                  %---------------------------------------------
                                %for thermal 
                                %---------------------------------------------
                                 %construct 2D matrix for nested t1, t2 integral with t1>t2 
                                 %-J_E_t_tau(:,1) is t1-t2, the array need to be fliped when multiply with the f_1(t2)
                                 %-J_E_t_tau(1,:) is t1
                                 % -------------------------------------------------------
                                    temp=obj.J_t(obj.N_t:(obj.N_t+iter_tau-MF.t_start-1));
                                    obj.J_E_t_tau(1:(iter_tau-MF.t_start),iter_tau-MF.t_start+1)=flip(temp);
                                   %---------------------------------------------
                                    %for squeezing
                                    %---------------------------------------------
                                     temp=obj.numerical_squeeze(C,MF,MF.t(iter_tau),MF.t(MF.t_start:(iter_tau-1)));
                                     obj.J_exp_squeeze(1:(iter_tau-MF.t_start),iter_tau-MF.t_start+1)=flip(temp);

                                  end

                                 obj.magic_operation=@(f1_HI,MF) obj.manipulate_short_thermal(f1_HI,MF);
                                obj.magic_operation_squeeze=@(C,f1_HI,MF) obj.manipulate_short_squeeze(C,f1_HI,MF);
            else
                         %for fine dt data points, to save memory, it goes back to
                         %brute force integration
                                  obj.magic_operation=@(f1_HI,MF) obj.manipulate_long_thermal(f1_HI,MF);
                                  obj.magic_operation_squeeze=@(C,f1_HI,MF) obj.manipulate_long_squeeze(C,f1_HI,MF);
            end


end