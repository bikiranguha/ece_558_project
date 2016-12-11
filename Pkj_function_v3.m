%% Recursive function to calculate the fraction of original load at each bus 


function frac = Pkj_function_v3(PDk, post_flow1, post_flow2, CapAva1, CapAva2, MWmax, frac_init, tr_out_no, LODF_all, PTDF_all)

% Inputs:
% PDk: vector of loads at each bus
% post_flow1: % Branch flows (including transmission outages) assuming all load is being supplied by Gen Bus 1
% post_flow2: % Branch flows (including transmission outages) assuming all load is being supplied by Gen Bus 2
% CapAva1: Capacity available at Gen Bus 1
% CapAva2: Capacity available at Gen Bus 2
% MWmax: vector of branch flow limits (MW)
% frac_init: Initial value of frac, generally 1
% tr_out_no: Transmission line number which is out
% LODF_all: LODF matrix
% PTDF_all: PTDF_all

% Output:
% frac: Ratio of PDk/PD_pk

Ltotal = sum(PDk);

success = gen_dispatch(Ltotal, post_flow1, post_flow2, CapAva1, CapAva2, MWmax);



if (success == 1)
        frac = frac_init;
     disp(' Excellent! Generation  and Transmission constraints have been met!')
else
    frac = frac_init - 0.1;
    PDk = frac*PDk; % Curtailing the total load by 0.1*PD_pk
    Ltotal= sum(PDk);
    
    % Getting the new Pinj vectors
    Pinj1 = -PDk;
    Pinj2 = - PDk;
    Pinj1(1,1) = Ltotal;
    Pinj2(2,1) = Pinj2(2,1)+ Ltotal;
   
    % Calculating the new flows given that any 1 Gen Bus is supplying
    % entire load
    pre_flow1 = PTDF_all*Pinj1;
    pre_flow2 = PTDF_all*Pinj2;
    
    % Considering any transmission line outage and calculating post flows
    if (tr_out_no == 0)
       post_flow1 = pre_flow1;
   else
       post_flow1 = pre_flow1 + LODF_all(:,tr_out_no)*pre_flow1(tr_out_no);
   end
   
   if (tr_out_no == 0)
       post_flow2 = pre_flow2;
   else
       post_flow2 = pre_flow2 + LODF_all(:,tr_out_no)*pre_flow2(tr_out_no);
   end

    
    if(frac > 0)
    frac = Pkj_function_v3(PDk, post_flow1, post_flow2, CapAva1, CapAva2, MWmax, frac, tr_out_no, LODF_all, PTDF_all);
    else
        frac = 0;
        disp(' Dispatch totally failed! No load can be served!')
    end

end

end