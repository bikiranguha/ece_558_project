% Appendix
% Final Project Bikiran Guha (A20378383)


%% Frequency and Duration Based Generator Capacity Availability Tables


% Bus 1

lam_b1 = 1.1; % lambda (failures/year) for each generator in Bus 1 
mu_b1 = 73;   % mu (repairs/year) for each generator in Bus 1 

u_b1 = (lam_b1)/(lam_b1 + mu_b1); % FOR for generators at Bus 1
c_b1 = 20; % Capacity of each generator unit in Bus 1

% CO_b1: Capa. out in Bus 1, CI_b1: Capa. in. Bus 1
CO_b1 = c_b1*[0:4];  % Since there are 4*20 MW generators
CI_b1 = 80 - CO_b1;  % Capacity in service values

% Getting individual prob values for generator outages in Bus 1
for(i=1:length(CO_b1))
    Pr_b1(i,1) = prob((i-1)*c_b1,u_b1,c_b1,4);
end


% Eliminate values of Pr_b1, CO_b1 and CI_b1 where Pr_b1 < 1e-6
todelete = []; % vector to save all the rows where Pr_b1 < 1e-6

for (i = 1:length(Pr_b1))
    if(Pr_b1(i,1) < 1e-6)
        todelete = [todelete i];
    end
end

Pr_b1(todelete) = [];
CO_b1(todelete) = [];
CI_b1(todelete) = [];


% Getting positive departure rates
for (i=1:length(CO_b1))
% List of  positive (repair) departure rate values (per year)
lmp_b1(i,1) = lambdap(CO_b1(i),u_b1,c_b1,4, mu_b1); 
        
end

% Getting negative departure rates
for (i=1:length(CO_b1))
 % List of  negative (failure) departure rate values (per year)
lmn_b1(i,1) = lambdan(CO_b1(i),u_b1,c_b1,4,lam_b1); 
        
end

% Frequency values
for (i=1:length(CO_b1))
    freq_b1(i,1) = Pr_b1(i,1)*(lmp_b1(i,1)+lmn_b1(i,1));
end
 
% Creating the frequency and duration table for Bus 1 Generators
tbl_b1(:,1) = [1:length(CO_b1)]';
tbl_b1(:,2) = [0:length(CO_b1) - 1]';
tbl_b1(:,3) = CI_b1;
tbl_b1(:,4) = Pr_b1;
tbl_b1(:,5) = lmp_b1;
tbl_b1(:,6) = lmn_b1;
tbl_b1(:,7) = freq_b1;


disp('Frequency and Duration Table for Bus 1 Generators')
disp('  State       Units   Capa In   Prob.    Positive    Negative   Frequency')
disp('              out                       dep. rate   dep. rate    (occ/yr)')
disp(tbl_b1)




% Bus 2

lam_b2 = 0.5; % lambda (failures/year) for each generator in Bus 2 
mu_b2 = 100;   % mu (repairs/year) for each generator in Bus 2


u_b2 = (lam_b2)/(lam_b2 + mu_b2); % FOR for generators at Bus 2


 % CO_b2: Capacity in service values
CO_b2 = 0:5:130;  % Bus 2: 7*5 MW, 15 MW, 4*20 MW, Total: 130 MW
CI_b2 = 130 - CO_b2;  % Capacity in service values

% List of FOR values for each generator in Bus 2
u_b2_list = u_b2*ones(1,12); 
c_b2_5 = 5*ones(1,7);
c_b2_20 = 20*ones(1,4);

% List of capa. values for each generator in Bus 2
c_b2_list = [c_b2_5 15 c_b2_20];  

% Getting individual prob values for generator outages in Bus 2
for(i=1:length(CO_b2))
    Pr_b2(i,1) = prob_v2(CO_b2(i),u_b2_list,c_b2_list);
end

% Eliminate values of Pr_b2, CO_b2 and CI_b2 where Pr_b2 < 1e-6
todelete = []; % vector to save all the rows where Pr_b2 < 1e-6

for (i = 1:length(Pr_b2))
    if(Pr_b2(i,1) < 1e-6)
        todelete = [todelete i];
    end
end

Pr_b2(todelete) = [];
CO_b2(todelete) = [];
CI_b2(todelete) = [];


% Getting positive departure rates
for (i=1:length(CO_b2))
    % List of  positive (repair) departure rate values (per year)
lmp_b2(i,1) = lambdap_v2(CO_b2(i),u_b2_list,c_b2_list, mu_b2); 
        
end

% Getting negative departure rates
for (i=1:length(CO_b2))
 % List of  negative (failure) departure rate values (per year)       
lmn_b2(i,1) = lambdan_v2(CO_b2(i),u_b2_list,c_b2_list,lam_b2); 
end


% Frequency
for (i=1:length(CO_b2))
    freq_b2(i,1) = Pr_b2(i,1)*(lmp_b2(i,1)+lmn_b2(i,1));
end



% Creating the frequency and duration table for Bus 2 Generators
tbl_b2(:,1) = [1:length(CI_b2)]';
tbl_b2(:,2) = CI_b2';
tbl_b2(:,3) = Pr_b2;
tbl_b2(:,4) = lmp_b2;
tbl_b2(:,5) = lmn_b2;
tbl_b2(:,6) = freq_b2;

disp('Frequency and Duration Table for Bus 2 Generators')
disp('  State       Capa In   Prob.    Positive    Negative   Frequency')
disp('                                 dep. rate   dep. rate    (occ/yr)')
disp(tbl_b2)




% Generator subsystem


% Getting a probability array for the generator subsystem
for (i = 1: length(Pr_b1))
    for (j=1:length(Pr_b2))
        
        Pr_b1_b2(i,j) = Pr_b1(i)*Pr_b2(j); % The probability array
    end
end


% The following 7 variables will contain all the info in the generator
% subsystem table
Pr_b1_b2_trunc = [];  % Individual prob.
CI_b1_trunc = [];     % Capa. at Bus 1
CI_b2_trunc = [];     % Capa. at Bus 2
lmp_b1_b2_trunc = []; % Positive departure rate 
lmn_b1_b2_trunc = []; % Negative departure rate
freq_b1_b2_trunc = []; % Frequency


% Getting all the values for the generator subsystem table
for (i = 1: length(Pr_b1))
    for (j=1:length(Pr_b2))
        
        if (Pr_b1_b2(i,j) >= 5e-7)
        
        CI_b1_trunc = [CI_b1_trunc; CI_b1(i)];
        CI_b2_trunc = [CI_b2_trunc; CI_b2(j)];
        Pr_b1_b2_trunc = [Pr_b1_b2_trunc; Pr_b1_b2(i,j)];
        lmp_b1_b2 = lmp_b1(i)+lmp_b2(j);
        lmn_b1_b2 = lmn_b1(i)+lmn_b2(j);
        lmp_b1_b2_trunc = [lmp_b1_b2_trunc; lmp_b1_b2];
        lmn_b1_b2_trunc = [lmn_b1_b2_trunc; lmn_b1_b2];
        freq_b1_b2_trunc = [freq_b1_b2_trunc; Pr_b1_b2(i,j)*(lmp_b1_b2 + lmn_b1_b2)];
        
        end
        
    end
end



% Creating the frequency and duration table for the generator subsystem
tbl_b1_b2(:,1) = [1:length(Pr_b1_b2_trunc)]';
tbl_b1_b2(:,2) = CI_b1_trunc;
tbl_b1_b2(:,3) = CI_b2_trunc;
tbl_b1_b2(:,4) = Pr_b1_b2_trunc;
tbl_b1_b2(:,5) = lmp_b1_b2_trunc;
tbl_b1_b2(:,6) = lmn_b1_b2_trunc;
tbl_b1_b2(:,7) = freq_b1_b2_trunc;


disp('Frequency and Duration Table for Generator Subsystem')
disp('  State    Gen. Capa.   Gen. Capa. Prob.  Positive    Negative   Frequency')
disp('           Bus 1        Bus 2             dep. rate   dep. rate  (occ/yr)')
disp(tbl_b1_b2)





%% Frequency and Duration Based Transmission Capacity Availability Table

% Failure rate for each branch (occ/yr)
lam_br1 = 1.5;
lam_br2 = 5.0;
lam_br3 = 4.0;
lam_br4 = 1.0;
lam_br5 = 1.0;
lam_br6 = 1.5;
lam_br7 = 5.0;
lam_br8 = 1.0;

% Repair rates for the branches (occ/yr)

mu_br = 1/(10/8760);


% FOR for each branch 
u_br1 = lam_br1/(lam_br1+mu_br);
u_br2 = lam_br2/(lam_br2+mu_br);
u_br3 = lam_br3/(lam_br3+mu_br);
u_br4 = lam_br4/(lam_br4+mu_br);
u_br5 = lam_br5/(lam_br5+mu_br);
u_br6 = lam_br6/(lam_br6+mu_br);
u_br7 = lam_br7/(lam_br7+mu_br);
u_br8 = lam_br8/(lam_br8+mu_br);


% Availability for each branch
a_br1 = 1-u_br1;
a_br2 = 1-u_br2;
a_br3 = 1-u_br3;
a_br4 = 1-u_br4;
a_br5 = 1-u_br5;
a_br6 = 1-u_br6;
a_br7 = 1-u_br7;
a_br8 = 1-u_br8;



% Individual probability for each outage state

Pr_br0 = a_br1*a_br2*a_br3*a_br4*a_br5*a_br6*a_br7*a_br8; % All branches in service
Pr_br1 = u_br1*a_br2*a_br3*a_br4*a_br5*a_br6*a_br7*a_br8; % Branch 1 out
Pr_br2 = a_br1*u_br2*a_br3*a_br4*a_br5*a_br6*a_br7*a_br8; % Branch 2 out and so on...
Pr_br3 = a_br1*a_br2*u_br3*a_br4*a_br5*a_br6*a_br7*a_br8;
Pr_br4 = a_br1*a_br2*a_br3*u_br4*a_br5*a_br6*a_br7*a_br8;
Pr_br5 = a_br1*a_br2*a_br3*a_br4*u_br5*a_br6*a_br7*a_br8;
Pr_br6 = a_br1*a_br2*a_br3*a_br4*a_br5*u_br6*a_br7*a_br8;
Pr_br7 = a_br1*a_br2*a_br3*a_br4*a_br5*a_br6*u_br7*a_br8;
Pr_br8 = a_br1*a_br2*a_br3*a_br4*a_br5*a_br6*a_br7*u_br8;

Pr_br = [ Pr_br0; Pr_br1; Pr_br2; Pr_br3; Pr_br4; Pr_br5; Pr_br6; Pr_br7 ; Pr_br8];


% Positive departure rate for each outage state

lmp_br0 = 0;     % Repair rate all branches are in service
lmp_br1 = mu_br; % Repair rate when branch 1 is out
lmp_br2 = mu_br; % Branch 2 out and so on....
lmp_br3 = mu_br;
lmp_br4 = mu_br;
lmp_br5 = mu_br;
lmp_br6 = mu_br;
lmp_br7 = mu_br;
lmp_br8 = mu_br;

lmp_br = [ lmp_br0; lmp_br1; lmp_br2; lmp_br3; lmp_br4; lmp_br5; lmp_br6; lmp_br7 ; lmp_br8];



% Negative departure rate for each outage state

 % Negative departure rate when all branches in service
lmn_br0 = (lam_br1 + lam_br2 + lam_br3 + lam_br4 + lam_br5 + lam_br6 + lam_br7 + lam_br8);

 % Since a max of only one branch can go out at any point of time
lmn_br = [lmn_br0; zeros(8,1)];



% Frequency for each state

freq_br = Pr_br.*(lmp_br+lmn_br);



% Creating the frequency and duration based transmission availability table
tbl_br(:,1) = [1:9]';
tbl_br(:,2) = [0:8]';
tbl_br(:,3) = Pr_br;
tbl_br(:,4) = lmp_br;
tbl_br(:,5) = lmn_br;
tbl_br(:,6) = freq_br;



disp('Frequency and Duration Based Transmission Capacity Availability Table')
disp('  State       Branch    Prob.    Positive    Negative   Frequency')
disp('   No.        Name (L#)         dep. rate   dep. rate   (occ/yr)')
disp(tbl_br)


%% Frequency and Duration Based Composite Generation and Transmission Availability Table


% The following 9 variables will contain all the info in the composite
% generator and transmission availability table

Pr_b12_br = [];     % Individual Probability
tr_out = [];        % Transmission line out name
g1ca = [];          % Capacity available in Bus 1 Generators
g2ca = [];          % Capacity available in Bus 1 Generators
lmp_b12_br = [];    % Positive departure rate
lmn_b12_br = [];    % Negative departure rate
freq_b12_br = [];   % Frequency
dur_b12_br = [];    % Duration


for (i = 1:length(Pr_br))
    for (j=1:length(Pr_b1_b2_trunc))
        
        tr_out = [tr_out; tbl_br(i,2)];
        g1ca = [g1ca; CI_b1_trunc(j)];
        g2ca = [g2ca; CI_b2_trunc(j)];
        
        Pr_temp = Pr_b1_b2_trunc(j)*Pr_br(i);
        Pr_b12_br = [Pr_b12_br; Pr_temp];
        
        lmp_temp = lmp_br(i)+lmp_b1_b2_trunc(j);
        lmp_b12_br = [lmp_b12_br; lmp_temp];
        
        lmn_temp = lmn_br(i)+lmn_b1_b2_trunc(j);
        lmn_b12_br = [lmn_b12_br; lmn_temp];
        
        freq_temp = Pr_temp*(lmp_temp + lmn_temp);
        freq_b12_br = [freq_b12_br; freq_temp];
        
        dur_temp = Pr_temp/freq_temp*8760; % in hours
        dur_b12_br = [dur_b12_br; dur_temp];
        
        
    end
end

        
% Creating the frequency and duration based composite generation and 
% transmission availability table
tbl_b12_br(:,1) = [1:length(tr_out)]';
tbl_b12_br(:,2) = tr_out;
tbl_b12_br(:,3) = g1ca;
tbl_b12_br(:,4) = g2ca;
tbl_b12_br(:,5) = Pr_b12_br;
tbl_b12_br(:,6) = lmp_b12_br;
tbl_b12_br(:,7) = lmn_b12_br;
tbl_b12_br(:,8) = freq_b12_br;
tbl_b12_br(:,9) = dur_b12_br;


disp('Frequency and Duration Based Composite Generation and Transmission Availability Table')
disp('  State          TrOut    G1CA    G2CA      Prob.    PosDepR    NegDepR   Freq.     Dur.')
disp('                  (L#)                               (r/yr)     (f/yr)  (occ/yr)    (hr)')
disp(tbl_b12_br)



%% Generating load point indices

load('PTDF_LODF.mat')

% Peak load at each bus. Rows correspond to buses
PD_pk = [0; 20; 85; 40; 10]; 
Lkj= []; % Matrix of load curtailed (MW) for each load bus for every state
CapaL = []; % Matrix of capacity available at each bus
ELC = []; % Expected Load Curtailed matrix
EENS = []; % EENS values for all the load buses and composite states


% Loop to generate reliability indices for each state
for (i=1:length(tbl_b12_br))
    frac_init = 1;
   MWmax = 71*ones(8,1); 
   Ltotal = sum(PD_pk);
   Pinj1 = [Ltotal; -20; -85; -40; -10]; % Assuming all load is being supplied by Gen 1
   Pinj2 = [0; Ltotal-20;  -85; -40; -10]; % Assuming all load is being supplied by Gen 2

    pre_flow1 = PTDF_all*Pinj1;
    pre_flow2 = PTDF_all*Pinj2;

   
   CapAva1 = g1ca(i);
   CapAva2 = g2ca(i);
   tr_out_no = tr_out(i);
   
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
   
   frac(i) = Pkj_function_v3(PD_pk, post_flow1, post_flow2, CapAva1, CapAva2, MWmax, frac_init, tr_out_no, LODF_all, PTDF_all);
   
   if (frac(i) == 1)
       Pkj(i) = 0; % Note that Pkj will be same for all load buses
   else
       Pkj(i) = 1;
   end
   
   PDk = frac(i)*PD_pk;
   
   % Rows of CapaL represent different buses, while the columns represent different composite states 
   CapaL = [CapaL PDk]; 
   % Rows of Lkj represent different buses, while the columns represent different composite states 
   Lkj = [Lkj (PD_pk-PDk)]; 
   % Rows of ELC represent different buses, while the columns represent different composite states 
   ELC = [ELC Lkj(:,i)*freq_b12_br(i,1)]; 
   
   if (Pkj(i) == 1)
       NLC(i,1) = freq_b12_br(i,1);  % NLC: No. of load curtailments
   else
       NLC(i,1) = 0;  % Note that NLC will be same for all load buses
   end
   
   % Rows of EENS represent different buses, while the columns represent different composite states 
   EENS = [EENS Lkj(:,i)*Pr_b12_br(i,1)*8760]; % Units: MWh
   
   if (Pkj(i) == 1)
       %EDLC: Expected duration of load curtailments
       EDLC(i,1) = Pr_b12_br(i,1)*8760;  
   else
       % Note that EDLC will be same for all load buses
       EDLC(i,1) = 0;                 
   end
   


end


% Building the table for Bus 2

tbl_lb2(:,1) = tbl_b12_br(:,1); % State number
tbl_lb2(:,2) = CapaL(2,:); % Capacity available or load served at the bus
tbl_lb2(:,3) = Pkj;        % Probability of inadequacy
tbl_lb2(:,4) = Lkj(2,:);   % Load curtailed
tbl_lb2(:,5) = ELC(2,:);   % Expected load curtailed
tbl_lb2(:,6) = NLC;        % No. of load curtailments
tbl_lb2(:,7) = EENS(2,:);  % Expected Energy Not Served
tbl_lb2(:,8) = EDLC;       % Expected duration of load curtailments


% Eliminating all rows where Pkj = 0 from the table

tbl_lb2_comp = tbl_lb2; 

todelete = []; % vector to save all the rows where Pkj = 0

for (i = 1:length(tbl_lb2(:,1)))
    if(Pkj(i) == 0)
        todelete = [todelete i];
    end
end

tbl_lb2_comp(todelete,:) = [];  % New more compact table

disp('Load Point Reliability Indices Table for Bus 2')
disp('  State      CapaAvaL     Pkj      Lkj        ELC       NLC      EENS      EDLC')
disp('               (MW)                (MW)       (MW)     (occ)    (MWh)     (hr)')
disp(tbl_lb2_comp)





% Building the table for Bus 3

tbl_lb3(:,1) = tbl_b12_br(:,1); % State number
tbl_lb3(:,2) = CapaL(3,:); % Capacity available or load served at the bus
tbl_lb3(:,3) = Pkj;        % Probability of inadequacy
tbl_lb3(:,4) = Lkj(3,:);   % Load curtailed
tbl_lb3(:,5) = ELC(3,:);   % Expected load curtailed
tbl_lb3(:,6) = NLC;        % No. of load curtailments
tbl_lb3(:,7) = EENS(3,:);  % Expected Energy Not Served
tbl_lb3(:,8) = EDLC;       % Expected duration of load curtailments


% Eliminating all rows where Pkj = 0 from the table

tbl_lb3_comp = tbl_lb3; 

todelete = []; % vector to save all the rows where Pkj = 0

for (i = 1:length(tbl_lb3(:,1)))
    if(Pkj(i) == 0)
        todelete = [todelete i];
    end
end

tbl_lb3_comp(todelete,:) = [];  % New more compact table

disp('Load Point Reliability Indices Table for Bus 3')
disp('  State      CapaAvaL     Pkj      Lkj        ELC       NLC      EENS      EDLC')
disp('               (MW)                (MW)       (MW)     (occ)    (MWh)     (hr)')
disp(tbl_lb3_comp)




% Building the table for Bus 4

tbl_lb4(:,1) = tbl_b12_br(:,1); % State number
tbl_lb4(:,2) = CapaL(4,:); % Capacity available or load served at the bus
tbl_lb4(:,3) = Pkj;        % Probability of inadequacy
tbl_lb4(:,4) = Lkj(4,:);   % Load curtailed
tbl_lb4(:,5) = ELC(4,:);   % Expected load curtailed
tbl_lb4(:,6) = NLC;        % No. of load curtailments
tbl_lb4(:,7) = EENS(4,:);  % Expected Energy Not Served
tbl_lb4(:,8) = EDLC;       % Expected duration of load curtailments


% Eliminating all rows where Pkj = 0 from the table

tbl_lb4_comp = tbl_lb4; 

todelete = []; % vector to save all the rows where Pkj = 0

for (i = 1:length(tbl_lb4(:,1)))
    if(Pkj(i) == 0)
        todelete = [todelete i];
    end
end

tbl_lb4_comp(todelete,:) = [];  % New more compact table


disp('Load Point Reliability Indices Table for Bus 4')
disp('  State      CapaAvaL     Pkj      Lkj        ELC       NLC      EENS      EDLC')
disp('               (MW)                (MW)       (MW)     (occ)    (MWh)     (hr)')
disp(tbl_lb4_comp)



% Building the table for Bus 5

tbl_lb5(:,1) = tbl_b12_br(:,1); % State number
tbl_lb5(:,2) = CapaL(5,:); % Capacity available or load served at the bus
tbl_lb5(:,3) = Pkj;        % Probability of inadequacy
tbl_lb5(:,4) = Lkj(5,:);   % Load curtailed
tbl_lb5(:,5) = ELC(5,:);   % Expected load curtailed
tbl_lb5(:,6) = NLC;        % No. of load curtailments
tbl_lb5(:,7) = EENS(5,:);  % Expected Energy Not Served
tbl_lb5(:,8) = EDLC;       % Expected duration of load curtailments


% Eliminating all rows where Pkj = 0 from the table

tbl_lb5_comp = tbl_lb5; 

todelete = []; % vector to save all the rows where Pkj = 0

for (i = 1:length(tbl_lb5(:,1)))
    if(Pkj(i) == 0)
        todelete = [todelete i];
    end
end

tbl_lb5_comp(todelete,:) = [];  % New more compact table
    

disp('Load Point Reliability Indices Table for Bus 5')
disp('  State      CapaAvaL     Pkj      Lkj        ELC       NLC      EENS      EDLC')
disp('               (MW)                (MW)       (MW)     (occ)    (MWh)     (hr)')
disp(tbl_lb5_comp)


% Summary of results

disp('Summary of load point reliability indices:')

ELC_b2 = sum(tbl_lb2_comp(:,5));
X = sprintf('Expected Load Curtailed at Bus 2 = %d MW',ELC_b2);
disp(X)
ELC_b3 = sum(tbl_lb3_comp(:,5));
X = sprintf('Expected Load Curtailed at Bus 3 = %d MW',ELC_b3);
disp(X)
ELC_b4 = sum(tbl_lb4_comp(:,5));
X = sprintf('Expected Load Curtailed at Bus 4 = %d MW',ELC_b4);
disp(X)
ELC_b5 = sum(tbl_lb5_comp(:,5));
X = sprintf('Expected Load Curtailed at Bus 5 = %d MW',ELC_b5);
disp(X)

NLC_sum_b1 = sum(tbl_lb2_comp(:,6));
X = sprintf('NLC for any load bus = %d',NLC_sum_b1);
disp(X)


EENS_b2 = sum(tbl_lb2_comp(:,7));
X = sprintf('EENS at Bus 2 = %d MWh',EENS_b2);
disp(X)
EENS_b3 = sum(tbl_lb3_comp(:,7));
X = sprintf('EENS at Bus 3 = %d MWh',EENS_b3);
disp(X)
EENS_b4 = sum(tbl_lb4_comp(:,7));
X = sprintf('EENS at Bus 4 = %d MWh',EENS_b4);
disp(X)
EENS_b5 = sum(tbl_lb5_comp(:,7));
X = sprintf('EENS at Bus 5 = %d MWh',EENS_b5);
disp(X)

EDLC_sum_b1 = sum(tbl_lb2_comp(:,8));
X = sprintf('EDLC for any load bus = %d hours',EDLC_sum_b1);
disp(X)



%% Testing recommendations

% Using a loop to make a table (fstate) which lists the values of TrOut, 
%G1CA, G2CA which cause a failure
for (i = 1:length(tbl_lb2_comp))

fstate(i,1) = tbl_lb2_comp(i);
fstate(i,2) = tr_out(tbl_lb2_comp(i));
fstate(i,3) = g1ca(tbl_lb2_comp(i));
fstate(i,4) = g2ca(tbl_lb2_comp(i));

end


disp(' Table showing details of composite states which cause load curtailment')
disp('State   TrOut  G1CA  G2CA')
disp(fstate)


% Testing Recommendation about increasing  branch 1 capacity

MWmax = 71*ones(8,1); 


% Outer loop gets values for different MWmax for branch 1
for (i = 1:10)
 Lkj= []; 
CapaL = []; 
ELC = []; 
EENS = [];
    
    
 % Inner loop gets the reliability indices for the MWmax value
 for (j=1:length(tbl_b12_br))
    frac_init = 1; 
   Ltotal = sum(PD_pk);
   Pinj1 = [Ltotal; -20; -85; -40; -10]; % Assuming all load is being supplied by Gen 1
   Pinj2 = [0; Ltotal-20;  -85; -40; -10]; % Assuming all load is being supplied by Gen 2

    pre_flow1 = PTDF_all*Pinj1;
    pre_flow2 = PTDF_all*Pinj2;

   
   CapAva1 = g1ca(j);
   CapAva2 = g2ca(j);
   tr_out_no = tr_out(j);
   
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
   
   frac(j) = Pkj_function_v3(PD_pk, post_flow1, post_flow2, CapAva1, CapAva2, MWmax, frac_init, tr_out_no, LODF_all, PTDF_all);
   
   if (frac(j) == 1)
       Pkj(j) = 0; % Note that Pkj will be same for all load buses
   else
       Pkj(j) = 1;
   end
   
   PDk = frac(j)*PD_pk;
   
   CapaL = [CapaL PDk]; % Rows of CapaL represent different buses, while the columns represent different composite states 
   Lkj = [Lkj (PD_pk-PDk)]; % Rows of lkj represent different buses, while the columns represent different composite states 
   ELC = [ELC Lkj(:,j)*freq_b12_br(j,1)]; % Rows of ELC represent different buses, while the columns represent different composite states 
   
   if (Pkj(j) == 1)
       NLC(j,1) = freq_b12_br(j,1);  % NLC: No. of load curtailments
   else
       NLC(j,1) = 0;                 % Note that NLC will be same for all load buses
   end
   
   % Rows of EENS represent different buses, while the columns represent different composite states 
   EENS = [EENS Lkj(:,j)*Pr_b12_br(j,1)*8760]; % Units: MWh
   
   if (Pkj(j) == 1)
       EDLC(j,1) = Pr_b12_br(j,1)*8760;  %EDLC: Expected duration of load curtailments
   else
       EDLC(j,1) = 0;                 % Note that EDLC will be same for all load buses
   end
   


 end

 MWmax_b1(i) = MWmax(1,1);
 MWmax(1,1) = MWmax(1,1) + 1; % Increasing the MWmax of branch 1 by 1 
 
 
 
 % Calculating total ELC (for all buses)
 ELC_row_sum = sum(ELC);
 ELC_sum_b1(i) = sum(ELC_row_sum); 
 
  % Calculating total NLC 
 NLC_sum_b1(i) = sum(NLC);
 
 % Calculating total ELC (for all buses)
 EENS_row_sum = sum(EENS);
 EENS_sum_b1(i) = sum(EENS_row_sum); 
 
   % Calculating total EDLC 
 EDLC_sum_b1(i) = sum(EDLC);
 
end


reco_table_b1(:,1) = MWmax_b1;
reco_table_b1(:,2) =  ELC_sum_b1;
reco_table_b1(:,3) = NLC_sum_b1;
reco_table_b1(:,4) = EENS_sum_b1;
reco_table_b1(:,5) = EDLC_sum_b1;


disp('Table showing effect of Branch 1 MWmax increase on load reliability indices')
disp('    MWmax      ELC       NLC       EENS      EDLC')
disp('               (MW)                (MWh)      (hr)')
disp(reco_table_b1)









 % Recommendation about increasing  branch 6 capacity
 
 MWmax = 71*ones(8,1); 
 
for (i = 1:10)
 Lkj= []; % This matrix will contain the values of load curtailed (MW) for each load bus for every state
CapaL = []; % This matrix will contain the values of capacity available at each bus
ELC = []; % Expected load curtailed Matrix
EENS = []; % This will contain the EENS values for all the load buses and composite states
    
    
    
 for (j=1:length(tbl_b12_br))
    frac_init = 1; 
   Ltotal = sum(PD_pk);
   Pinj1 = [Ltotal; -20; -85; -40; -10]; % Assuming all load is being supplied by Gen 1
   Pinj2 = [0; Ltotal-20;  -85; -40; -10]; % Assuming all load is being supplied by Gen 2

    pre_flow1 = PTDF_all*Pinj1;
    pre_flow2 = PTDF_all*Pinj2;

   
   CapAva1 = g1ca(j);
   CapAva2 = g2ca(j);
   tr_out_no = tr_out(j);
   
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
   
   frac(j) = Pkj_function_v3(PD_pk, post_flow1, post_flow2, CapAva1, CapAva2, MWmax, frac_init, tr_out_no, LODF_all, PTDF_all);
   
   if (frac(j) == 1)
       Pkj(j) = 0; % Note that Pkj will be same for all load buses
   else
       Pkj(j) = 1;
   end
   
   PDk = frac(j)*PD_pk;
   
   CapaL = [CapaL PDk]; % Rows of CapaL represent different buses, while the columns represent different composite states 
   Lkj = [Lkj (PD_pk-PDk)]; % Rows of lkj represent different buses, while the columns represent different composite states 
   ELC = [ELC Lkj(:,j)*freq_b12_br(j,1)]; % Rows of ELC represent different buses, while the columns represent different composite states 
   
   if (Pkj(j) == 1)
       NLC(j,1) = freq_b12_br(j,1);  % NLC: No. of load curtailments
   else
       NLC(j,1) = 0;                 % Note that NLC will be same for all load buses
   end
   
   % Rows of EENS represent different buses, while the columns represent different composite states 
   EENS = [EENS Lkj(:,j)*Pr_b12_br(j,1)*8760]; % Units: MWh
   
   if (Pkj(j) == 1)
       EDLC(j,1) = Pr_b12_br(j,1)*8760;  %EDLC: Expected duration of load curtailments
   else
       EDLC(j,1) = 0;                 % Note that EDLC will be same for all load buses
   end
   


 end

 MWmax_b6(i) = MWmax(6,1);
 MWmax(6,1) = MWmax(6,1) + 1; % Increasing the MWmax of branch 6 by 1 
 
 
 
 % Calculating total ELC (for all buses)
 ELC_row_sum = sum(ELC);
 ELC_sum_b6(i) = sum(ELC_row_sum); 
 
  % Calculating total NLC 
 NLC_sum_b6(i) = sum(NLC);
 
 % Calculating total ELC (for all buses)
 EENS_row_sum = sum(EENS);
 EENS_sum_b6(i) = sum(EENS_row_sum); 
 
   % Calculating total EDLC 
 EDLC_sum_b6(i) = sum(EDLC);
 
end


reco_table_b6(:,1) = MWmax_b6;
reco_table_b6(:,2) = ELC_sum_b6;
reco_table_b6(:,3) = NLC_sum_b6;
reco_table_b6(:,4) = EENS_sum_b6;
reco_table_b6(:,5) = EDLC_sum_b6;

disp('Table showing effect of Branch 6 MWmax increase on load reliability indices')
disp('    MWmax      ELC       NLC       EENS      EDLC')
disp('               (MW)                (MWh)      (hr)')
disp(reco_table_b6)

disp(' ') 
disp(' ') 
disp('End of project! Thanks for your time and attention! Happy Holidays!')    
    
    
    



