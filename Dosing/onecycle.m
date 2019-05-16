params = [0.398*24 0.86 175*24 42000*24 2400*24 2.09*24 106*24 7.56*24 81.3*24 14.9*24 2.73*24 0.31*24 6000*24 1.81*24 30300 140 29 1421 2.5 0.0037*24 339.3 0.50 1 0.4 3 (6*10^10)*24 0]';

%% Dosing Regimen 1: Full dose, 5 days ON, 23 days OFF

bolus = 1.3;
                                                                                                                                                                                                                                          % set to 0 for susceptible
% cycle 1
y0 = [bolus 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6]';

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
[T1,Y1] = ode15s(@TMZ_equations,[0 1],y0,options,params);

dim = size(T1);
totalD = repmat(bolus,dim(1),1); % initial drug in

% doses 2-5 

for i = 1:4
    
    y0 = Y1(end,:)';
    y0(1) = y0(1) + bolus;
    
    tspanON = [(i) (i+1)];
    [T1s_new,Y1s_new] = ode15s(@TMZ_equations,tspanON,y0,options,params);
    
    dimON = size(T1s_new);
   
    onD = totalD(end)+bolus; % add bolus
    onD = repmat(onD,dimON(1),1);
    totalD = vertcat(totalD,onD(1:end-1,:));

    T1 = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1 = vertcat(Y1,Y1s_new(1:end-1,:));
   
end

y0 = Y1(end,:)';

tspanOFF = [(5) (5+23)];
[T1s_off,Y1s_off] = ode15s(@TMZ_equations,tspanOFF,y0,options,params);

dimOFF = size(T1s_off);
offD = repmat(totalD(end),dimOFF(1),1);
totalD = vertcat(totalD,offD(1:end-1,:));

T1 = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
Y1 = vertcat(Y1,Y1s_off(1:end-1,:));

%% Dosing Regimen 2: 3/4 dose, 7 days ON, 8 days OFF

bolus2 = 1.3*(3/4);

% dose 1
y0 = [bolus2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6]';
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
[T2,Y2] = ode15s(@TMZ_equations,[0 1],y0,options,params);

dim = size(T2);
totalD2 = repmat(bolus2,dim(1),1); % drug in after first dose

% doses 2-7
y0 = Y2(end,:)'; % initial conditions for dose 2
y0(1) = y0(1) + bolus2; % add in bolus for dose 2

for i = 1:6
    
    tspanON = [(i) (i+1)];
    [T2s_new,Y2s_new] = ode15s(@TMZ_equations,tspanON,y0,options,params);
   
    dimON = size(T2s_new);
   
    onD2 = totalD2(end)+bolus2; % total drug in
    onD2 = repmat(onD2,dimON(1),1);
    totalD2 = vertcat(totalD2,onD2(1:end-1,:)); % add new drug to array
   
    T2 = vertcat(T2,T2s_new(1:end-1,:)); %goes to end-1 to prevent repeats
    Y2 = vertcat(Y2,Y2s_new(1:end-1,:));
    y0 = Y2(end,:)'; % reset initial conditions
    y0(1) = y0(1) + bolus2; % add bolus to initial conditions
   
end

tspanOFF = [7 15];
y0 = Y2(end,:)';
[T2s_off,Y2s_off] = ode15s(@TMZ_equations,tspanOFF,y0,options,params);

dimOFF = size(T2s_off);
offD2 = repmat(totalD2(end),dimOFF(1),1);
totalD2 = vertcat(totalD2,offD2(1:end-1,:));

T2 = vertcat(T2,T2s_off(1:end-1,:)); %goes to end-1 to prevent repeats
Y2 = vertcat(Y2,Y2s_off(1:end-1,:));

y0 = Y2(end,:)'; % reset initial conditions
y0(1) = y0(1) + bolus2; % add bolus to initial conditions


% second cycle

for i = 15:21

    tspanON = [(i) (i+1)];
    [T2s_new,Y2s_new] = ode15s(@TMZ_equations,tspanON,y0,options,params);
   
    dimON = size(T2s_new);
   
    onD2 = totalD2(end)+bolus2; % add bolus
    onD2 = repmat(onD2,dimON(1),1);
   
    T2 = vertcat(T2,T2s_new(1:end-1,:)); %goes to end-1 to prevent repeats
    Y2 = vertcat(Y2,Y2s_new(1:end-1,:));
    totalD2 = vertcat(totalD2,onD2(1:end-1,:));
    y0 = Y2(end,:)';
    y0(1) = y0(1) + bolus2;

end
     
tspanOFF = [22 28];
y0 = Y2(end,:)';
[T2s_off,Y2s_off] = ode15s(@TMZ_equations,tspanOFF,y0,options,params);

dimOFF = size(T2s_off);
offD2 = repmat(totalD2(end),dimOFF(1),1);
totalD2 = vertcat(totalD2,offD2(1:end-1,:));

T2 = vertcat(T2,T2s_off(1:end-1,:)); %goes to end-1 to prevent repeats
Y2 = vertcat(Y2,Y2s_off(1:end-1,:));

%% Dosing Regimen 3: 3/8 dose, 28 days on

bolus3 = 1.3*(3/8);
                                                                                                                                                                                                                                          % set to 0 for susceptible
% cycle 1
y0 = [bolus3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6]';

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
[T3,Y3] = ode15s(@TMZ_equations,[0 1],y0,options,params);

dim = size(T3);

% doses 2-28 
y0 = Y3(end,:)';
y0(1) = y0(1) + bolus3;
totalD3 = repmat(bolus3,dim(1),1); % initial drug in

for i = 1:27
    
    tspanON = [(i) (i+1)];
    [T3s_new,Y3s_new] = ode15s(@TMZ_equations,tspanON,y0,options,params);
    
    dimON = size(T3s_new);
   
    newD3 = totalD3(end)+bolus3; % add bolus
    newD3 = repmat(newD3,dimON(1),1);

    T3 = vertcat(T3,T3s_new(1:end-1,:)); %goes to end-1 to prevent repeats
    Y3 = vertcat(Y3,Y3s_new(1:end-1,:));
    totalD3 = vertcat(totalD3,newD3(1:end-1,:));
    y0 = Y3(end,:)';
    y0(1) = y0(1) + bolus3;
   
end

%% Mass balance 1

TMZm = Y1(:,1)+Y1(:,2)+Y1(:,3)+Y1(:,4)+Y1(:,5)+Y1(:,6)+Y1(:,7)+Y1(:,8);
MTICm = Y1(:,9)+Y1(:,10)+Y1(:,11)+Y1(:,12)+Y1(:,13)+Y1(:,14);
Cm = Y1(:,15)+Y1(:,16)+Y1(:,17);
DNAm = Y1(:,18)+Y1(:,19)+Y1(:,23)+Y1(:,24);
balance = totalD - (TMZm + MTICm + Cm + DNAm);

%% Mass balance 2

TMZm2 = Y2(:,1)+Y2(:,2)+Y2(:,3)+Y2(:,4)+Y2(:,5)+Y2(:,6)+Y2(:,7)+Y2(:,8);
MTICm2 = Y2(:,9)+Y2(:,10)+Y2(:,11)+Y2(:,12)+Y2(:,13)+Y2(:,14);
Cm2 = Y2(:,15)+Y2(:,16)+Y2(:,17);
DNAm2 = Y2(:,18)+Y2(:,19)+Y2(:,23)+Y2(:,24);
balance2 = totalD2 - (TMZm2 + MTICm2 + Cm2 + DNAm2);

%% Mass balance 3

TMZm3 = Y3(:,1)+Y3(:,2)+Y3(:,3)+Y3(:,4)+Y3(:,5)+Y3(:,6)+Y3(:,7)+Y3(:,8);
MTICm3 = Y3(:,9)+Y3(:,10)+Y3(:,11)+Y3(:,12)+Y3(:,13)+Y3(:,14);
Cm3 = Y3(:,15)+Y3(:,16)+Y3(:,17);
DNAm3 = Y3(:,18)+Y3(:,19)+Y3(:,23)+Y3(:,24);
balance3 = totalD3 - (TMZm3 + MTICm3 + Cm3 + DNAm3);

%% Plots

tV1 = Y1(:,20)+Y1(:,21); % tumor volume (changes wrt time)
tV2 = Y2(:,20)+Y2(:,21);
tV3 = Y3(:,20)+Y3(:,21);

% plot plasma TMZ
TMZ_plasma1 = Y1(:,2)/params(15);
TMZ_plasma2 = Y2(:,2)/params(15);
TMZ_plasma3 = Y3(:,2)/params(15);
figure;
plot(T1,TMZ_plasma1,'b',T2,TMZ_plasma2,'r',T3,TMZ_plasma3,'g')
title('Concentration of TMZ in Plasma')
xlabel('Time (days)')
ylabel('Concentration (M)')

% plot tumor TMZ
TMZ_tumor1 = Y1(:,7)./tV1;
TMZ_tumor2 = Y2(:,7)./tV2;
TMZ_tumor3 = Y3(:,7)./tV3;
figure;
plot(T1,TMZ_tumor1,'b',T2,TMZ_tumor2,'r',T3,TMZ_tumor3,'g')
title('Concentration of TMZ in Tumor')
xlabel('Time (days)')
ylabel('Concentration (M)')

% plot tumor MTIC
MTIC_tumor1 = Y1(:,14)./tV1;
MTIC_tumor2 = Y2(:,14)./tV2;
MTIC_tumor3 = Y3(:,14)./tV3;
figure;
plot(T1,MTIC_tumor1,'b',T2,MTIC_tumor2,'r',T3,MTIC_tumor3,'g')
title('Concentration of MTIC in Tumor')
xlabel('Time (days)')
ylabel('Concentration (M)')

% Mass balance
figure;
plot(T1,balance,'b',T2,balance2,'r',T3,balance3,'g')
title('Mass Balance')
xlabel('Time (days)')
ylabel('Balance (mmol)')

%% Save data to .mat file for visualization - fig 2

save singlecycle_reg1.mat T1 TMZ_plasma1 TMZ_tumor1 MTIC_tumor1 balance;
save singlecycle_reg2.mat T2 TMZ_plasma2 TMZ_tumor2 MTIC_tumor2 balance2;
save singlecycle_reg3.mat T3 TMZ_plasma3 TMZ_tumor3 MTIC_tumor3 balance3;