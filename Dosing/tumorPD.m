%% Dosing Regimen 1: Full dose, 5 days ON, 23 days OFF

bolus = 1.3;
params = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.0037*24*28 339.3 0.50 1 0.4 3/28 (6*10^10)*28*24 0]';

% cycle 1
y0 = [bolus 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6]';
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
[T1,Y1] = ode15s(@TMZ_equations,[0 1/28],y0,options,params);

for i = 1:4
    
    y0 = Y1(end,:)';
    y0(1) = y0(1) + bolus;
    
    tspanON = [(i*1/28) (1/28+i*1/28)];
    [T1s_new,Y1s_new] = ode15s(@TMZ_equations,tspanON,y0,options,params);
   
    T1 = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1 = vertcat(Y1,Y1s_new(1:end-1,:));
   
end

y0 = Y1(end,:)';

tspanOFF = [((1/28)*5) ((1/28)*5+23*(1/28))];
[T1s_off,Y1s_off] = ode15s(@TMZ_equations,tspanOFF,y0,options,params);

T1 = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
Y1 = vertcat(Y1,Y1s_off(1:end-1,:));

% additional cycles
for cycle = 1:5 

    for i = 1:5
        
        y0 = Y1(end,:)';
        y0(1) = y0(1) + bolus;
        
        tspan_cyclesON = [(cycle+(i-1)*1/28) (cycle+(i)*(1/28))];

        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,tspan_cyclesON,y0,options,params);

        T1 = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1 = vertcat(Y1,Y1s_new(1:end-1,:));

    end
    
    y0 = Y1(end,:)';
    tspan_cyclesOFF = [(cycle+(1/28)*5) (cycle+1)];

    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,tspan_cyclesOFF,y0,options,params);

    T1 = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1 = vertcat(Y1,Y1s_off(1:end-1,:));

end

%% Dosing Regimen 2: 3/4 dose, 7 days ON, 8 days OFF

params = [0.398*24 0.86 175*24 42000*24 2400*24 2.09*24 106*24 7.56*24 81.3*24 14.9*24 2.73*24 0.31*24 6000*24 1.81*24 30300 140 29 1421 2.5 0.0037*24 339.3 0.50 1 0.4 3 (6*10^10)*24 0]';
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

tspanOFF = [8 15];
y0 = Y2(end,:)';
[T2s_off,Y2s_off] = ode15s(@TMZ_equations,tspanOFF,y0,options,params);

dimOFF = size(T2s_off);
offD2 = repmat(totalD2(end),dimOFF(1),1);
totalD2 = vertcat(totalD2,offD2(1:end-1,:));

T2 = vertcat(T2,T2s_off(1:end-1,:)); %goes to end-1 to prevent repeats
Y2 = vertcat(Y2,Y2s_off(1:end-1,:));

y0 = Y2(end,:)'; % reset initial conditions
y0(1) = y0(1) + bolus2; % add bolus to initial conditions

% additional cycles
for j = 1:10
    for i = 1:7

        tspanON = [((j*15)+i) ((j*15)+(i+1))];
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

tspanOFF = [((15*j)+8) ((15*j)+15)];
y0 = Y2(end,:)';
[T2s_off,Y2s_off] = ode15s(@TMZ_equations,tspanOFF,y0,options,params);

dimOFF = size(T2s_off);
offD2 = repmat(totalD2(end),dimOFF(1),1);
totalD2 = vertcat(totalD2,offD2(1:end-1,:));

T2 = vertcat(T2,T2s_off(1:end-1,:)); %goes to end-1 to prevent repeats
Y2 = vertcat(Y2,Y2s_off(1:end-1,:));

y0 = Y2(end,:)'; % reset initial conditions
y0(1) = y0(1) + bolus2; % add bolus to initial conditions

end

for i = 165:167
    
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

T2 = T2./28;

%% Dosing Regimen 3: 3/8 dose, 28 days on

params = [0.398*24 0.86 175*24 42000*24 2400*24 2.09*24 106*24 7.56*24 81.3*24 14.9*24 2.73*24 0.31*24 6000*24 1.81*24 30300 140 29 1421 2.5 0.0037*24 339.3 0.50 1 0.4 3 (6*10^10)*24 0]';
bolus3 = 1.3*(3/8);
                                                                                                                                                                                                                                          
% cycle 1
y0 = [bolus3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6]';

options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
[T3,Y3] = ode15s(@TMZ_equations,[0 1],y0,options,params);

dim = size(T3);

% doses 2-28 
y0 = Y3(end,:)';
y0(1) = y0(1) + bolus3;
totalD3 = repmat(bolus3,dim(1),1); % initial drug in

for i = 1:167
    
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

T3 = T3./28;

%% Plot tumor adduct over 6 cycles

tV1 = Y1(:,20)+Y1(:,21); % tumor volume (changes wrt time)
tV2 = Y2(:,20)+Y2(:,21);
tV3 = Y3(:,20)+Y3(:,21);

% plot tumor adduct
adductT = Y1(:,19)./tV1;
adductT2 = Y2(:,19)./tV2;
adductT3 = Y3(:,19)./tV3;

figure;
plot(T1,adductT,'b')
hold on
plot(T2,adductT2,'r')
plot(T3,adductT3,'g')
title('Tumor Adduct over 6 Months')
xlabel('Number of Treatment Months')
ylabel('Concentration (M)')
hold off

%% Plot tumor growth (volume) over 6 cycles

susT = Y1(:,20);
resT = Y1(:,21);
totalT = Y1(:,20)+Y1(:,21);

susT2 = Y2(:,20);
resT2 = Y2(:,21);
totalT2 = Y2(:,20)+Y2(:,21);

susT3 = Y3(:,20);
resT3 = Y3(:,21);
totalT3 = Y3(:,20)+Y3(:,21);

figure;
plot(T1,totalT,'b')
hold on
plot(T2,totalT2,'r')
plot(T3,totalT3,'g')
title('Tumor Growth over 6 Cycles')
xlabel('Number of Treatment Cycles')
ylabel('Volume (mL)')
hold off

%% Save data to .mat files for visualization - supp fig 2

save tumorPD_5on23off.mat T1 adductT resT susT totalT;
save tumorPD_7on8off.mat T2 adductT2 resT2 susT2 totalT2;
save tumorPD_28on.mat T3 adductT3 resT3 susT3 totalT3;