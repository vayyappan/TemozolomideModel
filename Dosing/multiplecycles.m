%% Dosing Regimen 1: Full dose, 5 days ON, 23 days OFF

bolus = 1.3;
params = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.0037*24*28 339.3 0.50 1 0.4 3/28 (6*10^10)*28*24 0]';

% cycle 1
y0 = [bolus 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6]';
options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
[T1,Y1] = ode15s(@TMZ_equations,[0 1/28],y0,options,params);

dim = size(T1);
totalD = repmat(bolus,dim(1),1); % initial drug in

for i = 1:4
    
    y0 = Y1(end,:)';
    y0(1) = y0(1) + bolus;
    
    tspanON = [(i*1/28) (1/28+i*1/28)];
    [T1s_new,Y1s_new] = ode15s(@TMZ_equations,tspanON,y0,options,params);
   
    dimON = size(T1s_new);
   
    onD = totalD(end)+bolus; % add bolus
    onD = repmat(onD,dimON(1),1);
    totalD = vertcat(totalD,onD(1:end-1,:));
   
    T1 = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1 = vertcat(Y1,Y1s_new(1:end-1,:));
   
end

y0 = Y1(end,:)';

tspanOFF = [((1/28)*5) ((1/28)*5+23*(1/28))];
[T1s_off,Y1s_off] = ode15s(@TMZ_equations,tspanOFF,y0,options,params);

dimOFF = size(T1s_off);
offD = repmat(totalD(end),dimOFF(1),1);
totalD = vertcat(totalD,offD(1:end-1,:));

T1 = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
Y1 = vertcat(Y1,Y1s_off(1:end-1,:));

% additional cycles
for cycle = 1:5 

    for i = 1:5
        
        y0 = Y1(end,:)';
        y0(1) = y0(1) + bolus;
        
        tspan_cyclesON = [(cycle+(i-1)*1/28) (cycle+(i)*(1/28))];

        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,tspan_cyclesON,y0,options,params);
        
        dimON = size(T1s_new);
   
        onD = totalD(end)+bolus; % add bolus
        onD = repmat(onD,dimON(1),1);
        totalD = vertcat(totalD,onD(1:end-1,:));

        T1 = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1 = vertcat(Y1,Y1s_new(1:end-1,:));

    end
    
    y0 = Y1(end,:)';
    tspan_cyclesOFF = [(cycle+(1/28)*5) (cycle+1)];

    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,tspan_cyclesOFF,y0,options,params);
    
    dimOFF = size(T1s_off);
    offD = repmat(totalD(end),dimOFF(1),1);
    totalD = vertcat(totalD,offD(1:end-1,:));

    T1 = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1 = vertcat(Y1,Y1s_off(1:end-1,:));

end

%% Mass balance 1

TMZm = Y1(:,1)+Y1(:,2)+Y1(:,3)+Y1(:,4)+Y1(:,5)+Y1(:,6)+Y1(:,7)+Y1(:,8);
MTICm = Y1(:,9)+Y1(:,10)+Y1(:,11)+Y1(:,12)+Y1(:,13)+Y1(:,14);
Cm = Y1(:,15)+Y1(:,16)+Y1(:,17);
DNAm = Y1(:,18)+Y1(:,19)+Y1(:,23)+Y1(:,24);
balance = totalD - (TMZm + MTICm + Cm + DNAm);

%% Plots

tV1 = Y1(:,20)+Y1(:,21); % tumor volume (changes wrt time)

% plot plasma TMZ
plasmaT = Y1(:,2)/params(15);

figure;
plot(T1,plasmaT)
title('Concentration of TMZ in Plasma')
xlabel('Number of Treatment Cycles')
ylabel('Concentration (M)')

% plot tumor TMZ
tumorT = Y1(:,7)./tV1;

figure;
plot(T1,tumorT)
title('Concentration of TMZ in Tumor')
xlabel('Dosing Cycle')
ylabel('Concentration (M)')

% plot tumor MTIC
tumorM = Y1(:,14)./tV1;

figure;
plot(T1,tumorM)
title('Concentration of MTIC in Tumor')
xlabel('Dosing Cycle')
ylabel('Concentration (M)')

% Mass balance
figure;
plot(T1,balance)
title('Mass Balance')
xlabel('Dosing Cycle')
ylabel('Balance (mmol)')

%% Save data to .mat files for visualization - supp fig 2

save multiplecycles_out.mat T1 plasmaT tumorT tumorM balance;
