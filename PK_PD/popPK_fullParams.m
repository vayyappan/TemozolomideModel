clear all; %required to make sure don't mix up variable names!

personIndex = [1:1:100]';
params_mgmt = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0.00026*0]';
params_silenced = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0.00026*0]';

meanKa = 0.398*24*28;
stdKa = 0.398*35.2/100*24*28;
initpopulationKa = normrnd(meanKa, stdKa, 200, 1);
for val = 1:length(initpopulationKa)
    if initpopulationKa(val) < 0 % eliminating negative values
        while initpopulationKa(val) < 0
            initpopulationKa(val) = normrnd(meanKa, stdKa);
        end
    end
end

meanktm = 2.09*24*28;
stdktm = 2.09*0.0002*24*28;
initpopulationktm = normrnd(meanktm, stdktm, 200, 1);
for val = 1:length(initpopulationktm)
    if initpopulationktm(val) < 0 % eliminating negative values
        while initpopulationktm(val) < 0
            initpopulationktm(val) = normrnd(meanktm, stdktm);
        end
    end
end

meanktCL = 2.73*24*28;
stdktCL = 2.73*0.992*24*28;
initpopulationktCL = normrnd(meanktCL, stdktCL, 200, 1);
for val = 1:length(initpopulationktCL)
    if initpopulationktCL(val) < 0 % eliminating negative values
        while initpopulationktCL(val) < 0
            initpopulationktCL(val) = normrnd(meanktCL, stdktCL);
        end
    end
end

meankMC = 0.31*24*28;
stdkMC = 0.31*0.0026*24*28;
initpopulationkMC = normrnd(meankMC, stdkMC, 200, 1);
for val = 1:length(initpopulationkMC)
    if initpopulationkMC(val) < 0 % eliminating negative values
        while initpopulationkMC(val) < 0
            initpopulationkMC(val) = normrnd(meankMC, stdkMC);
        end
    end
end

meankDNA = 1.81*24*28;
stdkDNA = 1.81*0.33*24*28;
initpopulationDNA = normrnd(meankDNA, stdkDNA, 200, 1);
for val = 1:length(initpopulationDNA)
    if initpopulationDNA(val) < 0 % eliminating negative values
        while initpopulationDNA(val) < 0
            initpopulationDNA(val) = normrnd(meankDNA, stdkDNA);
        end
    end
end

meanGrowth = 0.089*28;
stdGrowth = 0.1494*28;
initpopulationGrowth = normrnd(meanGrowth, stdGrowth, 200, 1);
for val = 1:length(initpopulationGrowth)
    if initpopulationGrowth(val) < 0 % eliminating negative values
        while initpopulationGrowth(val) < 0
            initpopulationGrowth(val) = normrnd(meanGrowth, stdGrowth);
        end
    end
end

meanSize = 33.2;
stdSize = 29;
initpopulationSize = normrnd(meanSize, stdSize, 200, 1);
for val = 1:length(initpopulationSize)
    if initpopulationSize(val) < 0 % eliminating negative values
        while initpopulationSize(val) < 0
            initpopulationSize(val) = normrnd(meanSize, stdSize);
        end
    end
end

meanA = 0.5;
stdA = 0.15;
initpopulationA = normrnd(meanA, stdA, 200, 1);
for val = 1:length(initpopulationA)
    if initpopulationA(val) < 0 % eliminating negative values
        while initpopulationA(val) < 0
            initpopulationA(val) = normrnd(meanA, stdA);
        end
    end
end

meanRf = 0.5;
stdRf = 0.15;
initpopulationRf = normrnd(meanRf, stdRf, 200, 1);
for val = 1:length(initpopulationRf)
    if initpopulationRf(val) < 0 % eliminating negative values
        while initpopulationRf(val) < 0
            initpopulationRf(val) = normrnd(meanRf, stdRf);
        end
    end
end

plasmaTMZ_AUC_mgmt = zeros(100,1);
tumorTMZ_AUC_mgmt = zeros(100,1);
tumorAdduct_AUC_mgmt = zeros(100,1);
CtroughAdduct_mgmt = zeros(100,1);
finalVol_mgmt = zeros(100,1);

for val = 1:100
    y0 = [1.3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6 ]';
    
    params_mgmt(1) = initpopulationKa(val);
    params_mgmt(6) = initpopulationktm(val);
    params_mgmt(11) = initpopulationktCL(val);
    params_mgmt(12) = initpopulationkMC(val);
    params_mgmt(14) = initpopulationDNA(val);
    params_mgmt(20) = initpopulationGrowth(val);
    y0(20) = initpopulationSize(val)/2;
    y0(21) = initpopulationSize(val)/2; 
    params_mgmt(22) = initpopulationA(val);
    y0(20) = (1-initpopulationRf(val))*(y0(20)+y0(21));
    y0(21) = (initpopulationRf(val))*(y0(20)+y0(21));    
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0 1/28],y0,options,params_mgmt);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28 1/28+(i-1)*1/28],y0,options,params_mgmt);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
   
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5 (1/28)*5+23*(1/28)],y0,options,params_mgmt);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles
    for cycle = 1:5
        y0 = Y1(end,:)';
   
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28 cycle*(1/28)*28+(i)*(1/28)],y0,options,params_mgmt);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5 cycle*(1/28)*28+(1/28)*5+23*(1/28)],y0,options,params_mgmt);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
    plasmaTMZ_AUC_mgmt(val) = trapz(T1, (Y1(:,2)./params_mgmt(15)));
    tumorTMZ_AUC_mgmt(val) = trapz(T1, (Y1(:,7)./(Y1(:,20)+Y1(:, 21))));
    tumorAdduct_AUC_mgmt(val) = trapz(T1, (Y1(:,19)./(Y1(:,20)+Y1(:,21))));
    CtroughAdduct_mgmt(val) = Y1(end,19)/(Y1(end,20)+Y1(end, 21));
    finalVol_mgmt(val) = Y1(end,20) + Y1(end, 21);
end

save full_plasT_AUC_mgmt.mat plasmaTMZ_AUC_mgmt personIndex
save full_tumT_AUC_mgmt.mat tumorTMZ_AUC_mgmt personIndex
save full_tumorAdd_AUC_mgmt.mat tumorAdduct_AUC_mgmt personIndex
save full_CtroughAdduct_mgmt.mat CtroughAdduct_mgmt personIndex
save full_finalVol_mgmt.mat finalVol_mgmt personIndex

plasmaTMZ_AUC_silence = zeros(100,1);
tumorTMZ_AUC_silence = zeros(100,1);
tumorAdduct_AUC_silence = zeros(100,1);
CtroughAdduct_silence = zeros(100,1);
finalVol_silence = zeros(100,1);

for val2 = 1:100    
    y0 = [1.3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';
    
    params_mgmt(1) = initpopulationKa(val2+100);
    params_mgmt(6) = initpopulationktm(val2+100);
    params_mgmt(11) = initpopulationktCL(val2+100);
    params_mgmt(12) = initpopulationkMC(val2+100);
    params_mgmt(14) = initpopulationDNA(val2+100);
    params_mgmt(20) = initpopulationGrowth(val2+100);
    y0(20) = initpopulationSize(val2+100)/2;
    y0(21) = initpopulationSize(val2+100)/2; 
    params_mgmt(22) = initpopulationA(val2+100);
    y0(20) = (1-initpopulationRf(val2+100))*(y0(20)+y0(21));
    y0(21) = (initpopulationRf(val2+100))*(y0(20)+y0(21));    
    
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0 1/28],y0,options,params_silenced);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28 1/28+(i-1)*1/28],y0,options,params_silenced);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5 (1/28)*5+23*(1/28)],y0,options,params_silenced);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles
    for cycle = 1:5
        y0 = Y1(end,:)';
   
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28 cycle*(1/28)*28+(i)*(1/28)],y0,options,params_silenced);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5 cycle*(1/28)*28+(1/28)*5+23*(1/28)],y0,options,params_silenced);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
    plasmaTMZ_AUC_silence(val2) = trapz(T1, Y1(:,2)./params_silenced(15));
    tumorTMZ_AUC_silence(val2) = trapz(T1, Y1(:,7)./(Y1(:,20)+Y1(:, 21)));
    tumorAdduct_AUC_silence(val2) = trapz(T1, Y1(:,19)./(Y1(:,20)+Y1(:,21)));
    CtroughAdduct_silence(val2) = Y1(end,19)./(Y1(end,20)+Y1(end, 21));
    finalVol_silence(val2) = Y1(end,20) + Y1(end, 21);
end
save full_plasT_AUC_silence.mat plasmaTMZ_AUC_silence personIndex
save full_tumT_AUC_silence.mat tumorTMZ_AUC_silence personIndex
save full_tumorAdd_AUC_silence.mat tumorAdduct_AUC_silence personIndex
save full_CtroughAdduct_silence.mat CtroughAdduct_silence personIndex
save full_finalVol_silence.mat finalVol_silence personIndex

