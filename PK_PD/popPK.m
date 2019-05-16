% problem scaffolding --> generate 100 random values for each of the
% particular parameters we want (given means and standard deviations)
% do this for MGMT population and silenced population
personIndex = [1:1:100]';

%% Growth Rate
% Growth Rate Varies

params_mgmt = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0*0.00026]';
params_silenced = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0.00026*0]';
meanGrowth = 0.089*28;
stdGrowth = 0.1494*28;


initpopulation = normrnd(meanGrowth, stdGrowth, 200, 1);
for val = 1:length(initpopulation)
    if initpopulation(val) < 0 % eliminating negative values
        while initpopulation(val) < 0
            initpopulation(val) = normrnd(meanGrowth, stdGrowth);
        end
    end
end

plasmaTMZ_AUC_mgmt = zeros(100,1);
tumorTMZ_AUC_mgmt = zeros(100,1);
tumorAdduct_AUC_mgmt = zeros(100,1);
CtroughAdduct_mgmt = zeros(100,1);
finalVol_mgmt = zeros(100,1);

for val = 1:100
    params_mgmt(20) = initpopulation(val);
    
    y0 = [1.3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6]';
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
    plasmaTMZ_AUC_mgmt(val) = trapz(T1, Y1(:,2)/params_mgmt(15));
    tumorTMZ_AUC_mgmt(val) = trapz(T1, (Y1(:,7)./(Y1(:,20)+Y1(:, 21))));
    tumorAdduct_AUC_mgmt(val) = trapz(T1, (Y1(:,19)./(Y1(:,20)+Y1(:,21))));
    CtroughAdduct_mgmt(val) = Y1(end,19)./(Y1(end,20)+Y1(end, 21));
    finalVol_mgmt(val) = Y1(end,20) + Y1(end, 21);
end

save growth_plasT_AUC_mgmt.mat plasmaTMZ_AUC_mgmt personIndex
save growth_tumT_AUC_mgmt.mat tumorTMZ_AUC_mgmt personIndex
save growth_tumorAdd_AUC_mgmt.mat tumorAdduct_AUC_mgmt personIndex
save growth_CtroughAdduct_mgmt.mat CtroughAdduct_mgmt personIndex
save growth_finalVol_mgmt.mat finalVol_mgmt personIndex

plasmaTMZ_AUC_silence = zeros(100,1);
tumorTMZ_AUC_silence = zeros(100,1);
tumorAdduct_AUC_silence = zeros(100,1);
CtroughAdduct_silence = zeros(100,1);
finalVol_silence = zeros(100,1);

for val2 = 1:100
    params_silenced(20) = initpopulation(val2 + 100);
    
    y0 = [1.3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';
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
    plasmaTMZ_AUC_silence(val2) = trapz(T1, (Y1(:,2)./params_silenced(15)));
    tumorTMZ_AUC_silence(val2) = trapz(T1, (Y1(:,7)./(Y1(:,20)+Y1(:, 21))));
    tumorAdduct_AUC_silence(val2) = trapz(T1, (Y1(:,19)./(Y1(:,20)+Y1(:, 21))));
    CtroughAdduct_silence(val2) = Y1(end,19)./(Y1(end,20)+Y1(end, 21));
    finalVol_silence(val2) = Y1(end,20) + Y1(end, 21);
end
save growth_plasT_AUC_silence.mat plasmaTMZ_AUC_silence personIndex
save growth_tumT_AUC_silence.mat tumorTMZ_AUC_silence personIndex
save growth_tumorAdd_AUC_silence.mat tumorAdduct_AUC_silence personIndex
save growth_CtroughAdduct_silence.mat CtroughAdduct_silence personIndex
save growth_finalVol_silence.mat finalVol_silence personIndex

%% Tumor Size
% Tumor Size Varies
params_mgmt = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0.00026]';
params_silenced = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0.00026*0]';

meanSize = 33.2;
stdSize = 29;


initpopulation = normrnd(meanSize, stdSize, 200, 1);
for val = 1:length(initpopulation)
    if initpopulation(val) < 0 % eliminating negative values
        while initpopulation(val) < 0
            initpopulation(val) = normrnd(meanSize, stdSize);
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
    y0(20) = initpopulation(val)/2;
    y0(21) = initpopulation(val)/2;
    
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
    plasmaTMZ_AUC_mgmt(val) = trapz(T1, (Y1(:,2)/params_mgmt(15)));
    tumorTMZ_AUC_mgmt(val) = trapz(T1, (Y1(:,7)./(Y1(:,20)+Y1(:, 21))));
    tumorAdduct_AUC_mgmt(val) = trapz(T1, (Y1(:,19)./(Y1(:,20)+Y1(:, 21))));
    CtroughAdduct_mgmt(val) = Y1(end,19)./(Y1(end,20)+Y1(end, 21));
    finalVol_mgmt(val) = Y1(end,20) + Y1(end, 21);
end

save size_plasT_AUC_mgmt.mat plasmaTMZ_AUC_mgmt personIndex
save size_tumT_AUC_mgmt.mat tumorTMZ_AUC_mgmt personIndex
save size_tumorAdd_AUC_mgmt.mat tumorAdduct_AUC_mgmt personIndex
save size_CtroughAdduct_mgmt.mat CtroughAdduct_mgmt personIndex
save size_finalVol_mgmt.mat finalVol_mgmt personIndex

plasmaTMZ_AUC_silence = zeros(100,1);
tumorTMZ_AUC_silence = zeros(100,1);
tumorAdduct_AUC_silence = zeros(100,1);
CtroughAdduct_silence = zeros(100,1);
finalVol_silence = zeros(100,1);

for val2 = 1:100    
    y0 = [1.3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';
    y0(20) = initpopulation(val2+100)/2;
    y0(21) = initpopulation(val2+100)/2;
    
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
    plasmaTMZ_AUC_silence(val2) = trapz(T1, (Y1(:,2)./params_silenced(15)));
    tumorTMZ_AUC_silence(val2) = trapz(T1, (Y1(:,7)./(Y1(:,20)+Y1(:, 21))));
    tumorAdduct_AUC_silence(val2) = trapz(T1, (Y1(:,19)./(Y1(:,20)+Y1(:, 21))));
    CtroughAdduct_silence(val2) = Y1(end,19)./(Y1(end,20)+Y1(end, 21));
    finalVol_silence(val2) = Y1(end,20) + Y1(end, 21);
end
save size_plasT_AUC_silence.mat plasmaTMZ_AUC_silence personIndex
save size_tumT_AUC_silence.mat tumorTMZ_AUC_silence personIndex
save size_tumorAdd_AUC_silence.mat tumorAdduct_AUC_silence personIndex
save size_CtroughAdduct_silence.mat CtroughAdduct_silence personIndex
save size_finalVol_silence.mat finalVol_silence personIndex

%% Alpha
% Alpha Varies

params_mgmt = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0.00026]';
params_silenced = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0.00026*0]';

meanA = 0.5;
stdA = 0.15;


initpopulation = normrnd(meanA, stdA, 200, 1);
for val = 1:length(initpopulation)
    if initpopulation(val) < 0 % eliminating negative values
        while initpopulation(val) < 0
            initpopulation(val) = normrnd(meanSize, stdSize);
        end
    end
end

plasmaTMZ_AUC_mgmt = zeros(100,1);
tumorTMZ_AUC_mgmt = zeros(100,1);
tumorAdduct_AUC_mgmt = zeros(100,1);
CtroughAdduct_mgmt = zeros(100,1);
finalVol_mgmt = zeros(100,1);

for val = 1:100
    params_mgmt(22) = initpopulation(val);
    
    y0 = [1.3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6 ]';
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
    plasmaTMZ_AUC_mgmt(val) = trapz(T1, (Y1(:,2)/params_mgmt(15)));
    tumorTMZ_AUC_mgmt(val) = trapz(T1, (Y1(:,7)./(Y1(:,20)+Y1(:, 21))));
    tumorAdduct_AUC_mgmt(val) = trapz(T1, (Y1(:,19)./(Y1(:,20)+Y1(:, 21))));
    CtroughAdduct_mgmt(val) = Y1(end,19)./(Y1(end,20)+Y1(end, 21));
    finalVol_mgmt(val) = Y1(end,20) + Y1(end, 21);
end

save a_plasT_AUC_mgmt.mat plasmaTMZ_AUC_mgmt personIndex
save a_tumT_AUC_mgmt.mat tumorTMZ_AUC_mgmt personIndex
save a_tumorAdd_AUC_mgmt.mat tumorAdduct_AUC_mgmt personIndex
save a_CtroughAdduct_mgmt.mat CtroughAdduct_mgmt personIndex
save a_finalVol_mgmt.mat finalVol_mgmt personIndex

plasmaTMZ_AUC_silence = zeros(100,1);
tumorTMZ_AUC_silence = zeros(100,1);
tumorAdduct_AUC_silence = zeros(100,1);
CtroughAdduct_silence = zeros(100,1);
finalVol_silence = zeros(100,1);

for val2 = 1:100
    params_silenced(22) = initpopulation(val2+100);
    
    y0 = [1.3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';
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
    plasmaTMZ_AUC_silence(val2) = trapz(T1, (Y1(:,2)./params_silenced(15)));
    tumorTMZ_AUC_silence(val2) = trapz(T1, (Y1(:,7)./(Y1(:,20)+Y1(:, 21))));
    tumorAdduct_AUC_silence(val2) = trapz(T1, (Y1(:,19)./(Y1(:,20)+Y1(:, 21))));
    CtroughAdduct_silence(val2) = Y1(end,19)./(Y1(end,20)+Y1(end, 21));
    finalVol_silence(val2) = Y1(end,20) + Y1(end, 21);
end
save a_plasT_AUC_silence.mat plasmaTMZ_AUC_silence personIndex
save a_tumT_AUC_silence.mat tumorTMZ_AUC_silence personIndex
save a_tumorAdd_AUC_silence.mat tumorAdduct_AUC_silence personIndex
save a_CtroughAdduct_silence.mat CtroughAdduct_silence personIndex
save a_finalVol_silence.mat finalVol_silence personIndex

%% Resfrac
% Resistant Fraction Varies

params_MGMT = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0.00026]';
params_silenced = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*28 30300 140 29 1421 2.5 0.00401*24*28 339.3 0.50 0.5 0.8 3/28 (6*10^10)*28*24 0 0.00026*0]';

meanRf = 0.5;
stdRf = 0.15;


initpopulation = normrnd(meanRf, stdRf, 200, 1);
for val = 1:length(initpopulation)
    if initpopulation(val) < 0 % eliminating negative values
        while initpopulation(val) < 0
            initpopulation(val) = normrnd(meanSize, stdSize);
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
    y0(20) = (1-initpopulation(val))*(y0(20)+y0(21));
    y0(21) = (initpopulation(val))*(y0(20)+y0(21));
    
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
    tumorAdduct_AUC_mgmt(val) = trapz(T1, (Y1(:,19)./(Y1(:,20)+Y1(:, 21))));
    CtroughAdduct_mgmt(val) = Y1(end,19)./(Y1(end,20)+Y1(end, 21));
    finalVol_mgmt(val) = Y1(end,20) + Y1(end, 21);
end

save rf_plasT_AUC_mgmt.mat plasmaTMZ_AUC_mgmt personIndex
save rf_tumT_AUC_mgmt.mat tumorTMZ_AUC_mgmt personIndex
save rf_tumorAdd_AUC_mgmt.mat tumorAdduct_AUC_mgmt personIndex
save rf_CtroughAdduct_mgmt.mat CtroughAdduct_mgmt personIndex
save rf_finalVol_mgmt.mat finalVol_mgmt personIndex

plasmaTMZ_AUC_silence = zeros(100,1);
tumorTMZ_AUC_silence = zeros(100,1);
tumorAdduct_AUC_silence = zeros(100,1);
CtroughAdduct_silence = zeros(100,1);
finalVol_silence = zeros(100,1);

for val2 = 1:100    
    y0 = [1.3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';
    y0(20) = initpopulation(val2+100)/2;
    y0(21) = initpopulation(val2+100)/2;
    
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
    plasmaTMZ_AUC_silence(val2) = trapz(T1, (Y1(:,2)./params_silenced(15)));
    tumorTMZ_AUC_silence(val2) = trapz(T1, (Y1(:,7)./(Y1(:,20)+Y1(:, 21))));
    tumorAdduct_AUC_silence(val2) = trapz(T1, (Y1(:,19)./(Y1(:,20)+Y1(:, 21))));
    CtroughAdduct_silence(val2) = Y1(end,19)./(Y1(end,20)+Y1(end, 21));
    finalVol_silence(val2) = Y1(end,20) + Y1(end, 21);
end
save rf_plasT_AUC_silence.mat plasmaTMZ_AUC_silence personIndex
save rf_tumT_AUC_silence.mat tumorTMZ_AUC_silence personIndex
save rf_tumorAdd_AUC_silence.mat tumorAdduct_AUC_silence personIndex
save rf_CtroughAdduct_silence.mat CtroughAdduct_silence personIndex
save rf_finalVol_silence.mat finalVol_silence personIndex

