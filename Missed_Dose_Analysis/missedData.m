params_mgmt = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*30 30300 140 29 1421 2.5 0.0037*24*28 339.3 0.5 1 0.4 3/28 (6*10^10)*28*24 0]';
params_silenced = [0.398*24*28 0.86 175*24*28 42000*24*28 2400*24*28 2.09*24*28 106*24*28 7.56*24*28 81.3*24*28 14.9*24*28 2.73*24*28 0.31*24*28 6000*24*28 1.81*24*30 30300 140 29 1421 2.5 0.0037*24*28 339.3 0.5 1 0.4 3/28 (6*10^10)*28*24 0]';

%% Completely Missed Dose
% Dose is not taken at all

for val = 1
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6 ]';
   
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_mgmt);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_mgmt);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
   
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_mgmt);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles
    
    
    for cycle = 1:5
        y0 = Y1(end,:)';
        j = 1;
        if cycle == 3
            j = 0;
        end 
   
        for i = 1:5
            y0(1) = y0(1)+1.3*j;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28:0.001:cycle*(1/28)*28+(i)*(1/28)],y0,options,params_mgmt);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5:0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)],y0,options,params_mgmt);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end

T1m = T1;
Y1m = Y1;

for val2 = 1    
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';

    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_silenced);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_silenced);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_silenced);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles  
    
    for cycle = 1:5
        y0 = Y1(end,:)';
            j = 1;
        if cycle == 3
            j = 0;
        end 
        for i = 1:5
            y0(1) = y0(1)+1.3*j;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28:0.001:cycle*(1/28)*28+(i)*(1/28)],y0,options,params_silenced);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5:0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)],y0,options,params_silenced);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end

T1s = T1;
Y1s = Y1;

%% Cycle Taken at m/5

for val = 1
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6 ]';
 
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_mgmt);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_mgmt);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
   
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_mgmt);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles
    

    
    for cycle = 1:5
        y0 = Y1(end,:)';
        j = 1;
        q = 0;
        r = 0;
        m = 0;
    
        if cycle == 3
            j = 1;
            m = 1;
        end
    
        if cycle+1 == 3
            q = 1;
        end
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*28/(5*28)],y0,options,params_mgmt);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*28/(5*28)],y0,options,params_mgmt);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end

T2m = T1;
Y2m = Y1;

for val2 = 1    
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';

    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_silenced);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_silenced);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_silenced);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

% additional cycles
    
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 0;
        m = 0;
    end
    
    if cycle+1 == 3
        q = 1;
    end
    
    for cycle = 1:5
        y0 = Y1(end,:)';
    j = 1;
        q = 0;
        r = 0;
        m = 0;
    
        if cycle == 3
            j = 1;
            m = 1;
        end
    
        if cycle+1 == 3
            q = 1;
        end
   
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*28/(5*28)],y0,options,params_silenced);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*28/(5*28)],y0,options,params_silenced);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end
T2s = T1; 
Y2s = Y1;

%% Cycle Taken at 2m/5

for val = 1
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6 ]';

    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_mgmt);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_mgmt);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
   
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_mgmt);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles
    
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 0;
        m = 0;
    end
    
    if cycle+1 == 3
        q = 1;
    end
    
    for cycle = 1:5
        y0 = Y1(end,:)';
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 1;
        m = 1;
    end
    
    if cycle+1 == 3
        q = 1;
    end
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*2*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*2*28/(5*28)],y0,options,params_mgmt);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*2*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*2*28/(5*28)],y0,options,params_mgmt);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end

T3m = T1;
Y3m = Y1;

for val2 = 1   
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';

    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_silenced);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_silenced);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_silenced);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

% additional cycles
    
    
    for cycle = 1:5
        y0 = Y1(end,:)';
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 1;
        m = 1;
    end
    
    if cycle+1 == 3
        q = 1;
    end
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*2*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*2*28/(5*28)],y0,options,params_silenced);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*2*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*2*28/(5*28)],y0,options,params_silenced);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end
T3s = T1;
Y3s = Y1;

%% Cycle Taken at 3m/5

for val = 1
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6 ]'; 
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_mgmt);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_mgmt);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
   
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_mgmt);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles
    
    
    for cycle = 1:5
        y0 = Y1(end,:)';
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 1;
        m = 1;
    end
    
    if cycle+1 == 3
        q = 1;
    end
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*3*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*3*28/(5*28)],y0,options,params_mgmt);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*3*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*3*28/(5*28)],y0,options,params_mgmt);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end

T4m = T1;
Y4m = Y1;

for val2 = 1    
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';  
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_silenced);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_silenced);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_silenced);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

% additional cycles
    
    
    for cycle = 1:5
        y0 = Y1(end,:)';
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 1;
        m = 1;
    end
    
    if cycle+1 == 3
        q = 1;
    end
        
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*3*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*3*28/(5*28)],y0,options,params_silenced);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*3*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*3*28/(5*28)],y0,options,params_silenced);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end
T4s = T1;
Y4s = Y1;

%% Cycle Taken at 4m/5

for val = 1
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6 ]';  
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_mgmt);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_mgmt);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
   
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_mgmt);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles
    
    
    for cycle = 1:5
        y0 = Y1(end,:)';
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 1;
        m = 1;
    end
    
    if cycle+1 == 3
        q = 1;
    end
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*4*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*4*28/(5*28)],y0,options,params_mgmt);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*4*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*4*28/(5*28)],y0,options,params_mgmt);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end
T5m = T1;
Y5m = Y1;

for val2 = 1    
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';

    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_silenced);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_silenced);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_silenced);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

% additional cycles
    
    
    for cycle = 1:5
        y0 = Y1(end,:)';
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 1;
        m = 1;
    end
    
    if cycle+1 == 3
        q = 1;
    end
        for i = 1:5
            y0(1) = y0(1)+1.3;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*4*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*4*28/(5*28)],y0,options,params_silenced);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*4*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*4*28/(5*28)],y0,options,params_silenced);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end
T5s = T1;
Y5s = Y1; 

%% Double Cycle Later
for val = 1
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6 ]';   
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_mgmt);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_mgmt);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
   
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_mgmt);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles
    
    for cycle = 1:5
        y0 = Y1(end,:)';
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 0;
    end
    
    if cycle-1 == 3
        j = 2;
    end
        for i = 1:5
            y0(1) = y0(1)+1.3*j;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*28/(5*28)],y0,options,params_mgmt);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*28/(5*28)],y0,options,params_mgmt);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end
T6m = T1;
Y6m = Y1; 

for val2 = 1    
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';   
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_silenced);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_silenced);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_silenced);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

% additional cycles
    
    
    for cycle = 1:5
        y0 = Y1(end,:)';
    j = 1;
    q = 0;
    r = 0;
    m = 0;
    
    if cycle == 3
        j = 0;
    end
    
    if cycle-1 == 3
        j = 2;
    end
        for i = 1:5
            y0(1) = y0(1)+1.3*j;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28+m*28/(5*28):0.001:cycle*(1/28)*28+(i)*(1/28)+m*28/(5*28)],y0,options,params_silenced);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5+m*28/(5*28):0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)+q*28/(5*28)],y0,options,params_silenced);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end
T6s = T1;
Y6s = Y1;
%% Completely Missed Dose
% Dose is not taken at all

for val = 1
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 8*10^-6 ]';
   
    
    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_mgmt);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_mgmt);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
   
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_mgmt);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles
    
    
    for cycle = 1:5
        y0 = Y1(end,:)';
        j = 1;
        if cycle == 300
            j = 0;
        end 
   
        for i = 1:5
            y0(1) = y0(1)+1.3*j;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28:0.001:cycle*(1/28)*28+(i)*(1/28)],y0,options,params_mgmt);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5:0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)],y0,options,params_mgmt);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end

T7m = T1;
Y7m = Y1;

for val2 = 1    
    y0 = [1.300 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16.6 16.6 0 0 0 0*8*10^-6 ]';

    options = odeset('MaxStep',5e-2, 'AbsTol', 1e-5,'RelTol', 1e-5,'InitialStep', 1e-5);
    [T1,Y1] = ode15s(@TMZ_equations,[0:0.001:1/28],y0,options,params_silenced);
    % cycle 0 

    y0 = Y1(end,:)';
    y0(1) = y0(1)+1.3;
    for i = 2:5
    
        [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[0+(i-1)*1/28:0.001:1/28+(i-1)*1/28],y0,options,params_silenced);
   
        y0 = Y1s_new(end,:)';
        y0(1) = y0(1)+1.3;
   
        T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    
    end


    y0 = Y1(end,:)';
    [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[(1/28)*5:0.001:(1/28)*5+23*(1/28)],y0,options,params_silenced);
    T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
    Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
    T1 = T1s_new;
    Y1 = Y1s_new;

    % additional cycles  
    
    for cycle = 1:5
        y0 = Y1(end,:)';
            j = 1;
        if cycle == 300
            j = 0;
        end 
        for i = 1:5
            y0(1) = y0(1)+1.3*j;
            [T1s_new,Y1s_new] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(i-1)*1/28:0.001:cycle*(1/28)*28+(i)*(1/28)],y0,options,params_silenced);
            y0 = Y1s_new(end,:)';
            T1s_new = vertcat(T1,T1s_new(1:end-1,:)); %goes to end-1 to prevent repeats
            Y1s_new = vertcat(Y1,Y1s_new(1:end-1,:));
            T1 = T1s_new;
            Y1 = Y1s_new;
        end
        y0 = Y1(end,:)';
        [T1s_off,Y1s_off] = ode15s(@TMZ_equations,[cycle*(1/28)*28+(1/28)*5:0.001:cycle*(1/28)*28+(1/28)*5+23*(1/28)],y0,options,params_silenced);
 
        T1s_new = vertcat(T1,T1s_off(1:end-1,:)); %goes to end-1 to prevent repeats
        Y1s_new = vertcat(Y1,Y1s_off(1:end-1,:));
        T1 = T1s_new;
        Y1 = Y1s_new;
    end
end

T7s = T1;
Y7s = Y1;



%% Putting Everything Together

figure
subplot(3, 2, 1)
plot(T1m,Y1m(:,20) + Y1m(:,21), 'b', T2m, Y2m(:,21)+ Y1m(:,20), 'r', T3m, Y3m(:,21) + Y3m(:,20), 'g', T4m, Y4m(:,21) + Y4m(:,20), 'c', T5m, Y5m(:,21) + Y5m(:,20), 'm', T6m, Y6m(:,21) + Y6m(:,20), 'k', 'linewidth' , 2)
title('Tumor Volume Over 6 cycles - MGMT Expressed');
xlabel('time (Months)');
ylabel('Volume (mL)')
legend({'Missed', 'M/5', '2M/5', '3M/5', '4M/5', 'M'})
subplot(3, 2, 2)
plot(T1s,Y1s(:,20) + Y1s(:,21), 'b', T2s, Y2s(:,21)+ Y1s(:,20), 'r', T3s, Y3s(:,21) + Y3s(:,20), 'g', T4s, Y4s(:,21) + Y4s(:,20), 'c', T5s, Y5s(:,21) + Y5s(:,20), 'm', T6s, Y6s(:,21) + Y6s(:,20), 'k', 'linewidth' , 2)
title('Tumor Volume Over 6 cycles - MGMT Silenced');
xlabel('time (Months)');
ylabel('Volume (mL)')
legend({'Missed', 'M/5', '2M/5', '3M/5', '4M/5', 'M'})
subplot(3,2,3)
plot(T1m,Y1m(:,2), 'b', T2m,Y2m(:,2), 'r', T3m,Y3m(:,2), 'g', T4m,Y4m(:,2), 'c', T5m,Y5m(:,2), 'm', T6m,Y6m(:,2), 'k', 'linewidth' , 2)
title('Plasma Temozolomide (M) - MGMT Expressed');
xlabel('time (Months)');
ylabel('Concentration of Temozolomide (M)')
legend({'Missed', 'M/5', '2M/5', '3M/5', '4M/5', 'M'})
subplot(3,2,4)
plot(T1s,Y1s(:,2), 'b', T2s,Y2s(:,2), 'r', T3s,Y3s(:,2), 'g', T4s,Y4s(:,2), 'c', T5s,Y5s(:,2), 'm', T6s,Y6s(:,2), 'k', 'linewidth' , 2)
title('Plasma Temozolomide (M) - MGMT Silenced');
xlabel('time (Months)');
ylabel('Concentration of Temozolomide (M)')
legend({'Missed', 'M/5', '2M/5', '3M/5', '4M/5', 'M'})
subplot(3,2,5)
plot(T1m,Y1m(:,19)./(Y1m(:,20)+Y1m(:, 21)), 'b', T2m,Y2m(:,19)./(Y2m(:,20)+Y2m(:, 21)), 'r', T3m,Y3m(:,19)./(Y3m(:,20)+Y3m(:, 21)), 'g', T4m,Y4m(:,19)./(Y4m(:,20)+Y4m(:, 21)), 'c', T5m,Y5m(:,19)./(Y5m(:,20)+Y5m(:, 21)), 'm', T6m,Y6m(:,19)./(Y6m(:,20)+Y6m(:, 21)), 'k', 'linewidth' , 2)
title('Tumor DNA Adduct (M) - MGMT Expressed');
xlabel('time (Months)');
ylabel('Concentration of Adduct (M)')
legend({'Missed', 'M/5', '2M/5', '3M/5', '4M/5', 'M'})
subplot(3,2,6)
plot(T1s,Y1s(:,19)./(Y1s(:,20)+Y1s(:, 21)), 'b', T2s,Y2s(:,19)./(Y2s(:,20)+Y2s(:, 21)), 'r', T3s,Y3s(:,19)./(Y3s(:,20)+Y3s(:, 21)), 'g', T4s,Y4s(:,19)./(Y4s(:,20)+Y4s(:, 21)), 'c', T5s,Y5s(:,19)./(Y5s(:,20)+Y5s(:, 21)), 'm', T6s,Y6s(:,19)./(Y6s(:,20)+Y6s(:, 21)), 'k', 'linewidth' , 2)
title('Tumor DNA Adduct (M) - MGMT Silenced');
xlabel('time (Months)');
ylabel('Concentration of Adduct (M)')
legend({'Missed', 'M/5', '2M/5', '3M/5', '4M/5', 'M'})

tumVolM = [T1m, Y1m(:,20) + Y1m(:,21), Y2m(:,20) + Y2m(:,21), Y3m(:,20) + Y3m(:,21), Y4m(:,20) + Y4m(:,21), Y5m(:,20) + Y5m(:,21), Y6m(:,20) + Y6m(:,21), Y7m(:,20) + Y7m(:,21)];
tumVolS = [T1s, Y1s(:,20) + Y1s(:,21), Y2s(:,20) + Y2s(:,21), Y3s(:,20) + Y3s(:,21), Y4s(:,20) + Y4s(:,21), Y5s(:,20) + Y5s(:,21), Y6s(:,20) + Y6s(:,21), Y7s(:,20) + Y7s(:,21)];
plasTM = [T1m, Y1m(:,2), Y2m(:,2), Y3m(:,2), Y4m(:,2), Y5m(:,2), Y6m(:,2), Y7m(:,2)];
plasTS = [T1s, Y1s(:,2), Y2s(:,2), Y3s(:,2), Y4s(:,2), Y5s(:,2), Y6s(:,2), Y7s(:,2)];
tumDM = [T1m, Y1m(:,19)./(Y1m(:,20)+Y1m(:, 21)), Y2m(:,19)./(Y2m(:,20)+Y2m(:, 21)), Y3m(:,19)./(Y3m(:,20)+Y3m(:, 21)), Y4m(:,19)./(Y4m(:,20)+Y4m(:, 21)), Y5m(:,19)./(Y5m(:,20)+Y5m(:, 21)), Y6m(:,19)./(Y6m(:,20)+Y6m(:, 21)), Y7m(:,19)./(Y7m(:,20)+Y7m(:, 21))];
tumDS = [T1s, Y1s(:,19)./(Y1s(:,20)+Y1s(:, 21)), Y2s(:,19)./(Y2s(:,20)+Y2s(:, 21)), Y3s(:,19)./(Y3s(:,20)+Y3s(:, 21)), Y4s(:,19)./(Y4s(:,20)+Y4s(:, 21)), Y5s(:,19)./(Y5s(:,20)+Y5s(:, 21)), Y6s(:,19)./(Y6s(:,20)+Y6s(:, 21)), Y7s(:,19)./(Y7s(:,20)+Y7s(:, 21))];


save tumVolM.mat tumVolM
save tumVolS.mat tumVolS
save plasTM.mat plasTM
save plasTS.mat plasTS
save tumDM.mat tumDM
save tumDS.mat tumDS




