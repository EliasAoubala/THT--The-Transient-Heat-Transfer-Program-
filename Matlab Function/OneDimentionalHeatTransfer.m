% ONEDIMENTIONALHEATTRANSFER 
% By Elias Aoubala, 2432537A, 18/06/2021
% 
% Correct formating of Function;
%   -> ONEDIMENTIONALHEATTRANSFER(Tint, Tamb, CombP, matArray, BC);
%
% The Purpose of this script is to act as a preliminary assesment of the
% Heat Transfer experienced due to combustion in multiple configurations.
%
% **** Variable Identification *****
%
% Tint     - This is the Initialisaiton variable which specifies the initial
%            temperature of the Composite wall at time t=0.
%
% Tamb     - This is the Ambient Temperature of the Enviroment form which
%            the heat transfer process takes place.
%
% CombP    - This is an Array which includes all required combustion parameters
%            for the bartz analysis.
%            [Pr, Tc, gamma, mue, Cp, Pc, c_star, At, dt, M_study, A_study, d_study]
%               * Where the prefix study implies the section of the chamber
%                 to be studied
%
% matArray - This is the material Array containining all the material
%            property parameters and boundaries. Should the material be
%            composite, a new row is created for each layer.
%            [name, n, L, K, rho, Cp, T_limit]
%            
% BC       - This specifies the Boundary Conditions of the problem to be 
%            analysed
%            [Fo, t_end, ht, t1, t2, t3, etc...]

function OneDimentionalHeatTransfer(Tint, Tamb, CombP, matArray, BC)

% First We must intialise and quantify the key values including:
%   - Number of Layers in Composite Wall
%   - Number of time periods to investigate
%   - Initial Evaluation of Heat transfer values such as the Taw

% Error Checking will be Initially performe to determine the viability of
% the input values

if error_test(Tint, CombP, matArray, BC)==1
    
end

% We must identify Now the number of layers in the wall
[n_layers,~]=size(matArray);

% We must identify Now the number of time periods of Investigation
[~,t_int]=size(BC);
t_no=t_int-3;
t_study=BC(4:t_int);

% Deriving Cold Side Coefficient
hl=BC(3);

% Initialising the Displacement Vector,
x=[0];

% We must now start calculating and extracting data and step distances for
% work
for i=1:n_layers
    dx(i)=str2double(matArray(i,3))/(str2double(matArray(i,2))-1);
    alpha(i)=str2double(matArray(i,4))/(str2double(matArray(i,5))*str2double(matArray(i,6)));
    Names(i)=matArray(i,1);
    T_limit(i)=str2double(matArray(i,7));
    % We must Now define the minium time step required in thi simulation to
    % satisfy the following conditions
    dt_int(i)=((dx(i)^2)*BC(1))/alpha(i);
    
    if i==1
        % Initial dt value defined
        dt=dt_int(i);
        
        % Initialise L- distance- value
        L_total=str2double(matArray(i,3));
        
        % Initialise the total Number of Points
        n=str2double(matArray(i,2));
        
        
    elseif dt_int(i)<dt_int(i-1)
        % Should new dt value be smaller, it is recalculated
        dt=dt_int(i);
        
        % Total Distance added to that matrix 
        L_total=L_total+str2double(matArray(i,3));
        
        % Add subsequent number of Nodes
        n=n+str2double(matArray(i,2))-1;
    else
        % Total Distance added to that matrix
        L_total=L_total+str2double(matArray(i,3));
        
        % Add subsequent number of Nodes
        n=n+str2double(matArray(i,2))-1;
    end
    
    x_limit(i)=L_total;
    % Now we iterate the displacement vector 
    x=[x,((L_total-str2double(matArray(i,3))+dx(i)):dx(i):L_total)];
end


% Fourier Numbers must now be calculated for the new layers
Fo=alpha.*dt./(dx.^2);

% Final End-time defined
t_end=BC(2);

% We need to to now extract all necessary variables from the Combustion
% Array

Pr=CombP(1); Tc=CombP(2); gamma=CombP(3); mue=CombP(4); Cp=CombP(5);
Pc=CombP(6); c_star=CombP(7); At=CombP(8); dt_throat=CombP(9); M_study=CombP(10);
A_study=CombP(11); d_study=CombP(12);

fprintf("%d \n", A_study);
fprintf("%d \n", M_study);

% Now we calculate the Hot Gas Heat Transfer Parameters required for
% operation
r_var=Pr^0.33;
Taw=Tc*((1+r_var*((gamma-1)/2)*M_study^2)/(1+((gamma-1)/2)*M_study^2));

% Initialising the main Temperature array and counter 
T_actv=Tint*ones(1,n);
k=1;

% Reinitialisation of n
n=0;

% Now we need to start the major loop for solving for temperature
for t=0:dt:t_end
    % We must now equilibriate all the middle points for the composite wall
    % For each componenet of the wall
    
    for i=1:n_layers
        
        if i==1
           if n_layers>1
               % Equilibriate "Middle Points"
                for m= 2:1:str2double(matArray(i,2))
                    % This will equilibrate all the temperature Nodes
                    T_actv(m)=Fo(i)*(T_actv(m+1) + T_actv(m-1)) + (1-2*Fo(i))*T_actv(m);
                end
           else
                for m= 2:1:(str2double(matArray(i,2))-1)
                    % This will equilibrate all the temperature Nodes
                    T_actv(m)=Fo(i)*(T_actv(m+1) + T_actv(m-1)) + (1-2*Fo(i))*T_actv(m);
                end

           end
            % Define the Number Points in total
            n=str2double(matArray(i,2));
            
        else
            % Defining Outer Boundaries
            if n_layers-i>1
                for m= n:1:(n+str2double(matArray(i,2))-1)
                    if m==n
                        % This will average the new and original Data Point
                        T_two=Fo(i)*(T_actv(m+1) + T_actv(m-1)) + (1-2*Fo(i))*T_actv(m);
                        T_actv(m)=(T_actv(m)+T_two)/2;
                    else
                        % This will equilibrate all the temperature Nodes
                        T_actv(m)=Fo(i)*(T_actv(m+1) + T_actv(m-1)) + (1-2*Fo(i))*T_actv(m);
                    end
                end
            else
                for m= n:1:(n+str2double(matArray(i,2))-2)
                    if m==n
                        % This will average the new and original Data Point
                        T_two=Fo(i)*(T_actv(m+1) + T_actv(m-1)) + (1-2*Fo(i))*T_actv(m);
                        T_actv(m)=(T_actv(m)+T_two)/2;
                    else
                        % This will equilibrate all the temperature Nodes
                        T_actv(m)=Fo(i)*(T_actv(m+1) + T_actv(m-1)) + (1-2*Fo(i))*T_actv(m);
                    end
                end
            end
            
            n=n+str2double(matArray(i,2))-1;
        end 
    end
    
    % **** Data Collection Section ****
    Temp_out(k,:)=T_actv;
    tout(k)=t;
    
    % **** Boundary Condition Enforcement ****
    q_dot_hot=Bartz(T_actv(1), Tc, gamma, M_study, dt_throat, mue, Cp, Pr, Pc, c_star, At, A_study, Taw);
    
    % Update Hot-Side Temperature
    T_actv(1)=(q_dot_hot*dx(1)/str2double(matArray(1,4)))+T_actv(2);
    
    % Update Cold-Side Temperature
    Q_dot_cold=hl*(Tamb-T_actv(n));
    T_actv(n)=T_actv(n-1)-(Q_dot_cold*dx(n_layers)/str2double(matArray(1,4)));
    
    
    % Increase Counter by Unity
    k=k+1;
    
end

% **** Data Processing ****
Graphing(x, Temp_out, t_study, t_no, T_limit, Names, x_limit, dt, n_layers)



end

function [Q_dot_hot] = Bartz(Twg, Tc, gamma, M, dt, mue, Cp, Pr, Pc, c_star, At, A_study, Taw)
% This function's primaty use is for calculating the active heat Transfer
% from the combustion process into the chamber walls using the Bartz
% Approximation for the Heat Trasfer

% Firstly, we calculate the sigma parameter using the given values
% provided
sigma=1/(((0.5*(Twg/Tc)*(1+((gamma-1)/2)*M^2) + 0.5)^0.68)*((1+((gamma-1)/2)*M^2)^0.12));

% We then calculate the Heat Transfer Coefficient using Bartz
hg=(0.026/(dt^0.2))*((Cp*mue^0.2)/(Pr^0.6))*((Pc/c_star)^0.8)*((At/A_study)^0.9)*sigma;



% We can now calculate the Heat Transfer into the wall
Q_dot_hot=hg*(Taw-Twg);

end

function Graphing(x, Temp_out, t_study, n, T_limit, Names, x_limit,dt, n_layers)
% This function will generate the final plot from the results calculated
% for the composite wall.
% x = the displacement across the wall
% Temp_out= final computed temperatures at specified nodes
% T_melt= array of Limiting Temperatures for Analysis
% Names = Array of Text Names of Materials used.

% We intially apply hold function on to the graph
hold on;

plot(x,Temp_out(1,:),'DisplayName', "Intial Temperature at t=0 s");

% Plot Main graphs for analysis including Temperature Limits x- limits
for i=1:n
    name_graph=sprintf("Temperature at t= %.2f s", t_study(i));
    plot(x,Temp_out(round(t_study(i)/dt),:), 'DisplayName', name_graph);
end

% Plot Temperature Limits and x_limits

for i=1:n_layers
    yline(T_limit(i),'--', strcat("Maximum Operating Temperature of ", Names(i)), 'DisplayName', strcat("Maximum Operating Temperature of ", Names(i)));
end


xlabel("Distance along thickness of wall (m)");
ylabel("Temperature (K)")
title_name=sprintf("Temperature Variation across Chamber wall");
title(title_name);


ylim manual

legend('show','AutoUpdate','off'); % Legend is displayed for attributes
for i=1:n_layers
    if i<n_layers
    xline(x_limit(i),'-.',strcat(Names(i), '-',Names(i+1), ' Material border'));
    end
end

end


function r = error_test(Tint, CombP, matArray, BC)

[~,k]=size(CombP);
[~,p]=size(matArray);
[~,l]=size(BC);

if Tint<=0
   error("A Temperature of 0 Kelvin is Unviable")    
elseif k<12 || k>12
    error("Insufficient Number of Input Combustion Conditions")
    
elseif p<7 || p>7
    error("Insufficient Number of Input Material Conditions")   
elseif l<4
    error("Insufficient Number of Input Boundary Conditions")    
else
    r=1;
end

end

