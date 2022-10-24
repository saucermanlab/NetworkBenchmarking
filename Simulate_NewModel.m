%% If any inputs have changed, run the following section: 
clear all; clc; %close all; 

% % % model='Hypertrophy_model_update.xlsx';
model='RyallHypertrophy_3.xlsx';

% load('finalDrugOutputNetworkTargets.mat'); % Result from the webscraper NEED TO CLEAN UP WEBSCRAPER OUTPUT AND THEN CHANGE THIS TO BE READTABLE LIKE THE OTHER VARIABLES
warning off;
formattedReactions = table;
% Species/Reaction information from toy_model.xlsx or network, 'species'/'reactions' tab 
networkReactions = readtable(model, 'Sheet', 'reactions');
% Formats network reactions to only show the product/output. 
for i = 1:height(networkReactions)
    reaction = string(networkReactions{i,3});
    nodeOfReaction = extractAfter(reaction, '=>'); 
    formattedReactions{i,1} = strtrim(nodeOfReaction);
end
formattedReactions.Properties.VariableNames(1) = {'ReactionOutputNode'};
save('formattedReactions.mat', 'formattedReactions');

%% Generate the ODE and parameter files from NETFLUX
if exist([pwd '\ODEfun_fromNetflux.m'],'file') == 2
    delete('ODEfun_fromNetflux.m');
end
if exist([pwd '\ODEfun_loadParams_fromNetflux.m'],'file') == 2
    delete('ODEfun_loadParams_fromNetflux.m');
end
namepos = findstr('.xls', model); namestr = cellstr(model(1:namepos-1));

[specID,reactionIDs,~,paramList,ODElist,~,~] = util.xls2Netflux(namestr,model);
% % % commandODE = util.exportODE2(specID,paramList,ODElist);
[commandODE,commandPARAM,b] = util.exportODE(specID,paramList,ODElist,'ODEfun');
util.textwrite('ODEfun_fromNetflux.m',commandODE);
util.textwrite('ODEfun_loadParams_fromNetflux.m',commandPARAM);

%% Inputs for simulations/

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct

%% Inputs and initial simulations

inp = 1; %input to network; ISO=10; Stretch=15; PE=14
inp2 = 24; %PKA=84; Calcium=103; PI3K=83; Akt=16; CaMK=24

% Node parameters
[params,y0] = ODEfun_loadParams_fromNetflux;
[rpar,tau,ymax,speciesNames]=params{:};
w = rpar(1,:);
n = rpar(2,:);
EC50 = rpar(3,:);
%New input value
% % % w(inp)=0.1128;

rpar = [w;n;EC50];
params = {rpar,tau,ymax,speciesNames};
% Steady-state Control Simulation
tspan = [0 50]; options = [];
[t,y] = ode15s(@ODEfun_fromNetflux,tspan,y0,options,params);
yEnd = y(end,:)';
% Reset initial y values
y0 = real(yEnd); 

%New input values
% % % w(inp)=0.1128;
% % % ymax(inp2) = 0;
y0(inp2)=1; tau(inp2)=10^9;

rpar = [w;n;EC50];
params = {rpar,tau,ymax,speciesNames};
[t2,y2] = ode15s(@ODEfun_fromNetflux,tspan,y0,options,params); % Make sure y0 is correct here.
y2 = real(y2);
ySimEnd = y2(end,:);
% Plot results
figure; hold on
 leg = [];
 for i = [10 3 4 6 8 36]   
     plot(t2,y2(:,i)); ylim([0 1]);
     leg = [leg, speciesNames(i)];
 end
 legend(leg)

 %Write data to text file
X = categorical(speciesNames');
X = reordercats(X,string(X));
T_Fin = y2(end,:)'; T0 = y0;
T = table(X,T0,T_Fin,T_Fin-T0,'VariableNames',{'Species_name','T0','T_Fin','T_Diff'});
filename = strcat(char(speciesNames(inp2)),'_ISO','.txt');
writetable(T, filename)

 %% Mechanistic Subnetwork

 phen = 4;
 knockoutData_drug = table('Size', [(length(speciesNames)+1) 2], 'VariableTypes', {'string', 'double'});
 knockoutData_control = table('Size', [(length(speciesNames)+1) 2], 'VariableTypes', {'string', 'double'});

 %control knockouts
 % Node parameters
[params,y0] = ODEfun_loadParams_fromNetflux;
[rpar,tau,ymax,speciesNames]=params{:};
w = rpar(1,:);
n = rpar(2,:);
EC50 = rpar(3,:);
%New input value
% % % w(inp)=0.1128;

rpar = [w;n;EC50];
params = {rpar,tau,ymax,speciesNames};
% Steady-state Control Simulation
tspan = [0 50]; options = [];
[t,y] = ode15s(@ODEfun_fromNetflux,tspan,y0,options,params);
yEnd = y(end,:)';
y0 = real(yEnd); 

inputNodeW = num2cell(1:length(speciesNames)); % Nodes to test drug against
% simulate knockdowns
for r = 0:length(speciesNames)
    disp(num2str(r))
    ymax_Knockout = ymax;
    if r > 0
        ymax_Knockout(r) = 0;
    end

    params = {rpar,tau,ymax_Knockout,speciesNames};
    tspan = [0 50]; 
    options = []; 
    [t2,y2] = ode15s(@ODEfun_fromNetflux,tspan,y0,options,params); 

    % Calculate final values 
    if r == 0
        ySim = real(y2(:,phen))';
        [c2] = ySim;
        c2End_control = c2(end);
        c2End_knockout = c2(end);
        label = 'Control';
    else
        ySim = real(y2(:,phen))';
        [c2] = ySim;
        c2End_knockout = c2(end);
        label = speciesNames{r};
    end

    % Calculates the change in cell area from the knockout (in the
    % presence of stimulus alone) and control (also in presence
    % of stimulus alone)
    dataVector_control = {label, c2End_knockout-c2End_control};
    knockoutData_control(r+1, :) = dataVector_control;
end

 %treatment knockouts
 % simulate knockdowns
% % %  w(inp)=0.1128;
% % % ymax(inp2) = 0;
y0(inp2)=1; tau(inp2)=10^9;
for r = 0:length(speciesNames)
    disp(num2str(r));
    ymax_Knockout = ymax;
    if r > 0
        ymax_Knockout(r) = 0;
    end
    
    params = {rpar,tau,ymax_Knockout,speciesNames};
    tspan = [0 50]; 
    options = []; 
    [t2,y2] = ode15s(@ODEfun_fromNetflux,tspan,y0,options,params);

    % Calculate final values
    if r == 0
        T_Fin = y2(end,:)'; T0 = y2(1,:)';
        T_Diff = T_Fin-T0;
        ySim = real(y2(:,phen))';
        [c1] = ySim;
        c1End_control = c1(end);
        c1End_knockout = c1(end);
        label = 'Control';
    else
        ySim = real(y2(:,phen))';
        [c1] = ySim;
        c1End_knockout = c1(end);
        label = speciesNames{r};
    end
    
    % Calculates the change in cell area from the knockout (in the
    % presence of stimulus and drug) and control (also in presence
    % of stimulus and drug)
    dataVector = {label, c1End_knockout-c1End_control};
    knockoutData_drug(r+1, :) = dataVector;
end

% Calculate the 'difference of the difference' - This finds how the
% changes in cell area in each analysis differ from each other,
% resulting in a measurement of how each node contributes to the
% mechanism of the drug
knockoutData = knockoutData_drug{:,2}-knockoutData_control{:,2};

%% Figures

X = categorical(knockoutData_drug{:,1}); 
X = reordercats(X,string(X));

figure
bar(X,knockoutData_control{:,2})
% % % % % set(gca, 'xtick', 1:(length(speciesNames)+1));
set(gca, 'xticklabel', knockoutData_control{:,1},'FontSize', 14);
xtickangle(90)
ylabel(strcat('Change in ',{' '},char(speciesNames(phen)),' ((Stimulus+KD) - Stimulus)'))
title('No Drug; ISO')
xlabel('Knockdowns')

figure
bar(X,knockoutData_drug{:,2})
% % % % % set(gca, 'xtick', 1:(length(speciesNames)+1));
set(gca, 'xticklabel', knockoutData_drug{:,1},'FontSize', 14);
xtickangle(90)
ylabel(strcat('Change in ',{' '},char(speciesNames(phen)),' ((Stimulus+Drug+KD) - (Stimulus+Drug))'))
title(strcat(char(speciesNames(inp2)),'-i ; ISO'))
xlabel('Knockdowns')

%Filter by nodes that were affected by drug action
for ii=1:length(T_Diff)
    if abs(T_Diff(ii)) < 0.01
        knockoutData(ii+1) = 0;
    end
end

%Filter by effect size (identifies the most important nodes)
thresh = 0.1*max(abs(knockoutData));
for ii = 1:length(knockoutData)
    if abs(knockoutData(ii)) < thresh
        knockoutData(ii) = 0;
    end
end

figure
bar(X,knockoutData)
% % % % % set(gca, 'xtick', 1:(length(speciesNames)+1));
set(gca, 'xticklabel', knockoutData_drug{:,1},'FontSize', 14);
xtickangle(90)
ylabel(strcat('Difference of KD Effect on',{' '},char(speciesNames(phen))))
title(strcat(char(speciesNames(inp2)),' - No Drug; ISO'))
xlabel('Knockdowns')

%write data to text file
T = table(X,knockoutData,'VariableNames',{'Species_name','KDd'}); 
T(1,:) = []; filename = strcat('NewEdges_ISO',char(speciesNames(inp2)),'_KDdata.txt');
writetable(T, filename)