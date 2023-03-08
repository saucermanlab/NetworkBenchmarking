%Take the pathway commons database and look at the interaction types for
%the hypertrophy network interactions, the int list produced when we use
%all netowrk genes, and the overall database. 

%importing Pathway Commons 
clear all 
close all

run ../../networkInfExpression/findClusterInteractions/importPathwayCommons;
database = database(1:1546602,:);
database.INTERACTION_TYPE = cellstr(database.INTERACTION_TYPE(:));
database.PARTICIPANT_A = cellstr(database.PARTICIPANT_A(:));
database.PARTICIPANT_B = cellstr(database.PARTICIPANT_B(:));

intTypes = readtable("InteractionTypeKey.xlsx");
database = outerjoin(database, intTypes,'Type','Left','MergeKeys',true);

toRemove = find(database.KEEP == 0);
database(toRemove,:) = [];

databaseFunctional = database((database.DIRECTED == 1),:);
databasePhysical = database((database.DIRECTED == 0),:);

%%
%Perform the tabulation for all of Pathway Commons
PCIntTypeCountFunctional = tabulate(databaseFunctional.INTERACTION_TYPE);
PCIntTypeCountPhysical = tabulate(databasePhysical.INTERACTION_TYPE);

%%

%Read in the interactions from Karen's network 
networkInts = readtable('../refNetwork_RyallHypertrophy.csv');

%Check for forward matches 
inMatchForward = cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_A)),...
        networkInts.inGene, 'Un', 0);
outMatchForward= cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_B)),...
        networkInts.outGene, 'Un', 0);
 
forward = cellfun(@(x, y) intersect(x, y), inMatchForward,...
    outMatchForward,'Un', 0);
forward = vertcat(forward{:});

%Check for reverse matches 
inMatchReverse = cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_A)),...
        networkInts.outGene, 'Un', 0);
outMatchReverse = cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_B)),...
        networkInts.inGene, 'Un', 0);
 
reverse = cellfun(@(x, y) intersect(x, y), inMatchReverse,...
    outMatchReverse,'Un', 0);
reverse = vertcat(reverse{:});

RyallNetworkInteractions = database(unique([forward; reverse]), :);
RyallNetworkInteractionsFunctional = ...
    RyallNetworkInteractions((RyallNetworkInteractions.DIRECTED == 1),:);
RyallNetworkInteractionsPhysical = ...
    RyallNetworkInteractions((RyallNetworkInteractions.DIRECTED == 0),:);

RyallNetworkInteractionCountsFunctional = ...
    tabulate(RyallNetworkInteractionsFunctional.INTERACTION_TYPE);
RyallNetworkInteractionCountsPhysical = ...
    tabulate(RyallNetworkInteractionsPhysical.INTERACTION_TYPE);
%%
%Do the same calucaltions for all Ryall network genes 
networkGenes = unique([networkInts.inGene(:) ; networkInts.outGene(:)]);

gene1Check = cellfun(@(x) find(strcmpi(x,networkGenes)),...
    database.PARTICIPANT_A, 'Un', 0);
match1 = cellfun(@(x) ~isempty(x), gene1Check, 'Un', 0);
gene1MatchIndex = find(cell2mat(match1));

gene2Check = cellfun(@(x) find(strcmpi(x,networkGenes)), ...
    database.PARTICIPANT_B, 'Un', 0);
match2 = cellfun(@(x) ~isempty(x), gene2Check, 'Un', 0);
gene2MatchIndex = find(cell2mat(match2));

overlapIndex = intersect(gene1MatchIndex, gene2MatchIndex);

RyallGeneInts = database(overlapIndex,:);
RyallGeneIntsFunctional = RyallGeneInts((RyallGeneInts.DIRECTED == 1),:);
RyallGeneIntsPhysical = RyallGeneInts((RyallGeneInts.DIRECTED == 0),:);



RyallGeneIntsCountsFunctional = ...
    tabulate(RyallGeneIntsFunctional.INTERACTION_TYPE);
RyallGeneIntsCountsPhysical = ...
    tabulate(RyallGeneIntsPhysical.INTERACTION_TYPE);
%%

%Read in the interactions from Phil's network 
networkInts = readtable('../refNetwork_TanMechanosignaling.csv');

%Check for forward matches 
inMatchForward = cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_A)),...
        networkInts.inGene, 'Un', 0);
outMatchForward= cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_B)),...
        networkInts.outGene, 'Un', 0);
 
forward = cellfun(@(x, y) intersect(x, y), inMatchForward,...
    outMatchForward,'Un', 0);
forward = vertcat(forward{:});

%Check for reverse matches 
inMatchReverse = cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_A)),...
        networkInts.outGene, 'Un', 0);
outMatchReverse = cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_B)),...
        networkInts.inGene, 'Un', 0);
 
reverse = cellfun(@(x, y) intersect(x, y), inMatchReverse,...
    outMatchReverse,'Un', 0);
reverse = vertcat(reverse{:});

TanNetworkInteractions = database(unique([forward; reverse]), :);
TanNetworkInteractionsFunctional = ...
    TanNetworkInteractions((TanNetworkInteractions.DIRECTED == 1),:);
TanNetworkInteractionsPhysical = ...
    TanNetworkInteractions((TanNetworkInteractions.DIRECTED == 0),:);

TanNetworkInteractionCountsFunctional = ...
    tabulate(TanNetworkInteractionsFunctional.INTERACTION_TYPE);
TanNetworkInteractionCountsPhysical = ...
    tabulate(TanNetworkInteractionsPhysical.INTERACTION_TYPE);
%%

%Read in the interactions from Angela's network 
networkInts = readtable('../refNetwork_ZeiglerFibroblast.csv');

%Check for forward matches 
inMatchForward = cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_A)),...
        networkInts.inGene, 'Un', 0);
outMatchForward= cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_B)),...
        networkInts.outGene, 'Un', 0);
 
forward = cellfun(@(x, y) intersect(x, y), inMatchForward,...
    outMatchForward,'Un', 0);
forward = vertcat(forward{:});

%Check for reverse matches 
inMatchReverse = cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_A)),...
        networkInts.outGene, 'Un', 0);
outMatchReverse = cellfun(@(x) find(strcmpi(x, database.PARTICIPANT_B)),...
        networkInts.inGene, 'Un', 0);
 
reverse = cellfun(@(x, y) intersect(x, y), inMatchReverse,...
    outMatchReverse,'Un', 0);
reverse = vertcat(reverse{:});

ZeiglerNetworkInteractions = database(unique([forward; reverse]), :);
ZeiglerNetworkInteractionsFunctional = ...
    ZeiglerNetworkInteractions((ZeiglerNetworkInteractions.DIRECTED == 1),:);
ZeiglerNetworkInteractionsPhysical = ...
    ZeiglerNetworkInteractions((ZeiglerNetworkInteractions.DIRECTED == 0),:);

ZeiglerNetworkInteractionCountsFunctional = ...
    tabulate(ZeiglerNetworkInteractionsFunctional.INTERACTION_TYPE);
ZeiglerNetworkInteractionCountsPhysical = ...
    tabulate(ZeiglerNetworkInteractionsPhysical.INTERACTION_TYPE);
%%
%Do the same calucaltions for all Zeigler network genes 
networkGenes = unique([networkInts.inGene(:) ; networkInts.outGene(:)]);

gene1Check = cellfun(@(x) find(strcmpi(x,networkGenes)),...
    database.PARTICIPANT_A, 'Un', 0);
match1 = cellfun(@(x) ~isempty(x), gene1Check, 'Un', 0);
gene1MatchIndex = find(cell2mat(match1));

gene2Check = cellfun(@(x) find(strcmpi(x,networkGenes)), ...
    database.PARTICIPANT_B, 'Un', 0);
match2 = cellfun(@(x) ~isempty(x), gene2Check, 'Un', 0);
gene2MatchIndex = find(cell2mat(match2));

overlapIndex = intersect(gene1MatchIndex, gene2MatchIndex);

ZeiglerGeneInts = database(overlapIndex,:);
ZeiglerGeneIntsFunctional = ZeiglerGeneInts((ZeiglerGeneInts.DIRECTED == 1),:);
ZeiglerGeneIntsPhysical = ZeiglerGeneInts((ZeiglerGeneInts.DIRECTED == 0),:);



ZeiglerGeneIntsCountsFunctional = ...
    tabulate(ZeiglerGeneIntsFunctional.INTERACTION_TYPE);
ZeiglerGeneIntsCountsPhysical = ...
    tabulate(ZeiglerGeneIntsPhysical.INTERACTION_TYPE);
%%
%Save these variables
save('GeneInteractionAnalysis_PathwayCommons.mat', ...
    'RyallGeneIntsFunctional', ...
    'RyallGeneIntsPhysical', ...
    'ZeiglerGeneIntsFunctional', ...
    'ZeiglerGeneIntsPhysical', ...
    'RyallNetworkInteractionsFunctional', ...
    'RyallNetworkInteractionsPhysical', ...
    'TanNetworkInteractionsFunctional', ...
    'TanNetworkInteractionsPhysical',...
    'ZeiglerNetworkInteractionsFunctional', ...
    'ZeiglerNetworkInteractionsPhysical',...
    'databaseFunctional', 'databasePhysical');
%%
%Match up all the cells for Functional Interactions 
[x,~] = size(PCIntTypeCountFunctional);

allFunctional = table;
allFunctional.labels = PCIntTypeCountFunctional(:,1);
allFunctional.PCdatabase = cell2mat(PCIntTypeCountFunctional(:,3));
allFunctional.RyallNetworkGenes = zeros(x,1);
allFunctional.RyallInteractions = zeros(x,1);
allFunctional.TanInteractions = zeros(x,1);
allFunctional.ZeiglerInteractions = zeros(x,1);
allFunctional = sortrows(allFunctional,'labels','ascend');

%Karen's network interactions 
RyallIntTableFunctional = cell2table(RyallNetworkInteractionCountsFunctional);
RyallIntTableFunctional.Properties.VariableNames = {'labels', 'counts', 'percent'};

all2Functional = outerjoin(allFunctional,RyallIntTableFunctional,...
    'Type','Left','MergeKeys',true);
allFunctional.RyallInteractions = all2Functional.percent;

%Phil's network interactions 
TanIntTableFunctional = cell2table(TanNetworkInteractionCountsFunctional);
TanIntTableFunctional.Properties.VariableNames = {'labels', 'counts', 'percent'};

all2Functional = outerjoin(allFunctional,TanIntTableFunctional,...
    'Type','Left','MergeKeys',true);
allFunctional.TanInteractions = all2Functional.percent;

%Angela's network interactions 
ZeiglerIntTableFunctional = cell2table(ZeiglerNetworkInteractionCountsFunctional);
ZeiglerIntTableFunctional.Properties.VariableNames = {'labels', 'counts', 'percent'};

all2Functional = outerjoin(allFunctional,ZeiglerIntTableFunctional,...
    'Type','Left','MergeKeys',true);
allFunctional.ZeiglerInteractions = all2Functional.percent;

%Karen's network gene query
RyallGeneTableFunctional = cell2table(RyallGeneIntsCountsFunctional);
RyallGeneTableFunctional.Properties.VariableNames = {'labels', 'counts', 'percent'};

all2Functional = outerjoin(allFunctional,RyallGeneTableFunctional,...
    'Type','Left','MergeKeys',true);
allFunctional.RyallNetworkGenes = all2Functional.percent;
%%
%Match up all the cells for Physical Interactions 
[x,~] = size(PCIntTypeCountPhysical);

allPhysical = table;
allPhysical.labels = PCIntTypeCountPhysical(:,1);
allPhysical.PCdatabase = cell2mat(PCIntTypeCountPhysical(:,3));
allPhysical.RyallNetworkGenes = zeros(x,1);
allPhysical.RyallInteractions = zeros(x,1);
allPhysical.TanInteractions = zeros(x,1);
allPhysical.ZeiglerInteractions = zeros(x,1);
allPhysical = sortrows(allPhysical,'labels','ascend');

%Karen's network interactions 
RyallIntTablePhysical = cell2table(RyallNetworkInteractionCountsPhysical);
RyallIntTablePhysical.Properties.VariableNames = {'labels', 'counts', 'percent'};

all2Physical = outerjoin(allPhysical,RyallIntTablePhysical,...
    'Type','Left','MergeKeys',true);
allPhysical.RyallInteractions = all2Physical.percent;

%Phil's network interactions 
TanIntTablePhysical = cell2table(TanNetworkInteractionCountsPhysical);
TanIntTablePhysical.Properties.VariableNames = {'labels', 'counts', 'percent'};

all2Physical = outerjoin(allPhysical,TanIntTablePhysical,...
    'Type','Left','MergeKeys',true);
allPhysical.TanInteractions = all2Physical.percent;

%Angela's network interactions 
ZeiglerIntTablePhysical = cell2table(ZeiglerNetworkInteractionCountsPhysical);
ZeiglerIntTablePhysical.Properties.VariableNames = {'labels', 'counts', 'percent'};

all2Physical = outerjoin(allPhysical,ZeiglerIntTablePhysical,...
    'Type','Left','MergeKeys',true);
allPhysical.ZeiglerInteractions = all2Physical.percent;

%Karen's network gene query
RyallGeneTablePhysical = cell2table(RyallGeneIntsCountsPhysical);
RyallGeneTablePhysical.Properties.VariableNames = {'labels', 'counts', 'percent'};

all2Physical = outerjoin(allPhysical,RyallGeneTablePhysical,...
    'Type','Left','MergeKeys',true);
allPhysical.RyallNetworkGenes = all2Physical.percent;
%%
%Plotting the total number of physical and total number of functional
%interactions in each category 
lists = {'Pathway Commons Database'...
    ...%'CHSN Gene Interactions'...
    'Cardiac Hypertrophy'...
    'Mechano-Signaling'...
    'Cardiac Fibroblast' ...
    };
numbers = zeros(4,2);
numbers(1,:) = [size(databaseFunctional,1), size(databasePhysical,1)] ./ ...
    (size(databaseFunctional,1) + size(databasePhysical,1));
% numbers(2,:) = [size(RyallGeneIntsFunctional,1),...
%     size(RyallGeneIntsPhysical,1)] ./ ...
%     (size(RyallGeneIntsFunctional,1) + size(RyallGeneIntsPhysical,1));
numbers(2,:) = [size(RyallNetworkInteractionsFunctional,1),...
    size(RyallNetworkInteractionsPhysical,1)] ./ ...
    (size(RyallNetworkInteractionsFunctional,1) + ...
    size(RyallNetworkInteractionsPhysical,1));
numbers(3,:) = [size(TanNetworkInteractionsFunctional,1),...
    size(TanNetworkInteractionsPhysical,1)] ./ ...
    (size(TanNetworkInteractionsFunctional,1) + ...
    size(TanNetworkInteractionsPhysical,1));
numbers(4,:) = [size(ZeiglerNetworkInteractionsFunctional,1),...
    size(ZeiglerNetworkInteractionsPhysical,1)] ./ ...
    (size(ZeiglerNetworkInteractionsFunctional,1) + ...
    size(ZeiglerNetworkInteractionsPhysical,1));

figure('units', 'pixels', 'position', [0 0 1000 400]); 
c = barh(numbers,'stacked');

grid on 
box off
yticklabels(lists);
yt = get(gca, 'YTick');
set(gca, 'FontSize', 14, 'FontName', 'Arial')

legend({'Functional Interactions','Physical Interactions'}, ...
    'FontSize', 14, 'Location', 'eastoutside', 'FontName', 'Arial');

title('Fraction of Functional vs Physical Interactions Across Networks',...
    'FontSize', 16, 'FontName', 'Arial')
xlabel('Fraction of Total Interactions',...
    'FontSize', 14,...
    'FontName', 'Arial')
xlim([0 1])

colorMap = [0 0 0; 1 1 1];

c(1).FaceColor = colorMap(1,:);
c(2).FaceColor = colorMap(2,:);


%%
%Plotting the results - grouped bar graph 
%Functional Interactions 
figure('units', 'pixels', 'position', [0 0 800 800]); 
hold on 

labels = allFunctional.labels';
numbers = [allFunctional.PCdatabase, ...
    ...#allFunctional.RyallNetworkGenes,...
    allFunctional.RyallInteractions, ...
    allFunctional.TanInteractions,...
    allFunctional.ZeiglerInteractions];
lists = {'Pathway Commons Database'...
    ...%'CHSN Gene Interactions'...
    'Cardiac Hypertrophy'...
    'Mechano-Signaling'...
    'Cardiac Fibroblast' ...
    };

functionalBar = barh(numbers,'grouped');

grid on 
box off
yticks([1:1:length(labels)])
yticklabels(labels)
yt = get(gca, 'YTickLabel');
set(gca, 'FontSize', 14, 'FontName', 'Arial')

title('Functional Interaction Types Across Networks',...
    'FontSize', 16, 'FontName', 'Arial')
xlabel('Percent of All Functional Interactions',...
    'FontSize', 14,...
    'FontName', 'Arial')

functionalBar(4).FaceColor = [0 0 0];
functionalBar(3).FaceColor = [.25 .25 .25];
functionalBar(2).FaceColor = [.5 .5 .5];
functionalBar(1).FaceColor = [.75 .75 .75];
%%
%Physical Interactions
figure('units', 'pixels', 'position', [0 0 1000 800]); 
hold on 
grid on
box off

labels = allPhysical.labels';
numbers = [allPhysical.PCdatabase, ...
    ...#allPhysical.RyallNetworkGenes,...
    allPhysical.RyallInteractions, ...
    allPhysical.TanInteractions,...
    allPhysical.ZeiglerInteractions];
physicalBar = barh(numbers,'grouped');
yticks([1:1:length(labels)])

yticklabels(labels)
yt = get(gca, 'YTickLabel');
set(gca, 'FontSize', 14, 'FontName', 'Arial')


legend(lists, ...
    'FontSize', 14, ...
    'Location', 'eastoutside', ...
    'FontName', 'Arial');

title('Physical Interaction Types Across Networks',...
    'FontSize', 16, 'FontName', 'Arial')
xlabel('Percent of All Physical Interactions',...
    'FontSize', 14,...
    'FontName', 'Arial')

physicalBar(4).FaceColor = [0 0 0];
physicalBar(3).FaceColor = [.25 .25 .25];
physicalBar(2).FaceColor = [.5 .5 .5];
physicalBar(1).FaceColor = [.75 .75 .75];
%legend(lists, 'Location','best');

%%
%Subplotting the 3 groups seaprately (will be a subplot with 2 rows and 3
%columns )
figure 
hold on 

%The top row will be the functional interactions 
labels = allFunctional.labels';
numbers = [allFunctional.PCdatabase, ...
    ...%allFunctional.RyallNetworkGenes, ...
    allFunctional.RyallInteractions, ...
    allFunctional.TanInteractions,...
    allFunctional.ZeiglerInteractions];
lists = {'Pathway Commons Database'...
    ...%'CHSN Gene Interactions'...
    'Cardiac Hypertrophy'...
    'Mechano-Signaling'...
    'Cardiac Fibroblast' ...
    };

%lists = allFunctional.Properties.VariableNames(2:end);

for i = 1:length(lists)
    subplot(2,5,i)
    
    myBar = barh(numbers(:,i),  'FaceColor', 'black');
    grid on 
    box off
    if i == 1
        yticklabels(labels)
    else
        set(gca,'yticklabel',[])
    end

    title(['Functional Interactions: ' lists(i)])
    xlabel('Percent of Total Interactions')
    xlim([0 100])
end


%The bottom row will be the Physical interactions 
labels = allPhysical.labels';
numbers = [allPhysical.PCdatabase, ...
    ...%allPhysical.RyallNetworkGenes, ...
    allPhysical.RyallInteractions, ...
    allPhysical.TanInteractions,...
    allPhysical.ZeiglerInteractions];
lists = {'Pathway Commons Database'...
    ...%'CHSN Gene Interactions'...
    'Cardiac Hypertrophy'...
    'Mechano-Signaling'...
    'Cardiac Fibroblast' ...
    };

for i = 1:length(lists)
    subplot(2,5,i + 5)
    
    myBar = barh(numbers(:,i), 'FaceColor', 'black');
    grid on 
    box off
    if i == 1
        yticklabels(labels)
    else
        set(gca,'yticklabel',[])
    end

    title(['Physical Interactions: ' lists(i)])
    xlabel('Percent of Total Interactions')
    xlim([0 100])
end

