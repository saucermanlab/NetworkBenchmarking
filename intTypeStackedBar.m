%Load in the summarized interations for pathway commons and each netowrk to
%make the comparison figure with 3 categories - Percent Functional, percent
%physical, percent both 

%Plotting the total number of physical and total number of functional
%interactions in each category 
clear all 
close all

load('summarizedIneractions.mat')

lists = {'Cardiac Fibroblast' ...
    'Mechano-Signaling'...
    'Cardiac Hypertrophy'...
    ...%'CHSN Gene Interactions'...
    'Pathway Commons Database'...
    };
numbers = zeros(4,3);
numbers(4,:) = countTypes(databaseSummary);
% numbers(2,:) = [size(RyallGeneIntsFunctional,1),...
%     size(RyallGeneIntsPhysical,1)] ./ ...
%     (size(RyallGeneIntsFunctional,1) + size(RyallGeneIntsPhysical,1));
numbers(3,:) = countTypes(RyallNetworkInteractionsSummary);
numbers(2,:) = countTypes(TanNetworkInteractionsSummary);
numbers(1,:) = countTypes(ZeiglerNetworkInteractionsSummary);
RyallGeneNums = countTypes(RyallGeneInteractionsSummary);

figure('units', 'pixels', 'position', [0 0 1000 400]); 
c = barh(numbers,'stacked');

grid on 
box off
yticklabels(lists);
yt = get(gca, 'YTick');
set(gca, 'FontSize', 14, 'FontName', 'Arial')

legend({'Functional Interaction','Physical Interaction', 'Both'}, ...
    'FontSize', 14, 'Location', 'eastoutside', 'FontName', 'Arial');

title('Gene Product Interaction Annotations Across Curated Network Models',...
    'FontSize', 16, 'FontName', 'Arial')
xlabel('Fraction of Total Interactions', 'FontSize', 14,...
    'FontName', 'Arial')
xlim([0 1])

c(1).FaceColor = [186 228 188]./255;
c(2).FaceColor = [123 204 196]./255;
c(3).FaceColor = [43 140 190]./255;

%Also want to make a figure of 5 boxplots comparing the 
%distributions of PubMed IDs for each network + all PC + 
%Ryall Gene interactions 
x1 = databaseSummary.NumPubMedIDs(:);
x2 = RyallNetworkInteractionsSummary.NumPubMedIDs(:);
x3 = TanNetworkInteractionsSummary.NumPubMedIDs(:);
x4 = ZeiglerNetworkInteractionsSummary.NumPubMedIDs(:);
x5 = RyallGeneInteractionsSummary.NumPubMedIDs(:);

[~, ~, ib] = intersect(RyallNetworkInteractionsSummary(:,1:2),...
    RyallGeneInteractionsSummary(:,1:2),'rows');

catIdxAll = zeros((length(x2) + length(x5)),1);
catIdxAll(ib + (length(x2))) = 1;

%Generate tables that only contain interactions in 'both' section
network = RyallNetworkInteractionsSummary(...
    (RyallNetworkInteractionsSummary.totalFunctional > 0 & ...
    RyallNetworkInteractionsSummary.totalPhysical > 0),:);
genes = RyallGeneInteractionsSummary(...
    (RyallGeneInteractionsSummary.totalFunctional > 0 & ...
    RyallGeneInteractionsSummary.totalPhysical > 0),:);

[~, ~, ibBoth] = intersect(network(:,1:2),...
    genes(:,1:2),'rows');
catIdxBoth = zeros((length(x2Both) + length(x5Both)),1);
catIdxBoth(ibBoth + (length(x2Both))) = 1;

    
x2Both = network.NumPubMedIDs(:);
x5Both = genes.NumPubMedIDs(:);
x5Both(ibBoth) = [];

pubMedNums = [x1; x2; x3; x4; x5];

grouping = [zeros(length(x1),1); ones(length(x2),1); ...
    2 * ones(length(x3),1); 3 * ones(length(x4),1); ...
    4 * ones(length(x5Both),1)];

pubMedNums2 = [x2Both; x5Both];
grouping2 = [ones(length(x2Both),1); 4 * ones(length(x5Both),1)];

%Only CHSN gene interactions with both functional and physical annotations 
%CHSN gene interactions contains only ones not in the network already 
figure('units', 'pixels', 'position', [0 0 400 275]);
boxplot(pubMedNums2, grouping2, 'OutlierSize', .00001)
%plotSpread({x2Both, x5Both}, 'spreadWidth', .4)
ylim([-10 200])

grid on 
box off
label1 = '\begin{tabular}{c} CHSN \\ Interactions\end{tabular}';
label2 = '\begin{tabular}{c} PC Interactions \\ of CHSN Gene Products \\ Not Already In Network\end{tabular}';
set(gca,'xtick', [1 2], 'TickLabelInterpreter', 'latex');

xticklabels({label1,label2});

xt = get(gca, 'XTick');
set(gca, 'FontSize', 8, 'FontName', 'Arial')

yt = get(gca, 'YTick');
set(gca, 'FontSize', 8, 'FontName', 'Arial')

title('Number of PubMed IDs Supporting Interaction',...
    'FontSize', 10, 'FontName', 'Arial')

%{
%All CHSN Gene interactions 
pubMedNums2 = [x2; x5];
grouping2 = [ones(length(x2),1); 4 * ones(length(x5),1)];

figure('units', 'pixels', 'position', [0 0 400 275]);
boxplot(pubMedNums2, grouping2, 'OutlierSize', 0.01)
ylim([-10 200])

grid on 
box off
label1 = '\begin{tabular}{c} CHSN \\ Interactions\end{tabular}';
label2 = '\begin{tabular}{c} PC Interactions \\ of CHSN Gene Products\end{tabular}';
set(gca,'xtick', [1 2], 'TickLabelInterpreter', 'latex');

xticklabels({label1,label2});

xt = get(gca, 'XTick');
set(gca, 'FontSize', 8, 'FontName', 'Arial')

yt = get(gca, 'YTick');
set(gca, 'FontSize', 8, 'FontName', 'Arial')

title('Number of PubMed IDs Supporting Interactions',...
    'FontSize', 10, 'FontName', 'Arial')
%}
function results = countTypes(networkTable)

totalNum = size(networkTable, 1);

functionalFraction = length(find(networkTable.totalFunctional > 0 & ...
    networkTable.totalPhysical == 0)) ./ totalNum;
physicalFraction = length(find(networkTable.totalFunctional == 0 & ...
    networkTable.totalPhysical > 0)) ./ totalNum;
bothFraction = length(find(networkTable.totalFunctional > 0 & ...
    networkTable.totalPhysical > 0)) ./ totalNum;

results = [functionalFraction  physicalFraction, bothFraction];

end