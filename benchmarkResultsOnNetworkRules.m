%Read in the compare overlap spreadsheet and write the scoring results for
%Pathway Commons onto the network spreadsheet 
clear all
close all

load('summarizedIneractions.mat')

overallNums = zeros(3,4);

sheetNames = {'refNetwork_RyallHypertrophy',...
    'refNetwork_TanMechanosignaling',...
    'refNetwork_ZeiglerFibroblast'};
spreadsheetPaths = {'../RyallHypertrophy.xls',...
    '../TanMechanosignaling.xls',...
    '../ZeiglerFibroblast.xls'};
[overallNums(1,:), RyallResults] = ...
    merge(sheetNames{1}, spreadsheetPaths{1},...
    RyallNetworkInteractionsSummary);
[overallNums(2,:), TanResults] = ...
    merge(sheetNames{2}, spreadsheetPaths{2},...
    TanNetworkInteractionsSummary);
[overallNums(3,:), ZeiglerResults] = ...
    merge(sheetNames{3}, spreadsheetPaths{3},...
    ZeiglerNetworkInteractionsSummary);

[RyalledgeInts_merged] = checkSignor(...
    RyallResults, sheetNames{1});

[TanedgeInts_merged] = checkSignor(...
    TanResults, sheetNames{2});

[ZeigleredgeInts_merged] = checkSignor(...
    ZeiglerResults, sheetNames{3});
%%
lists = {'Cardiac Fibroblast' ...
    'Mechano-Signaling'...
    'Cardiac Hypertrophy'...
    };
figure('units', 'pixels', 'position', [0 0 1000 400]); 
c = barh(overallNums,'stacked');

grid on 
box off
yticklabels(lists);
yt = get(gca, 'YTick');
set(gca, 'FontSize', 14, 'FontName', 'Arial')

legend({'Functional Interaction','Physical Interaction', ...
    'Both', 'None'}, ...
    'FontSize', 14, 'Location', 'eastoutside', 'FontName', 'Arial');

title('Pathway Commons Annotations for Edges in Curated Network Models',...
    'FontSize', 16, 'FontName', 'Arial')
xlabel('Fraction of Total Interactions', 'FontSize', 14,...
    'FontName', 'Arial')
xlim([0 1])

c(1).FaceColor = [186 228 188]./255;
c(2).FaceColor = [123 204 196]./255;
c(3).FaceColor = [43 140 190]./255;
c(4).FaceColor = [240 249 232]./255;
%%
function [summaryNums, mergedResults] = merge(sheetName, spreadsheetPath, ...
    interactionSummary)
benchmarkResults = readtable('compareOverlap.xlsx', ...
    'Sheet', sheetName);

networkInts = readtable(spreadsheetPath, ...
    'Sheet', 'reactions');

mergedResults = table;
mergedResults.name = networkInts.ID;
mergedResults.rule = networkInts.Rule;
mergedResults.physicalPC = zeros(length(mergedResults.name),1);
mergedResults.functionalPC = zeros(length(mergedResults.name),1);
mergedResults.bothPC = zeros(length(mergedResults.name),1);
mergedResults.pubMedIDs = cell(length(mergedResults.name),1);
mergedResults.pubMedNum = zeros(length(mergedResults.name),1);
mergedResults.colorKey = zeros(length(mergedResults.name),1);

%%
ruleNames = cellfun(@(x) x(1:end-1), benchmarkResults.name, 'Un', 0);
%This for loop puts the physical and fucntional annotations in the merged
%Results table 
for i = 1:length(mergedResults.name)
    thisRule = mergedResults.name(i);
    
    present = sum(strcmp(thisRule, ruleNames));
    
    if present == 0
        mergedResults.colorKey(i) = -1;
    end
    
    interactions = find(strcmp(thisRule{1}, ruleNames));
    
    intScores = benchmarkResults(interactions,:);
    
    functionalScore = sum(cell2mat(cellfun(@(x) str2double(x), ...
        intScores.PathwayCommonsDirected, 'Un', 0)));
    physicalScore = sum(cell2mat(cellfun(@(x) str2double(x), ...
        intScores.PathwayCommonsUndirected, 'Un', 0)));
    bothScore = (physicalScore > 0) & (functionalScore > 0);
    
    mergedResults.bothPC(i) = bothScore;
    mergedResults.physicalPC(i) = physicalScore;
    mergedResults.functionalPC(i) = functionalScore;

    if bothScore > 0
        mergedResults.colorKey(i) = 3;
    else if physicalScore > 0
        mergedResults.colorKey(i) = 2;
        else if functionalScore > 0
            mergedResults.colorKey(i) = 1;
            
            end
        end
    end

end
%%
for i = 1:size(interactionSummary, 1)
%This for loop is going to put the pubMed ID's in the mergedResults table 

    inGene = interactionSummary.inGene{i};
    outGene = interactionSummary.outGene{i};
    
    forward = find(strcmp(inGene, benchmarkResults.inGene) & ...
        strcmp(outGene, benchmarkResults.outGene));
    reverse = find(strcmp(outGene, benchmarkResults.inGene) & ...
        strcmp(inGene, benchmarkResults.outGene));
    intIndex = unique([forward; reverse]);
    
    networkIndex = cellfun(@(x) find(strcmp(x, ...
        mergedResults.name)), ruleNames(intIndex), 'Un', 0);
    networkIndex = unique([networkIndex{:}]);
    
    for j = 1:length(networkIndex)
        annotations = unique([mergedResults.pubMedIDs{networkIndex(j)},...
            interactionSummary.PubMedIDs{i}]);
        mergedResults.pubMedIDs(networkIndex(j)) = {annotations};
        mergedResults.pubMedNum(networkIndex(j)) = length(annotations);
    end
end

summaryNums = zeros(1, 4);


functionalCount = sum(mergedResults.functionalPC > 0 &...
    mergedResults.bothPC == 0);
physicalCount = sum(mergedResults.physicalPC > 0 &...
    mergedResults.bothPC == 0);
bothCount = sum(mergedResults.bothPC > 0);

rxnCount = length(unique(ruleNames));

summaryNums(1:3) = [functionalCount, physicalCount, ...
    bothCount]./rxnCount;
summaryNums(4) = 1 - sum(summaryNums);

end


%Adding in Signor benchmarking 
function [edgeInts_merged] = checkSignor(...
    edgeInts_merged, sheetName)

    benchmarkResults = readtable('compareOverlap.xlsx', ...
        'Sheet', sheetName);
    geneInts_merged = benchmarkResults(:, 1:3);
    geneInts_merged.name = cellfun(@(x) x(1:end-1), geneInts_merged.name,...
        'Un', 0);
    
    database = readtable('Signor/SignorDatabases2.xlsx');

    intsKey = readtable('Signor/InteractionTypeKey.xlsx');

    edgeInts_merged.signorDirected = zeros(size(edgeInts_merged,1),1);
    edgeInts_merged.signorUndirected = zeros(size(edgeInts_merged,1),1);
    edgeInts_merged.signorDIRECT = zeros(size(edgeInts_merged,1),1);
    edgeInts_merged.signorPubMed = cell(size(edgeInts_merged,1),1);
    edgeInts_merged.signorSent = cell(size(edgeInts_merged,1),1);
    edgeInts_merged.matchPMID = zeros(size(edgeInts_merged,1),1);
    edgeInts_merged.combinedColorKey = zeros(size(edgeInts_merged,1),1);
    
    geneInts_merged.signorDirected = zeros(size(geneInts_merged,1),1);
    geneInts_merged.signorUndirected = zeros(size(geneInts_merged,1),1);
    geneInts_merged.signorDIRECT = zeros(size(geneInts_merged,1),1);
    geneInts_merged.signorPubMed = cell(size(geneInts_merged,1),1);
    geneInts_merged.signorSent = cell(size(geneInts_merged,1),1);
    
    database = innerjoin(database, intsKey);
    
    databaseDirected = database(database.DIRECTED == 1,:);
    databaseUndirected = database(database.DIRECTED == 0,:);
    %Loop through all interactions and score based on whether they're in
    %Signor as directed or undirected. Also produce a total score. Then, 
    %also copy over the sentence description + the direct column + PubMed 
    for iInt = 1:size(geneInts_merged,1)
        inGene = geneInts_merged.inGene(iInt);
        outGene = geneInts_merged.outGene(iInt);
        
        %Compare to directed and undirected networks 
        %Record presence, record PubMed ID, and sentence description  
        
        %Directed Matches
        directedMatches = intersect(find(strcmp(inGene{1}, ...
            databaseDirected.ENTITYA)), find(strcmp(outGene{1}, ...
            databaseDirected.ENTITYB)));
        if ~isempty(directedMatches)
            %Document info for the gene Matches
            geneInts_merged.signorDirected(iInt) = 1;
            
            geneInts_merged.signorPubMed(iInt) = ...
                {databaseDirected.PMID(directedMatches)};
            geneInts_merged.signorSent(iInt) = ...
                {databaseDirected.SENTENCE(directedMatches)};
            if sum(strcmp('YES', database.DIRECT(directedMatches))) > 0
                geneInts_merged.signorDIRECT(iInt) = 1;
            end  
        end
        %Undirected Matches
        undirectedMatches = [intersect(find(strcmp(inGene{1}, ...
            databaseUndirected.ENTITYA)), find(strcmp(outGene{1}, ...
            databaseUndirected.ENTITYB))) ; intersect(find(strcmp(inGene{1}, ...
            databaseUndirected.ENTITYA)), find(strcmp(outGene{1}, ...
            databaseUndirected.ENTITYB)))];
        if ~isempty(undirectedMatches)
            %Document info for the gene Matches
            geneInts_merged.signorUndirected(iInt) = 1;
            
            geneInts_merged.signorPubMed(iInt) = ...
                {unique([geneInts_merged.signorPubMed{iInt}; ...
                databaseUndirected.PMID(undirectedMatches)])};
            geneInts_merged.signorSent(iInt) = ...
                {unique([geneInts_merged.signorSent{iInt};...
                databaseDirected.SENTENCE(undirectedMatches)])};
            if sum(strcmp('YES', database.DIRECT(undirectedMatches))) > 0
                geneInts_merged.signorDIRECT(iInt) = 1;
            end  
        end
    end
   
    %Loop through the network edges and add the fields from the geneInts
    %table to the corresponding network edge 
    
    for iEdge = 1:size(edgeInts_merged,1)
        name = edgeInts_merged.name(iEdge);
    
        nameIndex = strcmp(name{1}, geneInts_merged.name);
        
        if sum(nameIndex) == 0
            edgeInts_merged.combinedColorKey(iEdge) = -1;
        end
        
        edgeInts_merged.signorDirected(iEdge) = ...
            sum(geneInts_merged.signorDirected(nameIndex));
        edgeInts_merged.signorUndirected(iEdge) = ...
            sum(geneInts_merged.signorUndirected(nameIndex));
        
        edgeInts_merged.signorDIRECT(iEdge) = ...
            sum(geneInts_merged.signorDIRECT(nameIndex));
        
        if (edgeInts_merged.signorDirected(iEdge) + ...
                edgeInts_merged.signorUndirected(iEdge)) > 0
            
            edgeInts_merged.signorPubMed(iEdge) = ...
                {unique(vertcat(geneInts_merged.signorPubMed{nameIndex}))};
            edgeInts_merged.signorSent(iEdge) = ...
                {unique(vertcat(geneInts_merged.signorSent{nameIndex}))};
        end
        
        if length(edgeInts_merged.pubMedIDs{iEdge}) > 0 && ...
                length(edgeInts_merged.signorPubMed{iEdge})
            matchingPMID = intersect(edgeInts_merged.pubMedIDs{iEdge}, ...
                edgeInts_merged.signorPubMed{iEdge});
        
            edgeInts_merged.matchPMID(iEdge) = length(matchingPMID);
        end
        
        if edgeInts_merged.functionalPC(iEdge) > 0 && ...
                edgeInts_merged.signorDirected(iEdge) > 0
            edgeInts_merged.combinedColorKey(iEdge) = 3;
        else if edgeInts_merged.functionalPC(iEdge) > 0
                edgeInts_merged.combinedColorKey(iEdge) = 2;
            else if edgeInts_merged.signorDirected(iEdge) > 0
                    edgeInts_merged.combinedColorKey(iEdge) = 1;
            end
            end
        end
        
        
    end
   
    
    
end