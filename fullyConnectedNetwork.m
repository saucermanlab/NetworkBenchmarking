%This script will be to edge benchmark the expanded interactions list

%Step 1 - generate the fully connected networks 
%fullyConnectedEdges
networkSymbolsPath = '../hypertrophyNetworkGeneSymbols.xlsx';
fullyConnectedEdges = generateNetwork(networkSymbolsPath);

%Step 2- convert the format so that it is gene interactions that can be
%benchmarked 
%convertFormat
fullyConnectedGenes = convertFormat(fullyConnectedEdges,...
    networkSymbolsPath, ...
    'refNetwork_RyallHypertrophyFullyConnected.csv');

%Step 3 - put the benchmarking results on the fully connected network 
%merge
load('summarizedIneractions.mat')
[RyallGeneNums, geneInts_merged, edgeInts_merged] = ...
    merge(fullyConnectedGenes, fullyConnectedEdges, ...
    RyallGeneInteractionsSummary);

%Step 4 - mark which interactions are already in the network 
networkGeneInts = readtable(...
    '../refNetwork_RyallHypertrophy.csv');
[geneInts_merged, edgeInts_merged] = markNetworkInts(...
    edgeInts_merged, geneInts_merged, ...
    networkGeneInts);

%Step 5 - count how many times each PubMed ID occurs to look into
%establishing this is a method for eliminating high throughput studies 
allPubMedCounts = tabulate([edgeInts_merged.pubMedIDs{:}]);

[~,index] = sort(cell2mat(allPubMedCounts(:,2)), 'descend');

allPubMedCounts = allPubMedCounts(index,:);

for i = 1:size(allPubMedCounts,1)
    thisID = allPubMedCounts(i,1);
    
    matches = cellfun(@(x) find(strcmp(x,thisID)),...
        edgeInts_merged.pubMedIDs, 'Un', 0);
    
    indexMatches = ...
        find(cell2mat(cellfun(@(x) ~isempty(x), matches, 'Un', 0)));
    
    allPubMedCounts(i, 4) = {edgeInts_merged.rule(indexMatches)};

end

occurCounts = tabulate(cell2mat(allPubMedCounts(:,2)));

%Step 6 - add in Signor benchmarking 
[geneInts_merged, edgeInts_merged] = checkSignor(...
    edgeInts_merged, geneInts_merged);

%Step 7 - look at the number of times that the JAK -> ERK12 citations are
%used - is there a way we can filter? 
index = find(strcmp('r5968', edgeInts_merged.name));
JAKERK12_PMIDs = edgeInts_merged.pubMedIDs(index);
JAKERK12_PMID = table;
JAKERK12_PMID.ID = transpose(JAKERK12_PMIDs{:});
JAKERK12_PMID.count = zeros(size(JAKERK12_PMID,1),1);

for i = 1:size(JAKERK12_PMID,1)
    thisID = JAKERK12_PMID.ID(i);
    
    index = find(strcmp(thisID{1}, allPubMedCounts(:,1)));
    
    JAKERK12_PMID.count(i) = cell2mat(allPubMedCounts(index,2));
    
end

%Make figures

%1 - making boxplots comparing number of PubMed IDs of network edges to
%others
networkEdgeInts = edgeInts_merged(edgeInts_merged.inNetwork == 1,:);
otherEdgeInts = edgeInts_merged(edgeInts_merged.inNetwork == 0,:);

x1 = networkEdgeInts.pubMedNum(:);
x2 = otherEdgeInts.pubMedNum(:);

pubMedNums2 = [x1; x2];
grouping2 = [ones(length(x1),1); 4 * ones(length(x2),1)];

figure('units', 'pixels', 'position', [0 0 400 275]);
boxplot(pubMedNums2, grouping2, 'OutlierSize', .00001)

ylim([-10 200])

grid on 
box off
label1 = 'Edges in Network';
label2 = 'Candidate Edges Not in Network';
set(gca,'xtick', [1 2]);

xticklabels({label1,label2});

xt = get(gca, 'XTick');
set(gca, 'FontSize', 8, 'FontName', 'Arial')

yt = get(gca, 'YTick');
set(gca, 'FontSize', 8, 'FontName', 'Arial')

title('Number of PubMed IDs Supporting Interaction',...
    'FontSize', 10, 'FontName', 'Arial')

%2 - a plot to help figuring out which PubMed IDs to consider
figure('units', 'pixels', 'position', [0 0 400 275]);

plot(occurCounts(:,1), occurCounts(:,3), '*k')

title('Number of PubMed IDs Supporting "x" Number of Edges',...
    'FontSize', 10, 'FontName', 'Arial')

xlabel('Number of Edges Supported',...
    'FontSize', 10, 'FontName', 'Arial')
ylabel('Percent of PubMed IDs',...
    'FontSize', 10, 'FontName', 'Arial')

%(function for step 1)
function fullyConnectedEdges = generateNetwork(networkSymbolsPath)

networkSymbols = readtable(networkSymbolsPath);
symbols = networkSymbols.node(:);

fullyConnectedEdges = table;
fullyConnectedEdges.ID = cell(length(symbols).^2,1);
fullyConnectedEdges.Rule = cell(length(symbols).^2,1);

letter = cell(length(symbols).^2,1);
letter(:) = {'r'};
numbers = cellstr(num2str(transpose((1:length(symbols)^2))));
rules = strcat(letter, numbers);
rules = cellfun(@(x) strrep(x, ' ', ''), rules, 'Un', 0);

fullyConnectedEdges.ID = rules;

for node1 = 1:length(symbols)
    firstNode = symbols(node1);
    
    for node2 = 1:length(symbols)
        secondNode = symbols(node2);
        fullyConnectedEdges.Rule((node1 - 1) * length(symbols) + ...
            node2) = ...
            {strcat(firstNode{1},  " => ", secondNode{1})};
    end
end
end

%(function for step 2)
function refNetwork = ...
    convertFormat(networkReactions, networkSymbolsPath, outputName)

networkSymbols = readtable(networkSymbolsPath);

rules = networkReactions.Rule;

rules = cellfun(@(x) strrep(x, '!', ''), rules, 'Un', 0);

rules = cellfun(@(x) strsplit(x, ' => '), rules, 'Un', 0);

size = cellfun(@(x) numel(x), rules, 'Un', 0);
size = find(cell2mat(size) == 1);
rules(size) = [];
networkReactions(size,:) = [];

refNetwork = table;
refNetwork.name = {};
refNetwork.inGene = {};
refNetwork.outGene = {};

x = length(rules);
%loop through each of the rules. Create a new array of inGene, 
%outGene, and rules 

for iRule = 1:x
    thisRule = rules{iRule};
    
    %Get the first node
    node1 = thisRule(1);
    node1 = strsplit(node1{1}, ' ');
    
    
    location1 = cellfun(@(x) strcmp(x, networkSymbols.node), ...
        node1, 'Un', 0);
    location1 = cellfun(@(x) find(x), location1, 'Un', 0);
    
    genes1 = networkSymbols.genesymbol(cell2mat(location1),1);
    genes1 = cellfun(@(x) strsplit(x, ','), genes1, 'Un', 0);
    genes1 = [genes1{:}];
    
    %Get the second node
    node2 = thisRule(2);
    node2 = strsplit(node2{1}, ' ');
    
    location2 = cellfun(@(x) strcmp(x, networkSymbols.node), ...
        node2, 'Un', 0);
    location2 = cellfun(@(x) find(x), location2, 'Un', 0);
    
    genes2 = networkSymbols.genesymbol(cell2mat(location2),1);
    genes2 = cellfun(@(x) strsplit(x, ','), genes2, 'Un', 0);
    genes2 = [genes2{:}];
    
    thisReaction = networkReactions.ID(iRule,1);
    
    for iGene = 1:numel(genes1)
        
        thisGene1 = genes1(iGene);
        
        for jGene = 1:numel(genes2)
            thisGene2 = genes2(jGene);
            refNetwork.name(end + 1,1) = thisReaction;
            refNetwork.inGene(end,1) = thisGene1;
            refNetwork.outGene(end,1) = thisGene2;
        end
    
    end
    
end

delete1 = cellfun(@(x) strcmp(x, ''), refNetwork.inGene, 'Un', 0);
delete1 = find(cell2mat(delete1) == 1);
refNetwork(delete1,:) = [];

delete2 = cellfun(@(x) strcmp(x, ''), refNetwork.outGene, 'Un', 0);
delete2 = find(cell2mat(delete2) == 1);
refNetwork(delete2,:) = [];

writetable(refNetwork, outputName);
end

%(function for step 3)
function [summaryNums, geneInts, mergedResults] = ...
    merge(geneInts, networkInts, ...
    interactionSummary)

mergedResults = table;
mergedResults.name = networkInts.ID;
mergedResults.rule = networkInts.Rule;
mergedResults.physicalPC = zeros(length(mergedResults.name),1);
mergedResults.functionalPC = zeros(length(mergedResults.name),1);
mergedResults.bothPC = zeros(length(mergedResults.name),1);
mergedResults.totalScore = zeros(length(mergedResults.name),1);
mergedResults.pubMedIDs = cell(length(mergedResults.name),1);
mergedResults.pubMedNum = zeros(length(mergedResults.name),1);
mergedResults.colorKey = zeros(length(mergedResults.name),1);

geneInts.functionalPC = zeros(size(geneInts,1),1);
geneInts.physicalPC = zeros(size(geneInts,1),1);
geneInts.bothPC = zeros(size(geneInts,1),1);
geneInts.totalScore = zeros(size(geneInts,1),1);

geneInts.PubMedIDs = cell(size(geneInts,1),1);

countDoubles = [];
%%
%First part = I need to go through and determine which interactions are
%present in the pathway commons table
for i = 1:size(interactionSummary,1)
    inGene = interactionSummary.inGene(i);
    outGene = interactionSummary.outGene(i);
    
    %Check for matches in geneInts
    matches = intersect(find(strcmp(inGene, geneInts.inGene)), ...
        find(strcmp(outGene, geneInts.outGene)));
    
    if length(matches) > 1
        countDoubles = [countDoubles; matches];
    end
    %Check for scores 
    if ~isempty(matches)
        geneInts.functionalPC(matches) = ...
            interactionSummary.totalFunctional(i);
        geneInts.physicalPC(matches) = ...
            interactionSummary.totalPhysical(i);
        
        annotations = unique([geneInts.PubMedIDs{matches},...
            interactionSummary.PubMedIDs{i}]);
        geneInts.PubMedIDs(matches) = {annotations};

        if interactionSummary.totalPhysical(i) > 0 
            reverseMatches = ...
                intersect(find(strcmp(outGene, geneInts.inGene)), ...
                find(strcmp(inGene, geneInts.outGene))); 
            geneInts.physicalPC(reverseMatches) = ...
                interactionSummary.totalPhysical(i);            
        end 
    end
    
end

bothIndex = intersect(find(geneInts.functionalPC > 0), ...
    find(geneInts.physicalPC > 0));
geneInts.bothPC(bothIndex) = 1;
geneInts.totalScore = geneInts.functionalPC + geneInts.physicalPC;

%%
ruleNames = geneInts.name;
%This for loop puts the physical and fucntional annotations in the merged
%Results table 
for i = 1:length(mergedResults.name)
    thisRule = mergedResults.name(i);
    
    present = sum(strcmp(thisRule, ruleNames));
    
    if present == 0
        mergedResults.colorKey(i) = -1;
    end
    
    interactions = find(strcmp(thisRule{1}, ruleNames));
    
    intScores = geneInts(interactions,:);
    
    functionalScore = sum(intScores.functionalPC);
    physicalScore = sum(intScores.physicalPC);
    bothScore = (physicalScore > 0) & (functionalScore > 0);
    mergedResults.bothPC(i) = bothScore;
    mergedResults.physicalPC(i) = physicalScore;
    mergedResults.functionalPC(i) = functionalScore;
    mergedResults.totalScore(i) = physicalScore + functionalScore;
   
    if size(intScores,1) > 0
        annotations = unique([intScores.PubMedIDs{:}]);
    
        mergedResults.pubMedIDs(i) = {annotations};
        mergedResults.pubMedNum(i) = length(annotations);
    end
end

summaryNums = zeros(1, 4);


functionalCount = sum(mergedResults.functionalPC > 0 &...
    mergedResults.bothPC == 0);
physicalCount = sum(mergedResults.physicalPC > 0 &...
    mergedResults.bothPC == 0);
bothCount = sum(mergedResults.bothPC > 0);

rxnCount = size(networkInts,1);

summaryNums(1:3) = [functionalCount, physicalCount, ...
    bothCount]./rxnCount;
summaryNums(4) = 1 - sum(summaryNums);

end


%(function for step 4)
function [geneInts_merged, edgeInts_merged] = markNetworkInts(...
    edgeInts_merged, geneInts_merged, ...
    networkGeneInts)

    geneInts_merged.inNetwork = zeros(size(geneInts_merged,1),1);
    edgeInts_merged.inNetwork = zeros(size(edgeInts_merged,1),1);
    
    common = cellfun(@(x,y) intersect(find(strcmp(x,geneInts_merged.inGene)),...
        find(strcmp(y, geneInts_merged.outGene))), networkGeneInts.inGene,...
        networkGeneInts.outGene, 'Un', 0);
    
    geneInts_merged.inNetwork(cell2mat(common)) = 1;
    rules = unique(geneInts_merged.name(cell2mat(common)));
    
    commonRuleIndex = ...
        cellfun(@(x) find(strcmp(x,edgeInts_merged.name)), ...
        rules, 'Un', 0);
    
    edgeInts_merged.inNetwork(cell2mat(commonRuleIndex)) = 1;
    
end


%(function for step 6)
function [geneInts_merged, edgeInts_merged] = checkSignor(...
    edgeInts_merged, geneInts_merged)

    database = readtable('Signor/SignorDatabases2.xlsx');

    intsKey = readtable('Signor/InteractionTypeKey.xlsx');

    edgeInts_merged.signorDirected = zeros(size(edgeInts_merged,1),1);
    edgeInts_merged.signorUndirected = zeros(size(edgeInts_merged,1),1);
    edgeInts_merged.signorDIRECT = zeros(size(edgeInts_merged,1),1);
    edgeInts_merged.signorPubMed = cell(size(edgeInts_merged,1),1);
    edgeInts_merged.signorNumPMID = zeros(size(edgeInts_merged,1),1);
    edgeInts_merged.signorSent = cell(size(edgeInts_merged,1),1);
    edgeInts_merged.matchPMID = zeros(size(edgeInts_merged,1),1);
    
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
        
        edgeInts_merged.signorNumPMID(iEdge) = ...
            length(edgeInts_merged.signorPubMed{iEdge});
        
        
                
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

