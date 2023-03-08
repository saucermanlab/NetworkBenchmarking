%%Just generating a bar graph for the workflow picture 
clear all
close all

%Bar graph with four outputs

allNums = [ .65 .75 .85 .45];
directedNums = [.8 .5 .65 .55];
indices = [0.5:1:length(directedNums)];

%all Nums
figure;
subplot(1,2,1)

hold on 
b = barh(indices, allNums,...
    'EdgeColor', 'black', 'FaceColor', 'flat');
hold on

ylim([indices(1) - 1, indices(end) + 1]);
xlim([0 1]);

set(gca,'xticklabel',{[]}) 
set(gca,'yticklabel',{[]}) 
ax1 = gca;
ax1.Visible = 'off';

b.CData(4,:) = [0 0 0];
b.CData(3,:) = [.25 .25 .25];
b.CData(2,:) = [.5 .5 .5];
b.CData(1,:) = [.75 .75 .75];

%directed Nums
subplot(1,2,2)

c = barh(indices, directedNums,...
    'EdgeColor', 'black', 'FaceColor', 'flat');
hold on

ylim([indices(1) - 1, indices(end) + 1]);
xlim([0 1]);

set(gca,'xticklabel',{[]}) 
set(gca,'yticklabel',{[]}) 
ax2 = gca;
ax2.Visible = 'off';

c.CData(4,:) = [0 0 0];
c.CData(3,:) = [.25 .25 .25];
c.CData(2,:) = [.5 .5 .5];
c.CData(1,:) = [.75 .75 .75];

hold off
