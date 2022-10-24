clear all; clc; close all; 
%Loop through different models (with edge additions)
%Note that you may need to edit the existing hypertrophy 
Modelfile_base ='Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_0.xlsx'; % Change the path to the file current path %RyallHypertrophy_0.xlsx
newEdge_path = 'Z:\Taylor\Code\Benchmarking_Extended\hypertrophyBenchmarking.xlsx';
[~, txt, ~] = xlsread(newEdge_path,'new ints table_TE');
newEdge = txt(2:end,1);
%make subset of newEdge
new(1) = newEdge(1);
% % % new(2) = newEdge(6);  
new(2) = newEdge(7);
% % % new(4) = newEdge(8);
% % % new(5) = newEdge(10);
% % % newEdge = new';

%make copies of base file
inputFolder = 'Z:\Taylor\Code\Benchmarking_Extended';
sourceFile = fullfile(inputFolder, 'RyallHypertrophy_0.xlsx');  
outputFolder = inputFolder; 
% Now make n copies, with different names, in the output folder.
for k = 1:length(newEdge)
    baseFileName = sprintf('RyallHypertrophy_%d.xlsx', k);
    outputFile = fullfile(outputFolder,baseFileName);
    copyfile(sourceFile, outputFile);
end
% write new edges to appropriate files
for i = 1:length(newEdge)
    xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(i),'.xlsx'),[{'middle'},{'r195'},newEdge(i),1,1.223,0.5],'reactions','A197'); 
end
%%%%%deleteme%%%%%
count = length(newEdge);
% % % for i = 1:length(newEdge)-1
% % %     for j = (i+1):length(newEdge)
% % %         count = count + 1;
% % %         xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(count),'.xlsx'),[{'middle'},{'r195'},newEdge(i),1,1.223,0.5],'reactions','A197'); 
% % %         xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(count),'.xlsx'),[{'middle'},{'r196'},newEdge(j),1,1.223,0.5],'reactions','A198'); 
% % %     end
% % % end
% % % for i = 1:length(newEdge)-2
% % %     for j = (i+1):length(newEdge)-1
% % %         for k = (j+1):length(newEdge)
% % %             count = count + 1;
% % %             xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(count),'.xlsx'),[{'middle'},{'r193'},newEdge(i),1,1.223,0.5],'reactions','A195'); 
% % %             xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(count),'.xlsx'),[{'middle'},{'r194'},newEdge(j),1,1.223,0.5],'reactions','A196'); 
% % %             xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(count),'.xlsx'),[{'middle'},{'r195'},newEdge(k),1,1.223,0.5],'reactions','A197'); 
% % %         end
% % %     end
% % % end
% % % for i = 1:length(newEdge)-3
% % %     for j = (i+1):length(newEdge)-2
% % %         for k = (j+1):length(newEdge)-1
% % %             for l = (k+1):length(newEdge)
% % %                 count = count + 1;
% % %                 xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(count),'.xlsx'),[{'middle'},{'r193'},newEdge(i),1,1.223,0.5],'reactions','A195'); 
% % %                 xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(count),'.xlsx'),[{'middle'},{'r194'},newEdge(j),1,1.223,0.5],'reactions','A196'); 
% % %                 xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(count),'.xlsx'),[{'middle'},{'r195'},newEdge(k),1,1.223,0.5],'reactions','A197'); 
% % %                 xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(count),'.xlsx'),[{'middle'},{'r196'},newEdge(l),1,1.223,0.5],'reactions','A198');
% % %             end
% % %         end
% % %     end
% % % end
% % % count = count+1;
% % % for i = count
% % %     xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(i),'.xlsx'),[{'middle'},{'r194'},newEdge(1),1,1.223,0.5],'reactions','A196'); %r194 & A196 for RyallHypertrophy
% % %     xlswrite(strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(i),'.xlsx'),[{'middle'},{'r195'},newEdge(2),1,1.223,0.5],'reactions','A197'); %r195 & A197 for RyallHypertrophy
% % % end
%%%%%deleteme%%%%%

%% Validate network files
for i = 0:count%length(newEdge) 
    Modelfile_path =strcat('Z:\Taylor\Code\Benchmarking_Extended\RyallHypertrophy_',num2str(i),'.xlsx'); % Change the path to the file current path
    Validationfile_path ='Z:\Taylor\Code\Benchmarking_Extended\Gold_Standard_Validation.xlsx';% Change the path to the file current path
    Int_time = 40;
    Steady_time = 40;
    Threshold = 0.1; % percent of changes
    Model_version = 1; % 1= Original 2= Modified
    [percentMatch, resultChart, BMatch, byClass] = Automated_Validation_V1(Modelfile_path, Validationfile_path, Int_time, Steady_time,Threshold, Model_version);
    percentMatch
    
    %Rename Validation results file
    movefile('Validation_Results.xlsx',strcat('Validation_Results_',num2str(i),'.xlsx'));
    
    Name(i+1) = cellstr(strcat('Hypertrophy_',num2str(i)));
% % %     if i == 0
% % %         Name(i+1) = {'Base Network'};
% % %     else
% % %         Name(i+1) = newEdge(i);
% % %     end
    Percent(i+1) = percentMatch;
end
Name(1)={'Base Network'};Name(2)={'PKC=>PKC'};Name(3)={'Akt=>Raf1'};Name(4)={'PKC=>PKC; Akt=>Raf1'};
EdgeValidations = table(Name',Percent'); 

%% Graph differences in Validation outcomes
P = zeros(450,count+1);
for i=0:count %length(newEdge)
    [num, txt, raw] = xlsread(strcat('Validation_Results_',num2str(i),'.xlsx'));
    mtch = txt(2:end,8);
    for j=1:length(mtch)
        mtch{j}=strcmpi('yes',mtch{j});
    end
    mtch = cell2mat(mtch);
    P(:,i+1) = mtch;
end
selection = [];
for i=1:length(P)
    if sum(P(i,:))==0
        selection = [selection, i];
    elseif sum(P(i,:))==count+1
        selection = [selection, i]; 
% % %     elseif P(i,1)==1 
% % %         selection = [selection, i];   
    end
end
P(selection,:)=[];
lbl = [" "];
for i=2:length(txt)
    lbl = vertcat(lbl,convertCharsToStrings(strjoin(txt(i,2:4),' ')));
end
lbl(selection,:)=[];

figure;
set(gca, 'Visible', 'on');
% % % cmaprange = [0.5:0.005:1];
% % % blank = [zeros(1,length(cmaprange) - 20), 0.01:0.05:1];
% % % myrgbcmap = [blank',blank',cmaprange';1 1 1; flipud(cmaprange'),flipud(blank'),flipud(blank')];
colormap(gray);
imagesc(P,[-1,1]);
set(gca,'XAxisLocation','bottom');
set(gca,'XTick',1:length(Name));
set(gca,'XTickLabel',Name,'fontsize',10);
xlabel('Edge Added','fontsize',20);
set(gca,'YTick',1:length(P));
set(gca,'YTickLabel',lbl,'fontsize',10);
ylabel('Experimental Validations','fontsize',16);
% % % title([strcat('Change in ',{' '},specs2{dpos})]);
xtickangle(90)
% % % hcb=colorbar;
% % % set(get(hcb,'label'),'string','Predicted Outcome');
% % % set(get(hcb,'label'),'fontsize',16); set(get(hcb,'label'),'rotation',90);

newEdge2 = [{'Original Network'},newEdge'];
x = categorical(newEdge2);
x = reordercats(x,newEdge2);
y = EdgeValidations{:,2}';
figure
h1=axes;
bb=barh(x,y);
xl = xline(EdgeValidations.Var2(1),'--r');
bb.FaceColor = '#999999';
xlabel('Validation Percentage')
ylabel('Edge Addition')
xlim([30 80])
set(h1,'Ydir','reverse')