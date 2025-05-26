function Pub_GDFlex_subsequences(DataIdStrCell,RangeMatrix,Xinterval,ZnormOn,AnomalyOn,SubplotOn)
 %{

 DirName = 'smallScale'; FileName = '3_7_83';
 DataIdStrCell = {'003','027','083'};
 RangeMatrix =   {[46568,46988],[4112,4248],[143279,143541]};  
 Xinterval = [30,30];
 ZnormOn = true;
 AnomalyOn = true;
 SubplotOn = true; 
 Pub_GDFlex_subsequences(DataIdStrCell,RangeMatrix,Xinterval,ZnormOn,AnomalyOn,SubplotOn);
 
 DirName = 'largeScale'; FileName = '206';
 DataIdStrCell = {'206','206','206'; '206','206','206'};
 RangeMatrix =   {[29081,29111],[29086,29115],[29081,29115];
                  [29081,29111],[29111,29152],[29081,29152]};  
 Xinterval = [10,30;30+40,70];
 ZnormOn = false;
 AnomalyOn = false;
 SubplotOn = false; 
 Pub_GDFlex_subsequences(DataIdStrCell,RangeMatrix,Xinterval,ZnormOn,AnomalyOn,SubplotOn);
 
 DirName = 'preparation_PVnum'; FileName = '206_43';
 DataIdStrCell = {'206','206','206','206'; '043','043','043','043'};
 RangeMatrix =   {[29081,29086],[29081,29111],[29081,29115],[29081,29152];[19260,19309],[19260,19334],[19260,19385],[19260,19404]};  % 19197
 Xinterval = [20+30,20+22,18+58;30,30,30];
 ZnormOn = false;
 AnomalyOn = false;
 SubplotOn = false; 
 Pub_GDFlex_subsequences(DataIdStrCell,RangeMatrix,Xinterval,ZnormOn,AnomalyOn,SubplotOn);
 
 RangeMatrix =   {[29081,29086],[29081,29111],[29081,29115],[29081,29152];[19260,19309],[19260,19334],[19260,19385],[19260,19404]};  % 19197
 R11 = RangeMatrix{1,1}; R12 = RangeMatrix{1,2}; R13 = RangeMatrix{1,3};
 R21 = RangeMatrix{2,1}; R22 = RangeMatrix{2,2}; R23 = RangeMatrix{2,3};
 (R21(2) - R21(1)) -(R11(2) - R11(1))   % 44
 (R22(2) - R22(1)) -(R12(2) - R12(1))   % 44
 (R23(2) - R23(1)) -(R13(2) - R13(1))   % 91
 %}
 close all;

 %% Parameter
 LineWidth = 1.5;
 %% Directory
 TSAD2021_Dir = '../data';

 %% argument
 if ~exist('AnomalyOn','var')
   AnomalyOn = false;
 end
 if ~exist('ZnormOn','var')
   ZnormOn = false;
 end
if ~exist('SubplotOn','var')
   SubplotOn = false;
 end

 %%
 fig_id = figure(); Row= size(DataIdStrCell,1); 
 if ~SubplotOn
  Column = 1;    
  t = tiledlayout(Row, Column, 'TileSpacing', 'compact', 'Padding', 'none');
  % t = tiledlayout(Row, Column, 'TileSpacing', 'none', 'Padding', 'none');
 else
  Column = size(DataIdStrCell,2);    
 end
 for i = 1:Row
   if ~SubplotOn
     nexttile;   
   end
   start = 1;
   for j = 1:length(DataIdStrCell(i,:))
     dataIdStr_ij = DataIdStrCell{i,j};
     [D_ij, Window_length_ij, Anomaly_region_ij,name_ij] = loadData(TSAD2021_Dir,dataIdStr_ij);
     Range_ij = RangeMatrix{i,j};
     ClippingLength_ij = Range_ij(2) - Range_ij(1) +1;
     ClipingTime_ij = (Range_ij(1):Range_ij(2));

     RangeX_ij = (start:start+ClippingLength_ij-1);
     if ZnormOn 
       D_ij_Clip = normalize(D_ij(ClipingTime_ij));  
     else
       D_ij_Clip = D_ij(ClipingTime_ij);
     end

     if SubplotOn
       subplot(Row,Column,(i-1)*Column+j);  
     end
     plot(RangeX_ij,D_ij_Clip,'k','LineWidth',LineWidth); hold on;
     
     if AnomalyOn
       startRel_anomaly_ij = max([Range_ij(1),Anomaly_region_ij(1)]) - Range_ij(1) + 1;
       start_anomaly_ij = start + startRel_anomaly_ij -1;
       length_Anomaly_ij = min([Anomaly_region_ij(2),Range_ij(2)]) - max([Range_ij(1),Anomaly_region_ij(1)]) + 1; 
       RangeXA_ij = (start_anomaly_ij: start_anomaly_ij+length_Anomaly_ij-1);    
       TimeAnomaly_ij = (startRel_anomaly_ij:startRel_anomaly_ij+length_Anomaly_ij-1);
       plot(RangeXA_ij,D_ij_Clip(TimeAnomaly_ij),'b','LineWidth',LineWidth); 
     end

     axis off; 
     if ~SubplotOn && j < length(DataIdStrCell(i,:)) 
       start = start + ClippingLength_ij - 1 + Xinterval(i,j); 
     % elseif j == length(DataIdStrCell(i,:))
     % axis off; 
     end
   end
 end

 %{
 FigFile = [OutputDir,'/',FileName,'.fig'];
 saveas(fig_id,FigFile);
 
 FigFile = [OutputDir,'/',FileName,'.svg'];
 saveas(fig_id,FigFile);
 %}
end

function [D, Window_length, Anomaly_region,DataName] = loadData(TSAD2021_Dir,dataIdStr_ij)
  RawDataQuery = [TSAD2021_Dir,'/',dataIdStr_ij,'*.mat'];
  RawDataList = dir(RawDataQuery); 
  RawFile =  [TSAD2021_Dir,'/',RawDataList(1).name];
  RawData = load(RawFile);
  DataName = RawData.name;
  Window_length = RawData.window_length;
  D = RawData.data;
  Anomaly_region =  RawData.anomaly_region;
end
