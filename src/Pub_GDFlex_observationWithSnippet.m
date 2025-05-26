function Pub_GDFlex_observationWithSnippet(dataIdStr,ClippingRange,NNpair,Xinterval,Name,ZnormOn,LineWidth_GT)
  
 %{
 Name = 'Observation_scale02'; 
 dataIdStr = '146';
 ClippingRange = [27601,28100];
 NNpair = [27895,27932,2325,2362;...
          5155,5179,212,236];
 Xinterval = [50]
 ZnormOn = true;
 LineWidth_GT = [1.5,1.5];
 Pub_GDFlex_observationWithSnippet(dataIdStr,ClippingRange,NNpair,Xinterval,Name,ZnormOn,LineWidth_GT);
  % arg=(5155-5179) NN=(212-236)

 Name = 'Observation_localizedMismatch'; 
 dataIdStr = '247';
 ClippingRange = [120001,123000];
 NNpair = [121735,122113,15648,16026;...
           145658,145988,2836,3166];
 Xinterval = [200]
 ZnormOn = true;
 LineWidth_GT = [2.5,2];
 Pub_GDFlex_observationWithSnippet(dataIdStr,ClippingRange,NNpair,Xinterval,Name,ZnormOn,LineWidth_GT);

 Name = 'Observation_PV3noise'; 
 dataIdStr = '102';
 ClippingRange = [18501,20500];
 NNpair = [19197,19385,9511,9699;...
           21417,21569,1684,1836];
 Xinterval = [50]
 ZnormOn = true;
 LineWidthGT = [1.5,2];
 Pub_GDFlex_observationWithSnippet(dataIdStr,ClippingRange,NNpair,Xinterval,Name,ZnormOn,LineWidth_GT);

 Name = 'Observation_Degeneration_old'; 
 dataIdStr = '057';
 ClippingRange = [15001,17000];
 NNpair = [15977,16164,88,275;...
           36817,36931,4453,4567];
 Xinterval = [20]
 ZnormOn = true;
 LineWidth_GT = [2.5,2];
 Pub_GDFlex_observationWithSnippet(dataIdStr,ClippingRange,NNpair,Xinterval,Name,ZnormOn,LineWidth_GT);

 Name = 'Observation_scale34';
 dataIdStr = '206';
 ClippingRange = [28501,29500];
 NNpair = [29081,29151,7604,7674;...
          29081,29111,6560,6590;...
          29111,29150,23313,23352]; % 29111,29152,23313,23354];
 Xinterval = [20,10];
 ZnormOn = true;
 LineWidth_GT = [1.5,1.5];
 Pub_GDFlex_observationWithSnippet(dataIdStr,ClippingRange,NNpair,Xinterval,Name,ZnormOn,LineWidth_GT);

 Name = 'Observation_scale02_old'; 
 dataIdStr = '038';
 ClippingRange = [27501,28200];
 NNpair = [27895,27932,2325,2362;...
           19150,19176,3074,3100];
 Xinterval = [20]
 ZnormOn = true;
 LineWidth_GT = [1.5,1.5];
 Pub_GDFlex_observationWithSnippet(dataIdStr,ClippingRange,NNpair,Xinterval,Name,ZnormOn,LineWidth_GT);

 Name = 'Observation_scale01';
 dataIdStr = '191';
 ClippingRange = [140001,146000];
 NNpair = [143270,143572,34147,34449;...
          129843,130203,21861,22221;...
          143389,143572,28171,28354];
 Xinterval = [150,150];
 ZnormOn = true;
 LineWidth_GT = [1.5,1.5];
 Pub_GDFlex_observationWithSnippet(dataIdStr,ClippingRange,NNpair,Xinterval,Name,ZnormOn,LineWidth_GT);
 %}

 close all;

 %% Parameter
 LineWidth = 1.25;
 if exist('LineWidth_GT','var')
  LineWidth_GT_CLIP = LineWidth_GT(1);   
  LineWidth_GT_Pattern = LineWidth_GT(2);
 else
  LineWidth_GT_CLIP = 1.5;   
  LineWidth_GT_Pattern = 1.5;
 end
 %% Directory
 TSAD2021_Dir = '../data';

 %%
 RawDataQuery = [TSAD2021_Dir,'/',dataIdStr,'*.mat']; 
 RawDataList = dir(RawDataQuery); 
 RawFile =  [TSAD2021_Dir,'/',RawDataList(1).name];
 RawData = load(RawFile);
 DataName = RawData.name;
 Window_length = RawData.window_length;
 D = RawData.data;
 Anomaly_region =  RawData.anomaly_region;

 TargetNum = length(D) - Window_length + 1;
 fprintf('Data length=%d Clipped length=%d \n',length(D),ClippingRange(2)-ClippingRange(1)+1);
 fprintf('Discord window-length=%d trgNum=%d \n',Window_length,TargetNum);
 %%
 fig_id = figure(); Row= 2; Column = 1; 
 % t = tiledlayout(Row, Column, 'TileSpacing', 'compact', 'Padding', 'none');
 t = tiledlayout(Row, Column, 'TileSpacing', 'none', 'Padding', 'none');
 nexttile;
    ClippingTimes = (ClippingRange(1):ClippingRange(end));
    plot(ClippingTimes,D(ClippingTimes),'k'); xlim(ClippingRange);
    hold on;
    AnomalyTimes = (Anomaly_region(1):Anomaly_region(end));
    plot(AnomalyTimes,D(AnomalyTimes),'b','LineWidth',LineWidth_GT_CLIP);
    axis off; 

 nexttile;
 start = 1;
 for i = 1:size(NNpair,1)  
     from_i = NNpair(i,1);
     to_i = NNpair(i,2);
     length_i = to_i - from_i + 1;
     if ZnormOn
       D_trg_plot_i = normalize(D(from_i:to_i));
     else
       D_trg_plot_i = D(from_i:to_i);
     end
     [OVERLAP_i, overlap_length_i,overlap_range_i] = overlapping_range([from_i, to_i],Anomaly_region);
  
     fromNN_i = NNpair(i,3);
     toNN_i = NNpair(i,4);
     if ZnormOn
       D_NN_plot_i = normalize(D(fromNN_i:toNN_i));
     else
       D_NN_plot_i = D(fromNN_i:toNN_i); 
     end

    RangeX_i = (start:start+length_i-1);
    plot(RangeX_i,D_trg_plot_i,'k','LineWidth',LineWidth); axis off; hold on;
    if ~isempty(overlap_range_i)
      GroundT_i = (overlap_range_i(1):overlap_range_i(2)); 
      ground_start_i = overlap_range_i(1) - from_i + 1;
      ground_end_i = ground_start_i + length(GroundT_i) - 1;
      GroundT_relative_i = (ground_start_i:ground_end_i);
      plot(GroundT_relative_i+start-1,D_trg_plot_i(GroundT_relative_i),'b','LineWidth',LineWidth_GT_Pattern); axis off;
    end

    plot(RangeX_i,D_NN_plot_i,'r','LineWidth',LineWidth); axis off;

    if i < size(NNpair,1) 
     start = start + length_i - 1 + Xinterval(i); 
    end
 end

 %{
 FigFileStem = [dataIdStr,'_',DataName];
 FigFile = [OutputDir,'/',FigFileStem,'.fig'];
 saveas(fig_id,FigFile);
 
 FigFile = [OutputDir,'/',FigFileStem,'.svg'];
 saveas(fig_id,FigFile);
 %}
end

function [OVERLAP, overlap_length,overlap_range] = overlapping_range(range_1,range_2)
  from_1 = range_1(1);
  to_1 = range_1(2);
  from_2 = range_2(1);
  to_2 = range_2(2);

  overlap_range = [];
  if (from_2 <= from_1 && from_1 <= to_2)     % range_2 proceeds
    OVERLAP = true;  
    overlap_range = [from_1, min(to_1,to_2)];
  elseif (from_1 <= from_2 && from_2 <= to_1) % range_1 proceeds
    OVERLAP = true;  
    overlap_range = [from_2,min(to_1,to_2)];
  elseif (from_2 <= from_1 && to_1 <= to_2)      % range_2 includes range_1
    OVERLAP = true;  
    overlap_range = [from_1,to_1];
  elseif (from_1 <= from_2 && to_2 <= to_1)      % range_1 includes range_2
    OVERLAP = true;  
    overlap_range = [from_2,to_2];
  else
    OVERLAP = false;
  end

  if isempty(overlap_range)
    overlap_length = 0;
  else
    overlap_length = overlap_range(2) - overlap_range(1) + 1;
  end
end