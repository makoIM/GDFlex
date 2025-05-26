function Pub_GDFlex_observation_Row2(dataIdStrList,NNpair,Name)
  
 %{
 Name = 'localizedMismatch';
 dataIdStrList = {'053','078','057','232'};
 NNpair = [2765,2996,217,448;...
          110208,110552,70309,70653;...
          15977,16164,88,275;
          57626,57764,440,578];
 Pub_GDFlex_observation_Row2(dataIdStrList,NNpair,Name);

 Name = 'degeneratePattern_old';
 dataIdStrList = {'053','078','232','247'};
 NNpair = [2765,2996,217,448;...
          110208,110552,70309,70653;...
          57626,57764,440,578;...
          121735,122113,15648,16026];
 Pub_GDFlex_observation_Row2(dataIdStrList,NNpair,Name);
 
 Name = 'Observation_znorm';
 dataIdStrList = {'077','145','145','116'};
 NNpair = [110799,110803,54649,54653;...
          17226,17243,2057,2074;...
          17201,17233,4180,4212;...
          5550,5561,823,834];
 Pub_GDFlex_observation_Row2(dataIdStrList,NNpair,Name);
 %}

 close all;

 %% Parameter
 LineWidth = 1.25;
 LineWidth_GT = 1.5;

 %% Directory
 TSAD2021_Dir = '../data';
 
 %%
 dataIdStr = [];
 for i = 1:length(dataIdStrList)
  if i == 1
    dataIdStr = dataIdStrList{i};
  else
    dataIdStr = [dataIdStr,'_',dataIdStrList{i}];  
  end
 end

 %%
 fig_id = figure(); start = 1;
 Row = 2; Column = size(NNpair,1);
 t = tiledlayout(Row, Column, 'TileSpacing', 'compact', 'Padding', 'compact');
 for i = 1:size(NNpair,1)  
   dataIdStr_i = dataIdStrList{i};  
   RawDataQuery_i = [TSAD2021_Dir,'/',dataIdStr_i,'*.mat']; 
   RawDataList_i = dir(RawDataQuery_i); 
   RawFile_i =  [TSAD2021_Dir,'/',RawDataList_i.name];
   RawData_i = load(RawFile_i);
   DataName_i = RawData_i.name;
   D_i = RawData_i.data;
   Anomaly_region_i =  RawData_i.anomaly_region;

     from_i = NNpair(i,1);
     to_i = NNpair(i,2);
     length_i = to_i - from_i + 1;
     D_trg_norm_i = normalize(D_i(from_i:to_i));
     D_trg_raw_i = D_i(from_i:to_i);
     
     [OVERLAP_i, overlap_length_i,overlap_range_i] = overlapping_range([from_i, to_i],Anomaly_region_i);
  
     fromNN_i = NNpair(i,3);
     toNN_i = NNpair(i,4);
     D_NN_norm_i = normalize(D_i(fromNN_i:toNN_i));
     D_NN_raw_i = D_i(fromNN_i:toNN_i); 
   
   %% norm
   % subplot(Row,Column,i);
    nexttile(i);
    RangeX_i = (start:start+length_i-1)';
    plot(RangeX_i,D_trg_norm_i,'b','LineWidth',LineWidth); axis off; hold on;
    if ~isempty(overlap_range_i)
      GroundT_i = (overlap_range_i(1):overlap_range_i(2)); 
      ground_start_i = overlap_range_i(1) - from_i + 1;
      ground_end_i = ground_start_i + length(GroundT_i) - 1;
      GroundT_relative_i = (ground_start_i:ground_end_i);
      plot(GroundT_relative_i+start-1,D_trg_norm_i(GroundT_relative_i),'b','LineWidth',LineWidth_GT); axis off;
    end
    plot(RangeX_i,D_NN_norm_i,'r','LineWidth',LineWidth); axis off;

   %% raw
   % subplot(Row,Column,i+Column);
   nexttile(i+Column);
    RangeX_i = (start:start+length_i-1);
    plot(RangeX_i,D_trg_raw_i,'b','LineWidth',LineWidth); axis off; hold on;
    if ~isempty(overlap_range_i)
      GroundT_i = (overlap_range_i(1):overlap_range_i(2)); 
      ground_start_i = overlap_range_i(1) - from_i + 1;
      ground_end_i = ground_start_i + length(GroundT_i) - 1;
      GroundT_relative_i = (ground_start_i:ground_end_i);
      plot(GroundT_relative_i+start-1,D_trg_raw_i(GroundT_relative_i),'b','LineWidth',LineWidth_GT); axis off;
    end
    plot(RangeX_i,D_NN_raw_i,'r','LineWidth',LineWidth); axis off;



    %{
    if i < size(NNpair,1) 
     start = start + length_i - 1 + Xinterval(i); 
    end
    %}
 end

 %{
 FigFile = [OutputDir,'/',dataIdStr,'.fig'];
 saveas(fig_id,FigFile);
 FigFile = [OutputDir,'/',dataIdStr,'.svg'];
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