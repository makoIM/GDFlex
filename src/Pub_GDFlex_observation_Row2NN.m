function Pub_GDFlex_observation_Row2NN(dataIdStrList,RangeCell)
  
 %{ 
 dataIdStrList = {'065','067','210'...
                 ;'065','067','210'};  
 RangeCell =     {[8824,8842,308,326],[5572,5596,4709,4733],[70097,70201,34957,35061]...
                 ;[7004,7028,2036,2060],[7013, 7037,2436,2460],[74901, 74998, 4429, 4526]};  ;
 Pub_GDFlex_observation_Row2NN(dataIdStrList,RangeCell);

 %}

 close all;

 %% Parameter
 LineWidth = 1.5;
 LineWidth_GT = 1.25;

 %% Directory
 TSAD2021_Dir = '../data';

 %%
 dataIdStr = [];
 for i = 1:size(dataIdStrList,2)
  if i == 1
    dataIdStr = dataIdStrList{1,i};
  else
    dataIdStr = [dataIdStr,'_',dataIdStrList{1,i}];  
  end
 end

 %%
 fig_id = figure(); 
 Row = size(RangeCell,1); Column = size(RangeCell,2);
 t = tiledlayout(Row, Column, 'TileSpacing', 'compact', 'Padding', 'compact');
 for i = 1:Row
   for j = 1:Column
     nexttile((i-1)*Column + j);
     dataIdStr_ij = dataIdStrList{i,j};  
     RawDataQuery_ij = [TSAD2021_Dir,'/',dataIdStr_ij,'*.mat']; 
     RawDataList_ij = dir(RawDataQuery_ij); 
     RawFile_ij =  [TSAD2021_Dir,'/',RawDataList_ij.name];
     RawData_ij = load(RawFile_ij);
     DataName_ij = RawData_ij.name;
     D_i = RawData_ij.data;
     Anomaly_region_ij =  RawData_ij.anomaly_region;

     NNpair_ij = RangeCell{i,j}; 
     from_ij = NNpair_ij(1);
     to_ij = NNpair_ij(2);
     length_i = to_ij - from_ij + 1;
     D_trg_norm_i = normalize(D_i(from_ij:to_ij));
       
     fromNN_i = NNpair_ij(3);
     toNN_i = NNpair_ij(4);
     D_NN_norm_i = normalize(D_i(fromNN_i:toNN_i));

     Range_ij = (1:length_i)';
     plot(Range_ij,D_trg_norm_i,'k','LineWidth',LineWidth); hold on;
     plot(Range_ij,D_NN_norm_i,'r','LineWidth',LineWidth);

     [OVERLAP_ij, overlap_length_ij,overlap_range_ij] = overlapping_range([from_ij, to_ij],Anomaly_region_ij);

     if ~isempty(overlap_range_ij)
      GroundT_i = (overlap_range_ij(1):overlap_range_ij(2)); 
      ground_start_i = overlap_range_ij(1) - from_ij + 1;
      ground_end_i = ground_start_i + length(GroundT_i) - 1;
      GroundT_relative_i = (ground_start_i:ground_end_i);
      plot(GroundT_relative_i,D_trg_norm_i(GroundT_relative_i),'b','LineWidth',LineWidth_GT); 
     end
    
     axis off;
   end
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