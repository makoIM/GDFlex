function Pub_GDFlex_segStatistics_chap4_91(InfoFile, ClippingRange,AlignOption)
  %{
 ClippingRange = [26001,30000];
 AlignOption = "left";% "left" "center"
 InfoFile = '../sample/InfoDemo_91.mat';  % InfoDemo_91.mat is generated as an output of Pub_GDFlex_execute(91, "interNoise").
 Pub_GDFlex_segStatistics_chap4_91(InfoFile, ClippingRange, AlignOption);
  %}

  close all;
  % data
  Name = '091';
  TSAD2021_Dir = '../data';
  RawDataQuery = [TSAD2021_Dir,'/',Name,'*.mat']; 
  RawDataList = dir(RawDataQuery); 
  RawFile =  [TSAD2021_Dir,'/',RawDataList.name];
  RawData = load(RawFile);
  DataName = RawData.name;
  Training_region = RawData.training_region;
  Anomaly_region =  RawData.anomaly_region;
  
  MatLoad = load(InfoFile);
  Info = MatLoad.InfoPosNeg;

  % magnitude distribution
  Mag = Info.MagInfo.magnitude;
  D = Info.MagInfo.data_org;
  MatThres_Pos = Info.Param.SegStatistics.positive2.MaT;
  MatThres_Neg = Info.Param.SegStatistics.negative2.MaT;
  
  % fprintf('%s neg=%f pos=%f\n',dataIdStr,MatThres_Neg,MatThres_Pos);
  %% target segment
  Names = Info.SegInfo.SegDB_names;
  SegTrg = Info.SegInfo.SegDB;
  SegTrain = Info.SegInfo_train.SegDB;
 
  %% total plot 
  fig_id_1 = figure();
  Row = 2; Column = 1;
  t = tiledlayout(Row, Column, 'TileSpacing', 'none', 'Padding', 'none');
  nexttile; % subplot(Row,Column,1);
  TrainRegion = (1:Training_region);
  TrgRegion =  (Training_region+1:length(D));
  fprintf('train=%d length=%d\n',Training_region,length(D));
  plot(TrainRegion,D(TrainRegion),'r'); hold on;
  plot(TrgRegion,D(TrgRegion),'k');
  xline(Training_region,'k--');
  xlim([1, length(D)]);
    ax = gca;         
    ax.FontSize = 14;  
  axis off;

  % clipping plot 
  nexttile; % subplot(Row,Column,2);
  ClippingTimes = [ClippingRange(1):ClippingRange(end)];
  plot(ClippingTimes,D(ClippingTimes),'r'); xlim(ClippingRange);
  hold on;
  AnomalyTimes = [Anomaly_region(1):Anomaly_region(end)];
  plot(AnomalyTimes,D(AnomalyTimes),'b','LineWidth',2);
  ax = gca;         
  ax.FontSize = 14;  
  axis off;

  %{
  OutputFile = [OutputDir,'/','91_data']; 
  saveas(fig_id_1,[OutputFile,'.fig']);
  saveas(fig_id_1,[OutputFile,'.svg']);
  %}
  %%  magnitude
  fig_id_2 = figure();
  Row = 3; Column = 2;
  subplot(Row,Column,1);
    h1 = histogram(Mag); xlim([-3000,3000]); 
    xline(MatThres_Pos,'r--'); xline(MatThres_Neg,'r--');
    ax = gca;         
    ax.FontSize = 14; 
  subplot(Row,Column,2);
    Mag_signal = Mag(Mag >MatThres_Pos | Mag < MatThres_Neg); 
    edges1 = (MatThres_Neg+250:-500:-5000);
    edges1 = fliplr(edges1);
    edges2 = (MatThres_Neg+250:500:MatThres_Pos-500);
    edges3 = (MatThres_Pos:500:5000);
    edges = [edges1,edges2, edges3];
    h2 = histogram(Mag_signal, 'BinEdges', edges);
    xline(MatThres_Pos+500,'r--'); xline(MatThres_Neg+220,'r--');
    ax = gca;         
    ax.FontSize = 14; 
  % training
  SegTrain_restriction = SegTrain(SegTrain(:,Names=="length")>= 200,:);

  AnomalyRegion = [0,0];
  subplot(Row,Column,3);
    
    plot_motifCandidatesWithAnomalyRegion(SegTrain_restriction,Names,2,AlignOption,D,AnomalyRegion);
    ax = gca;         
    ax.FontSize = 14; 
    xlim([0 560]);
  subplot(Row,Column,4);
    plot_motifCandidatesWithAnomalyRegion(SegTrain_restriction,Names,-2,AlignOption,D,AnomalyRegion)
    ax = gca;         
    ax.FontSize = 14;  
    xlim([0 560]);

  % target 
  SegTrg_restriction = SegTrg(SegTrg(:,Names=="length")>= 200,:);
  subplot(Row,Column,5);
    plot_motifCandidatesWithAnomalyRegion(SegTrg_restriction,Names,2,AlignOption,D,AnomalyRegion);
    ax = gca;    
    ax.FontSize = 14; 
    xlim([0 560]);
  subplot(Row,Column,6);
    plot_motifCandidatesWithAnomalyRegion(SegTrg_restriction,Names,-2,AlignOption,D,AnomalyRegion); 
    Index_anomaly = SegTrg_restriction(:,Names=="length") >= 400; 
    from = SegTrg_restriction(Index_anomaly,Names=="from");
    to = SegTrg_restriction(Index_anomaly,Names=="to");
    plot((1:to(1)-from(1)+1),D(from(1):to(1)),'b','LineWidth',2);
    ax = gca;         
    ax.FontSize = 14; 
    xlim([0 560]);

  %{
  OutputFile = [OutputDir,'/','91_stats'];
  saveas(fig_id_2,[OutputFile,'.fig']);
  saveas(fig_id_2,[OutputFile,'.svg']);  
  %}
end

% AlignOption = "center"; "left"
function plot_motifCandidatesWithAnomalyRegion(SegDB,Names,Type,AlignOption,D,AnomalyRegion) 
    SegDB = SegDB(SegDB(:,Names == "type")==Type,:);
    
    %% plot
    Overlap_id = []; Overlap_rate = [];
    if AlignOption == "center"
      min_t = inf; max_t = -inf; min_y = 0; max_y = 0;  
    elseif AlignOption == "left"
      min_t = 0;   max_t = 0; min_y = inf; max_y = -inf;
    end

    for i = 1:size(SegDB,1)
             id_i     = SegDB(i,Names=="id");
             center_i = SegDB(i,Names=="center");
             from_i   = SegDB(i,Names=="from"); 
             to_i     =   SegDB(i,Names=="to"); 
             mag_i    =  SegDB(i,Names=="magnitude");
             % if center_i == 19386;fprintf('debug\n');end
             if exist('AnomalyRegion','var')
               [OVERLAP,overlap_length] = overlapping_range(AnomalyRegion,[from_i,to_i]); 
               if OVERLAP
                 Overlap_id = [Overlap_id,i];
                 motif_length_i = to_i-from_i+1; 
                 anomaly_length = AnomalyRegion(2)-AnomalyRegion(1)+1;
                 Overlap_rate = [Overlap_rate,max(overlap_length/motif_length_i,overlap_length/anomaly_length)];
               end
             end
             data_oi = D(from_i:to_i); 
             if AlignOption == "center"
               adjustment_oi = adjustment_plot(data_oi,mag_i);
               plot(from_i-center_i:to_i-center_i,data_oi-adjustment_oi);
             elseif AlignOption == "left"
               plot(1:to_i-from_i+1,data_oi);  
             end

             if i == 1; hold on; end;
             if AlignOption == "center"
               min_t = min(from_i-center_i,min_t); max_t = max(to_i-center_i,max_t);  
               max_y = max(max_y,max(data_oi-adjustment_oi)); min_y = min(min_y,min(data_oi-adjustment_oi));
             else
               max_t = max(max_t,to_i-from_i+1);
               max_y = max(max_y,max(data_oi)); min_y = min(min_y,min(data_oi));
             end
    end
           if max_t ~= -inf && min_t ~= inf
             % LIM_X = lim_margin([min_t max_t],0.1); xlim(LIM_X);
             LIM_X = lim_margin([min_t max_t],0.01); xlim(LIM_X);
             LIM_Y = lim_margin([min_y max_y],0.1); ylim(LIM_Y);
           end
     if ~isempty(Overlap_id) && i == size(SegDB,1)
             for oi = 1:length(Overlap_id) 
               overlap_id_oi = Overlap_id(oi);  
               center_i = SegDB(overlap_id_oi, Names=="center"); 
               from_i = SegDB(overlap_id_oi,Names=="from"); 
               to_i = SegDB(overlap_id_oi,Names=="to"); 
               mag_i = SegDB(overlap_id_oi,Names=="magnitude"); 
               if Overlap_rate(oi) >= 0.5
                 Color = 'k';
               else
                 Color = 'b';  
               end
               data_oi = D(from_i:to_i);
               if AlignOption == "center"
                 adjustment_oi = adjustment_plot(data_oi,mag_i);
                 plot(from_i-center_i:to_i-center_i,data_oi-adjustment_oi,Color,'LineWidth',2);
               elseif AlignOption == "left"
                 plot(1:to_i-from_i+1,data_oi,Color,'LineWidth',2);  
               end
             end
     end
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

function adjustment_oi = adjustment_plot(data_oi,mag_i)
   if mag_i >= 0
     adjustment_oi = max(data_oi);
   elseif mag_i < 0
     adjustment_oi = min(data_oi);    
   end
end


function LIM = lim_margin(D,Rate)
  if ~exist('Rate','var')
    Rate = 0.2;
  end

  Margin = Rate*(max(D) - min(D));
  if Margin ~= 0
    LIM = [min(D)-Margin, max(D)+Margin];
  else
    LIM = [];
  end
end