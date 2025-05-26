function Pub_GDFlex_segStatistics_chap1_91(InfoFile,DiscordPair,SpikeletPair)
  %{
  InfoFile = '../sample/InfoDemo_91.mat';  % InfoDemo_91 is obtained by Pub_GDFlex_execute(91, "interNoise") 
  DiscordPair = [123360,123609,12796,13045];
  SpikeletPair = [113872,114433,25973,26534];
  Pub_GDFlex_segStatistics_chap1_91(InfoFile,DiscordPair,SpikeletPair);
  %}

  close all;
  
  %%
  TSAD2021_Dir = '../data';
  dataIdStr = '091'; AlignOption = "left"; 

  %% InfoPosNeg
  MatLoad = load(InfoFile);
  Info = MatLoad.InfoPosNeg;

  % magnitude distribution
  Mag = Info.MagInfo.magnitude;
  D = Info.MagInfo.data_org;
  MatThres_Pos = Info.Param.SegStatistics.positive2.MaT;
  SuppMediam_Pos = Info.Param.SegStatistics.positive2.Supp_median;
  SuppStd_Pos = Info.Param.SegStatistics.positive2.Supp_std;
  MatThres_Neg = Info.Param.SegStatistics.negative2.MaT;
  SuppMediam_Neg = Info.Param.SegStatistics.negative2.Supp_median;
  SuppStd_Neg = Info.Param.SegStatistics.negative2.Supp_std;
  
  %% target segment
  Names = Info.SegInfo.SegDB_names;
  SegTrg = Info.SegInfo.SegDB;
  SegTrain = Info.SegInfo_train.SegDB;

  fprintf('Spikelet lenght_MediamStd(%d+-%d,%d+-%d) ThresMat(%d,%d) trgNum=%d\n'...
          ,round(SuppMediam_Pos),round(SuppStd_Pos)...
          ,round(SuppMediam_Neg),round(SuppStd_Neg)...
          ,round(MatThres_Neg),round(MatThres_Pos)...
          ,size(SegTrg,1));

  %% source data
  ClippingRange = [113601,115000];ClippingRange = [112001,116000];
  RawDataQuery = [TSAD2021_Dir,'/',dataIdStr,'*.mat']; 
  RawDataList = dir(RawDataQuery); 
  RawFile =  [TSAD2021_Dir,'/',RawDataList.name];
  RawData = load(RawFile);
  D = RawData.data;
  DataName = RawData.name;
  Window_length = RawData.window_length;
  Training_region = RawData.training_region;
  Anomaly_region =  RawData.anomaly_region;

  TargetNum = length(D) - Window_length + 1;
  fprintf('Data length=%d Clipped length=%d \n',length(D),ClippingRange(2)-ClippingRange(1)+1);
  fprintf('Discord window-length=%d trgNum=%d \n',Window_length,TargetNum);
  %%
  [Znorm_trg_discord,Znorm_NN_discord,times_trg_discord,times_NN_discord] = ...
    extractZnormData(D,DiscordPair);

  [Znorm_trg_spikelet,Znorm_NN_spikelet,times_trg_spikelet,times_NN_spikelet] = ...
    extractZnormData(D,SpikeletPair);

  Length = max([length(Znorm_trg_discord),length(Znorm_trg_spikelet)]);
  
  %% total plot 
  fig_id = figure();
  Row = 3; Column = 2;
  % t = tiledlayout(Row, Column, 'TileSpacing', 'compact', 'Padding', 'compact');
  t = tiledlayout(Row, Column, 'TileSpacing', 'compact', 'Padding', 'none');

  % subplot(Row,Column,[1,2]);
  nexttile([1 2]);
    ClippingTimes = [ClippingRange(1):ClippingRange(end)];
    plot(ClippingTimes,D(ClippingTimes),'k'); xlim(ClippingRange);
    hold on;
    AnomalyTimes = (Anomaly_region(1):Anomaly_region(end));
    plot(AnomalyTimes,D(AnomalyTimes),'b','LineWidth',2);
    axis off; 

  %% target
  AnomalyRegion = [0,0];

  % subplot(Row,Column,3);
  nexttile(5);
    SegTrg_restriction = SegTrg(SegTrg(:,Names=="from")>= ClippingRange(1) & SegTrg(:,Names=="to") <= ClippingRange(2),:);
    plot_motifCandidatesWithAnomalyRegion(SegTrg_restriction,Names,2,AlignOption,D,AnomalyRegion);
    xlim([1 Length]);
    axis off; 

  % subplot(Row,Column,4);
  nexttile(6);
    SegTrg_restriction = SegTrg(SegTrg(:,Names=="from")>= ClippingRange(1) & SegTrg(:,Names=="to") <= ClippingRange(2) & SegTrg(:,Names=="length")>= 100,:);
    plot_motifCandidatesWithAnomalyRegion(SegTrg_restriction,Names,-2,AlignOption,D,AnomalyRegion); 
    Index_anomaly = SegTrg_restriction(:,Names=="length") >= 400; 
    from = SegTrg_restriction(Index_anomaly,Names=="from");
    to = SegTrg_restriction(Index_anomaly,Names=="to");
    plot((1:to(1)-from(1)+1),D(from(1):to(1)),'k','LineWidth',2);
    xlim([1 Length]);
    axis off; 

  % subplot(Row,Column,5);
  nexttile(3);
    plot((1:length(times_trg_discord)),Znorm_NN_discord,'r','LineWidth',1.5);hold on;
    plot((1:length(times_trg_discord)),Znorm_trg_discord,'k','LineWidth',1.5); 
    xlim([1 Length]); axis off; 
  % subplot(Row,Column,6);
  nexttile(4);
    plot((1:length(times_trg_spikelet)),Znorm_NN_spikelet,'r','LineWidth',1.5);hold on;
    plot((1:length(times_trg_spikelet)),Znorm_trg_spikelet,'k','LineWidth',2); 
    xlim([1 Length]); axis off; 

    startAR = Anomaly_region(1) - times_trg_spikelet(1) + 1;
    endAR = startAR + (Anomaly_region(2) -  Anomaly_region(1) + 1) -1 ;
    plot((startAR:endAR),Znorm_trg_spikelet(startAR:endAR),'b','LineWidth',2);

  %{
  OutputFile = [OutputDir,'/','91_top1'];
  saveas(fig_id,[OutputFile,'.fig']);
  saveas(fig_id,[OutputFile,'.svg']);  
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

function [Znorm_trg_discord,Znorm_NN_discord,times_trg_discord,times_NN_discord] = extractZnormData(D,DiscordPair)
  range_trg_discord = [DiscordPair(1),DiscordPair(2)];
  range_NN_discord = [DiscordPair(3),DiscordPair(4)];
  times_trg_discord = (range_trg_discord(1):range_trg_discord(end));
  times_NN_discord = (range_NN_discord(1):range_NN_discord(end));
  Raw_trg = D(times_trg_discord); 
  Raw_NN =  D(times_NN_discord); 
  Znorm_trg_discord = normalize(Raw_trg);
  Znorm_NN_discord = normalize(Raw_NN);
end