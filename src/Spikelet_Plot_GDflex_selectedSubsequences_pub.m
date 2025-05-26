function Spikelet_Plot_GDflex_selectedSubsequences_pub(Info_selected)
        
        SegInfo = Info_selected.SegInfo;
        SegDB_org = SegInfo.SegDB_org;
        SegDB  = SegInfo.SegDB;
        Names = SegInfo.SegDB_names;
        D = SegInfo.data_org;
        RawData = Info_selected.RawData; 
        
        Row = 3; Column = 4; 
        Rate_X = 0.05; Rate_Y = 0.1; ROUND = 3;

        t = tiledlayout(Row,Column);
        %% plot histogram
         nexttile([1 3]);
        X = (1:length(D))'; Y = D;
        plot(X,Y); xlim([1 length(D)]);
        if exist('RawData','var') && isfield(RawData,'anomaly_region')
          Anomaly_region = RawData.anomaly_region;
          from = Anomaly_region(1); to = Anomaly_region(2);
          training_region = RawData.training_region;
          hold on;
          X_a = (from:to);
          Y_a = D(X_a);
          plot(X_a,Y_a,'r','LineWidth',2);
          xline(training_region,'--');
          hold off;
        end
         name = RawData.name;
         c22 = RawData.C22;
         damp = RawData.DAMP;
         window_length = RawData.window_length;
         Title = [name, ' (C22,damp)=(',num2str(c22),' ',num2str(damp),')',' win-len(', num2str(window_length),')'];
         title(Title);
      nexttile(4)
        if exist('RawData','var') && isfield(RawData,'anomaly_region')
          Anomaly_region = RawData.anomaly_region;
          from = Anomaly_region(1); to = Anomaly_region(2);
          len = to - from + 1;
          from_plot = max(from-len,1);  to_plot = min(to+len,length(D));
          X = (from:to);
          Y = D(X);
          X_plot = (from_plot:to_plot);
          Y_plot = D(X_plot); 
          plot(X_plot,Y_plot); hold on; plot(X,Y,'r','LineWidth',2);
          LIM_X = lim_margin(X_plot,Rate_X); if ~isempty(LIM_X); xlim(LIM_X); end
          LIM_Y = lim_margin(Y_plot,Rate_Y); if ~isempty(LIM_Y); ylim(LIM_Y); end 
          title(['anomaly region ', num2str(from),'-',num2str(to)]);
        end
        HistP = Column*1;
        nexttile(HistP+1)
          [H_afterPos_mag, Title] = histogram_MagSuppStatistics(SegDB,Names,2,"magnitude",ROUND);
          title(Title)
        nexttile(HistP+2)
          [H_afterPos_supp, Title] = histogram_MagSuppStatistics(SegDB,Names,2,"length",ROUND);
          title(Title)
        nexttile(HistP+3)
          [H_afterNeg_mag, Title] = histogram_MagSuppStatistics(SegDB,Names,-2,"magnitude",ROUND);
          title(Title)
        nexttile(HistP+4)
          [H_afterNeg_supp, Title] = histogram_MagSuppStatistics(SegDB,Names,-2,"length",ROUND);
          title(Title)
        %% plot instant
        InstP = Column*2;
        nexttile(InstP+1)
          plot_motifCandidatesWithAnomalyRegion(SegDB,Names,2,"center",D,Anomaly_region);
          title('positive after');
        nexttile(InstP+2)
          plot_motifCandidatesWithAnomalyRegion(SegDB,Names,2,"left",D,Anomaly_region);
          title('positive after');
        nexttile(InstP+3)
          plot_motifCandidatesWithAnomalyRegion(SegDB,Names,-2,"center",D,Anomaly_region);
          title('negative after');
        nexttile(InstP+4)
          plot_motifCandidatesWithAnomalyRegion(SegDB,Names,-2,"left",D,Anomaly_region);
          title('negative after');
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

% AlignOption = "center"; "left"
function plot_motifCandidatesWithAnomalyRegion(SegDB,Names,Type,AlignOption,D,AnomalyRegion)
    SegDB_org = SegDB;
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

function [H, Title] = histogram_MagSuppStatistics(SegDB,Names,Type,Feature,ROUND)
  Index = SegDB(:,Names=="type") == Type;
  SegDB_plot = SegDB(Index,:);
  ValList = SegDB(Index,Names==Feature);
  
  H = histogram(ValList);
  Title = magSuppStatistics_title(SegDB_plot,Names,Feature,ROUND);
end

% Option: magnitude,  support
function Title = magSuppStatistics_title(SegDB,Names,Option,ROUND)
  Info = statistics_SegDB(SegDB,Names);
  if Option == "magnitude"
   Title = [num2str(size(SegDB,1))...
           ,' mag(mean,std,median) = ' ...
           ,num2str(round_N(Info.magMean,ROUND)) ...
           ,' ',num2str(round_N(Info.magStd,ROUND)) ...
           ,' ',num2str(round_N(Info.magMedian,ROUND))];
  elseif Option == "length"
    Title = [' length(mean,std,median)  = ' ...
             ,' ',num2str(round_N(Info.suppMean,ROUND)) ...
             ,' ',num2str(round_N(Info.suppStd,ROUND)) ...
             ,' ',num2str(round_N(Info.suppMedian,ROUND))];
  end
end

function Info = statistics_SegDB(SegDB,Names)

  Info = struct;
  Info.suppLen = SegDB(:, Names=="length");
  Info.suppMean = mean(Info.suppLen);
  Info.suppStd = std(Info.suppLen);
  Info.suppMedian = median(Info.suppLen);
  Info.magnitude = SegDB(:, Names=="magnitude");
  Info.magMean = mean(Info.magnitude);
  Info.magStd = std(Info.magnitude);
  Info.magMedian = median(Info.magnitude);
end

function rounded_x = round_N(x,N)
  rounded_x = round(x, N - 1 - floor(log10(abs(x))));
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

