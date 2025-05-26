function Demo_Spikelet5_KDD2025(TargetData_option,Ablation)
  %{
  Demo_Spikelet5_KDD2025(66,"interNoise");
  Demo_Spikelet5_KDD2025(91,"interNoise");
  Demo_Spikelet5_KDD2025("all","interNoise");

  ABlist=["baseline", "znormBias","intra","locMis","inter","interNoise"]  
  for i = 1:length(ABlist)
    Demo_Spikelet5(66,ABlist(i));
  end
  %}

  %% directory 
  FuncName = 'Demo_Spikelet5_KDD2025';
  Data_Com_Dir = 'C:/Users/imamura/0_data';
  TSAD2021_Dir = [Data_Com_Dir,'/2_open/UCR_TimeSeriesAnomalyDatasets2021/mat'];
  Apl_baseDir = ['C:/Users/imamura/1_Development/Matlab_Exec/application/Spikelet_5/Work/',FuncName]; FILE_createDir(Apl_baseDir);

  if isnumeric(TargetData_option)
    PLOT = true;  
    InfoSave = true;
  else
    PLOT = false;  
    InfoSave = false;
  end

  %% Parameters for spikelet decomposition
  Arg = struct;
  Arg.dist = 'Euc';
  Arg.subseq_norm = 'Znorm';
  Arg.SuppCutMaTr = 0.6;
  Arg.LegSuppCutMaTrFM = 0.9;

  %% Parameters for spikelet selection
  Arg.baseline.EXCLUSION_WINDOW_RATE = 0.2;
  Arg.baseline.spikeSelection.SR_ByMagRatio = 0.75;
  Arg.baseline.spikeSelection.SR_BySuppRatio = 0.75; 
  Arg.baseline.spikeSelection.MmedianSupp = 3;
  Arg.baseline.spikeSelection.MmedianSupp_upperLimit = 10;
  
  %% znormBias
  Arg.znormBias.FlacToFlac_0 = [0.35, 2.5];
  Arg.znormBias.FlacToFlac_1 = [0.4, 2]; 
  Arg.znormBias.FlacToFlac_2 = [0.75,1.5];
  Arg.znormBias.FlacToFlac_3 = [0.75,2.5];
  Arg.znormBias.FlacToFlac_5 = [0.5,2]; 
  Arg.znormBias.MatToMat_2 = [0.5,2];  

  %% Intra
  Arg.intra.Const_MagBandRate = 0.125;     
  Arg.intra.Const_LengthRate = 0.2; % 0.5  
  Arg.intra.Opp_NoiseThresRate = 0.2;    
  Arg.intra.Opp_VibLengthRate = 0.5;     
  Arg.intra.Opp_ConsecuteMagRate = 0.5;  
  Arg.intra.Opp_VibPVnum = 5;           

  %% Localized Mismatch
  Arg.locMis.SubseqDiffMax_Thr = 0.5;
  Arg.locMis.SuppMedian_ratio = 1.5;
  Arg.locMis.DtwEuc_ratio = 0.5;  % Epsiron_2: 0.3 is also OK
  Arg.locMis.BoudaryRaw_absDiffLr_ratio = 0.3; % Epsiron_2:0.1, 0.89 OK
  Arg.locMis.MaxDiffFromCenter_ratio = 0.2;% Epsiron_6
  Arg.locMis.DistTop1_ratio = 0.5;             % Epsiron_7
  Arg.locMis.DistanceTop_Thr = 1;     % Epsiron_7
  Arg.locMis.LocMisFluxRatio_Thr = 0.75; % Epsiron_7pre3 (pre2 1)

  %% ablation study
  if Ablation == "baseline"
    Arg.ablation.znormBias = false;  
    Arg.ablation.intra = false;  % Intra
    Arg.ablation.inter = false;  % Inter
    Arg.ablation.locMis = false;    % LocMis 1: twoSpikeWithSameSign 
    Arg.ablation.interNoise = false;   % InterNoise 
  elseif Ablation == "znormBias"
    Arg.ablation.znormBias = true;  
    Arg.ablation.intra = false; 
    Arg.ablation.inter = false; 
    Arg.ablation.locMis = false;   
    % Arg.ablation.mitigateBias_Znorm = false;         
    Arg.ablation.interNoise = false;  
  elseif Ablation == "intra"
    Arg.ablation.znormBias = true;  
    Arg.ablation.intra = true; 
    Arg.ablation.inter = false; 
    Arg.ablation.locMis = false;          
    Arg.ablation.internoise = false;  
  %{  
  elseif Ablation == "PV2tc_addPV34"  
    ArgPlus.ablation.testSetConstraint = true;  
    ArgPlus.ablation.expansion_shortTimeScale = false; 
    ArgPlus.ablation.expansion_largeTimeScale = true;  
    ArgPlus.ablation.mitigateBias_euclidean = false;           
    ArgPlus.ablation.largeTimeScale_noiseCut = false;  
  %}
  elseif Ablation == "locMis"
    Arg.ablation.znormBias = true;  
    Arg.ablation.intra = true;  
    Arg.ablation.inter = false; 
    Arg.ablation.locMis = true;    
    % Arg.ablation.mitigateBias_Znorm = false;         
    Arg.ablation.interNoise = false; 
  %{
  elseif Ablation == "PV2tc_addPV0134"
    ArgPlus.ablation.testSetConstraint = true;  
    ArgPlus.ablation.expansion_shortTimeScale = true;  
    ArgPlus.ablation.expansion_largeTimeScale = true; 
    ArgPlus.ablation.mitigateBias_euclidean = false;           
    ArgPlus.ablation.largeTimeScale_noiseCut = false;  
  %}
  elseif Ablation == "inter"
    Arg.ablation.znormBias = true;  
    Arg.ablation.intra = true;  
    Arg.ablation.inter = true; 
    Arg.ablation.locMis = true;         
    Arg.ablation.interNoise = false;  
  elseif Ablation == "interNoise"  
    Arg.ablation.znormBias = true;  
    Arg.ablation.intra = true; 
    Arg.ablation.inter = true; 
    Arg.ablation.locMis = true;         
    Arg.ablation.interNoise = true;  
  end
     
   %% Success/Failure table
   RawDataList = dir([TSAD2021_Dir,'/*.mat']);  SF = NaN(250,2); SF_names = ["DAMP","C22"];
   for i = 1:length(RawDataList)
     rawData_i = [TSAD2021_Dir,'/',RawDataList(i).name];  
     RawData_i = load(rawData_i);
     id_i = RawData_i.id;
     C22_i = RawData_i.C22;
     damp_i = RawData_i.DAMP;
     if id_i ~= i
       fprintf('[error] data id inconsistent in creating SF table. \n'); 
       return;
     else
       SF(i,:) = [damp_i, C22_i];  
     end
   end
   Success_both = find( SF(:,SF_names=="DAMP") == true & SF(:,SF_names=="C22") == true);
   Fail_both    = find( SF(:,SF_names=="DAMP") == false & SF(:,SF_names=="C22") == false);
   OnlySuccess_DAMP= find( SF(:,SF_names=="DAMP") == true & SF(:,SF_names=="C22") == false);
   OnlySuccess_C22 = find( SF(:,SF_names=="DAMP") == false & SF(:,SF_names=="C22") == true);

   if isstring(TargetData_option) && TargetData_option == "all"
       TargetData = (1:size(SF,1));  % TargetData = (1:2);  
       TargetData_str = TargetData_option;
   else
     TargetData = TargetData_option; 
     TargetData_str = string([num2str(TargetData_option(1)),'-',num2str(length(TargetData_option))]);
   end
   
   %% argument
   close all;
   ParamAprox = Arg2SpikeletParam(Arg);
  
   %% main
   RawDataList = dir([TSAD2021_Dir,'/*.mat']);

   %% anomaly detection
   TimePN = zeros(length(TargetData),1);
   RSLT_names = ["data_id","spike_id","spikelet","inRange","rank","damp","c22","anomaly_label_from","anomaly_label_to","window_label"...
                ,"spikelet_type","spikelet_from","spikelet_to","from_NN","to_NN","window_spikelet","distLengthNorm"...
                ,"numSpikelet","numSlidingWindow","time","total_length","testNumRate","adjust"];
   RSLT = zeros(length(TargetData),length(RSLT_names));
   % for i = 1:1
   for i = 1:length(TargetData)  % 12, 25, 43
     % load data  
     dataId_i = TargetData(i);
     rawData_i = [TSAD2021_Dir,'/',RawDataList(dataId_i).name];  

     RawData_i = load(rawData_i);
     D_i = RawData_i.data;
     %  figure(); X = [19001:20000]; plot(X,D_i(X));
     name_i = RawData_i.name;
     id_i = RawData_i.id;
     if dataId_i ~= id_i
       fprintf('[error] data id inconsistent in main loop. \n');    
       return;  
     end
     window_length_i = RawData_i.window_length;
     C22_i = RawData_i.C22;
     damp_i = RawData_i.DAMP;
     anomaly_region =  RawData_i.anomaly_region;
     training_region = RawData_i.training_region;

     %% print
     dataTitle_i = ['[',num2str(dataId_i),'] ',name_i, ' (C22,damp)=(',num2str(C22_i),' ',num2str(damp_i),')'...
             ,' (winLen,train)=(', num2str(window_length_i),' ',num2str(training_region),')'...
             ,' a-region(',num2str(anomaly_region(1)),' ',num2str(anomaly_region(2)),')'];
     fprintf('\n %s\n',dataTitle_i);

     %% Anomaly detection
     startTime_PN_i = tic; 
     MagInfo_i = Spikelet_aproximation_ver_05(D_i, ParamAprox);
     InfoPosNeg = Spikelet_MagInfo_anomalyDetection_ver_02(MagInfo_i,Arg,RawData_i);   

     if InfoSave
       OutputDir = [Apl_baseDir,'/',num2str(dataId_i)];
       OutputDir = Apl_baseDir;
       FILE_createDir(OutputDir);
       OutputFile = [OutputDir,'/','InfoDemo_',num2str(dataId_i),'.mat']; 
       % InfoPosNeg_org = InfoPosNeg;
       % InfoPosNeg = rmfield(InfoPosNeg, 'figId_total');
       save(OutputFile, 'InfoPosNeg');
       % InfoPosNeg = InfoPosNeg_org;
     end

     %% result
     Anomaly_names = InfoPosNeg.anomalyInfo.names;
     spike_id = InfoPosNeg.anomalyInfo.total(1,Anomaly_names=="id");
     top1_from = InfoPosNeg.anomalyInfo.total(1,Anomaly_names=="from");
     top1_to = InfoPosNeg.anomalyInfo.total(1,Anomaly_names=="to");
     top1_type = InfoPosNeg.anomalyInfo.total(1,Anomaly_names=="type");
     top1_fromNN = InfoPosNeg.anomalyInfo.total(1,Anomaly_names=="from_NN");
     top1_toNN = InfoPosNeg.anomalyInfo.total(1,Anomaly_names=="to_NN");
     top1_window = top1_to - top1_from + 1; 
     top1_distLengthNorm = InfoPosNeg.anomalyInfo.total(1,Anomaly_names=="distLengthNorm");
     
     %% UCR_AnomalyDataSets.ppt p.6,p.7 
     % mean of the region
     % center_top1 = (top1_from+top1_to)/2;
     % nn neighbor weighted center of the region
     center_top1 = weighted_center_L2(D_i,InfoPosNeg.anomalyInfo.total(1,:),Anomaly_names);
     if max(1,anomaly_region(1)-100) <= center_top1 && center_top1 <= min(anomaly_region(2)+100,length(D_i))
       SF_spikelet_i = true;
     else
       SF_spikelet_i = false;
     end
     % Evaluation of Interval Overlaps
     % SF_spikelet_i = Spikelet_overlapSegment([top1_from, top1_to], anomaly_region);
     FromTo_all = InfoPosNeg.anomalyInfo.total(:,[find(Anomaly_names=="from"),find(Anomaly_names=="to")]);
     [inRange_i, inRange_rank_i, inRange_rate_i] = overlapSegment_map(FromTo_all,anomaly_region);
     %{
     RSLT_names = ["data_id","spike_id","spikelet","inRange","rank","damp","c22","anomaly_label_from","anomaly_label_to","window_label"...
                ,"spikelet_type","spikelet_from","spikelet_to","from_NN","to_NN","window_spikelet","distLengthNorm"...
                ,"numSpikelet","numSlidingWindow","time","testNumRate"];
     %}
     numSpikelet_i = size(InfoPosNeg.anomalyInfo.total,1);
     testNum_rate_i = InfoPosNeg.anomalyInfo.total(1,Anomaly_names=="testNumRate");
     adjust_i = InfoPosNeg.anomalyInfo.total(1,Anomaly_names=="adjust");
     numSlidingWindow_i = (length(D_i) - training_region + 1) - window_length_i + 1;
     time_PN_i = toc(startTime_PN_i);

     RSLT(i,:) =  [dataId_i, spike_id, SF_spikelet_i,inRange_i, inRange_rank_i, damp_i,C22_i,anomaly_region(1),   anomaly_region(2), window_length_i...
                ,top1_type,       top1_from,     top1_to,     top1_fromNN,  top1_toNN,   top1_window,      top1_distLengthNorm...
                ,numSpikelet_i, numSlidingWindow_i, time_PN_i, length(D_i),testNum_rate_i,adjust_i];
 
     
     % time
     TimePN(i) = time_PN_i;

     %% plot
     % PLOT_ANOMLY = false;
     if  PLOT
       PlotOnTopK_total = false; PlotOnTopK_4 = false; PlotOnTopK_3 = false; PlotOnTopK_pos2 = false;  PlotOnTopK_neg2 = false; 

       %% AdInfo
       AdInfo_names = InfoPosNeg.anomalyInfo.names;
       AdInfo_Total = InfoPosNeg.anomalyInfo.total;
       AdInfo_Pos = InfoPosNeg.anomalyInfo.positive;
       AdInfo_Neg = InfoPosNeg.anomalyInfo.negative;
       AdInfo_Expanded = InfoPosNeg.anomalyInfo.expanded;

       AdInfo_names_pre = AdInfo_names;
       [AdInfo_Total,AdInfo_names] = AdInfo_addAnomalyRegion(AdInfo_Total,AdInfo_names_pre,anomaly_region,id_i);
       [AdInfo_Pos,AdInfo_names] = AdInfo_addAnomalyRegion(AdInfo_Pos,AdInfo_names_pre,anomaly_region,id_i);
       [AdInfo_Neg,AdInfo_names] = AdInfo_addAnomalyRegion(AdInfo_Neg,AdInfo_names_pre,anomaly_region,id_i);
       [AdInfo_Expanded,AdInfo_names] = AdInfo_addAnomalyRegion(AdInfo_Expanded,AdInfo_names_pre,anomaly_region,id_i);

       %% plot total figure
       Info_selected = InfoPosNeg.Info_selected;
       fig_id_total = figure();
         Spikelet_Plot_GDflex_selectedSubsequences_pub(Info_selected);
         Title = dataTitle_i;
         Fig_figtitle(strrep(Title,'_','\_'));    

       fig_id = figure();
         fprintf('total\n');
         plot_anomaly(AdInfo_Total,AdInfo_names,D_i,anomaly_region);
         Title = ['total ',dataTitle_i];
         Fig_figtitle(strrep(Title,'_','\_'));

        fig_id = figure();
         fprintf('positive\n');
         plot_anomaly(AdInfo_Pos,AdInfo_names,D_i,anomaly_region);
         Title = ['positive ',dataTitle_i];
         Fig_figtitle(strrep(Title,'_','\_'));

        fig_id = figure();
         fprintf('negative\n');
         plot_anomaly(AdInfo_Neg,AdInfo_names,D_i,anomaly_region);
         Title = ['negative ', dataTitle_i];
         Fig_figtitle(strrep(Title,'_','\_'));
     
        fig_id = figure();
         fprintf('expanded\n');
         plot_anomaly(AdInfo_Expanded,AdInfo_names,D_i,anomaly_region);
         Title = ['expanded ', dataTitle_i];
         Fig_figtitle(strrep(Title,'_','\_'));

     end
       
   end
  
  %% sum
  % RSLTrate_names = ["id_num","spikelet","inRange","damp","c22"];
  RSLTrate_names = ["id_num","spikelet","damp","c22"];
  RSLTrate = NaN(1,length(RSLTrate_names));
  RSLTrate(:,RSLTrate_names=="id_num") = size(RSLT,1);
  RSLTrate(:,RSLTrate_names=="spikelet") = sum(RSLT(:,RSLT_names=="spikelet"))/size(RSLT,1);
  % RSLTrate(:,RSLTrate_names=="inRange") = sum(RSLT(:,RSLT_names=="inRange"))/size(RSLT,1);
  RSLTrate(:,RSLTrate_names=="c22") = sum(RSLT(:,RSLT_names=="c22"))/size(RSLT,1);
  RSLTrate(:,RSLTrate_names=="damp") = sum(RSLT(:,RSLT_names=="damp"))/size(RSLT,1);

  %% save total output 
  TargetData_char = char(TargetData_str);
  % RSLT
  % OutputDir = [Apl_baseDir,'/result']; FILE_createDir(OutputDir);
  OutputDir = Apl_baseDir;
  T = array2table(RSLT, 'VariableNames', cellstr(RSLT_names));
  outputFile = [OutputDir,'/result_',TargetData_char,'.csv'];
  writetable(T, outputFile);
  
  % RSLTrate
  T = array2table(RSLTrate, 'VariableNames',cellstr(RSLTrate_names));
  outputFile = [OutputDir,'/rate_',TargetData_char,'.csv'];
  writetable(T, outputFile);
end

function  center_i = weighted_center_L2(D_i,top1,Names)
   trg_from_i =  top1(Names=="from");
   trg_to_i =  top1(Names=="to");
   length_i = trg_to_i - trg_from_i + 1;
   nn_from_i = top1(Names=="from_NN");
   nn_to_i = top1(Names=="to_NN");
   trg_znorm_i = normalize(D_i(trg_from_i:trg_to_i));
   nn_znorm_i = normalize(D_i(nn_from_i:nn_to_i));
   diff_i = trg_znorm_i - nn_znorm_i;
   center_rel_i = sum((diff_i.^2).*(1:length_i)')/sum(diff_i.^2);
   center_i = trg_from_i + center_rel_i - 1;
end

function [AdInfo_new,AdInfo_names_new] = AdInfo_addAnomalyRegion(AdInfo,AdInfo_names,anomaly_region,id_i)
   AdInfo_names_new = ["data_id","SF","anomaly_label_from","anomaly_label_to",AdInfo_names];
   AdInfo_new = NaN(size(AdInfo,1),length(AdInfo_names_new));

   for i = 1:size(AdInfo,1)
     from_i = AdInfo(i,AdInfo_names=="from");
     to_i = AdInfo(i,AdInfo_names=="to");

     [OVERLAP, overlap_length,overlap_range] = overlapping_range([from_i to_i],anomaly_region);
     AdInfo_new(i,:) = [id_i,OVERLAP,anomaly_region(1),anomaly_region(2),AdInfo(i,:)];
     %{
     if OVERLAP
       fromOv_i = overlap_range(1);toOv_i = overlap_range(2); lenOv_i = toOv_i - fromOv_i + 1;
       len_i = to_i - from_i + 1;  
       AdInfo_new(i,:) = [id_i,anomaly_region(1),anomaly_region(2),AdInfo(i,:)];
     else
       AdInfo_new(i,:) = [id_i,anomaly_region(1),anomaly_region(2),0,0,0];
     end
     %}
   end
end
function [inRange, rank, rate] = overlapSegment_map(FromTo_all,anomaly_region)
  inRange = false; rank = 0; rate = 0;
  for i =1:size(FromTo_all,1)
    [inRange_i, InRange_length_i] = Spikelet_overlapSegment([FromTo_all(i,1),FromTo_all(i,2)],anomaly_region);
    if inRange_i == true
      inRange  = true;
      rank = i;
      rate = InRange_length_i/(anomaly_region(2)-anomaly_region(1)+1);
      return;
    end
  end
end

function plot_anomaly(AdInfo,AdInfo_names,D_all, anomaly_region)
  R = 3;
  ROUND = ['%.',num2str(R),'g'];
  % plot
  ROW = 4; COLUMN = 5;  RANK = COLUMN; distLengthZNorm_org = zeros(RANK,1); distLengthZNorm = zeros(RANK,1); distLength = zeros(RANK,1);
  for i = 1:min(RANK,size(AdInfo,1))
      type_i = AdInfo(i,AdInfo_names=="type");
      from = AdInfo(i,AdInfo_names=="from"); to = AdInfo(i,AdInfo_names=="to");
      fromNN = AdInfo(i,AdInfo_names=="from_NN"); toNN = AdInfo(i,AdInfo_names=="to_NN");
      distLengthZNorm(i) = AdInfo(i,AdInfo_names=="distLengthNorm");
      T = (from:to); Tnn = (fromNN:toNN); 
      TestD_norm = normalize(D_all(T)); TrainD_normNN = normalize(D_all(Tnn));
      TestD = D_all(T); TrainD = D_all(Tnn);

      distLengthZNorm_org(i) = norm(TestD_norm-TrainD_normNN)/sqrt(length(T));
      distLength(i) = norm(TestD-TrainD)/sqrt(length(T));
      % fprintf('%f <- %f',distLengthZNorm(i),distLength(i));

      [OVERLAP, overlap_length,overlap_range] = overlapping_range([from to],anomaly_region);
      
      subplot(ROW,COLUMN,i);
      plot(T,TestD_norm); hold on; plot(T,TrainD_normNN); xlim([from to]);

      Title = ['top ', num2str(i),': ',num2str(type_i),' ',num2str(length(T)),' [',num2str(from), ' ',num2str(to),'] [',num2str(fromNN),' ',num2str(toNN),'] ',num2str(distLengthZNorm_org(i),ROUND),' ',num2str(AdInfo(i,AdInfo_names=="distLengthNorm"),ROUND)]; 
      if OVERLAP
        title(Title,'Color', 'b', 'FontWeight', 'bold');  
      else
        title(Title);  
      end

      subplot(ROW,COLUMN,COLUMN+i);
      plot(T,TestD); hold on; plot(T,TrainD); xlim([from to]);
      if ~isempty(overlap_range)
        Tol = (overlap_range(1):overlap_range(2)); 
        plot(Tol,D_all(Tol),'b','LineWidth',2);
      end
      Title = ['top ', num2str(i),': ',num2str(type_i),' ',num2str(length(T)),' [',num2str(from), ' ',num2str(to),'] [',num2str(fromNN),' ',num2str(toNN),'] ',num2str(distLength(i),ROUND)]; 
      if OVERLAP
        title(Title,'Color', 'b', 'FontWeight', 'bold');  
      else
        title(Title);  
      end
  end
  disp(distLengthZNorm);
  %%
  Ov_list = find(AdInfo(:,AdInfo_names=="SF"));
  for i = 1:min(length(Ov_list),COLUMN)
      rank_i = Ov_list(i);   
      type_i = AdInfo(rank_i,AdInfo_names=="type");
      from = AdInfo(rank_i,AdInfo_names=="from"); to = AdInfo(rank_i,AdInfo_names=="to");
      fromNN = AdInfo(rank_i,AdInfo_names=="from_NN"); toNN = AdInfo(rank_i,AdInfo_names=="to_NN");
      T = (from:to); Tnn = (fromNN:toNN); 
      TestD_norm = normalize(D_all(T)); TrainD_normNN = normalize(D_all(Tnn));
      TestD = D_all(T); TrainD = D_all(Tnn);

      distLengthZNorm_i = norm(TestD_norm-TrainD_normNN)/sqrt(length(T));
      distLength_i = norm(TestD-TrainD)/sqrt(length(T));
      
      [OVERLAP, overlap_length,overlap_range] = overlapping_range([from to],anomaly_region);

    subplot(ROW,COLUMN,COLUMN*2+i)
      plot(T,TestD_norm); hold on; plot(T,TrainD_normNN); xlim([from to]);
      Title = ['top ', num2str(rank_i),': ',num2str(type_i),' ',num2str(length(T)),' [',num2str(from), ' ',num2str(to),'] [',num2str(fromNN),' ',num2str(toNN),'] ',num2str(distLengthZNorm_i,ROUND),' ',num2str(AdInfo(rank_i,AdInfo_names=="distLengthNorm"),ROUND)]; 
      if OVERLAP
        title(Title,'Color', 'b', 'FontWeight', 'bold');  
      else
        title(Title);  
      end

    subplot(ROW,COLUMN,COLUMN*3+i)
      plot(T,TestD); hold on; plot(T,TrainD); xlim([from to]);
      if ~isempty(overlap_range)
        Tol = (overlap_range(1):overlap_range(2)); 
        plot(Tol,D_all(Tol),'b','LineWidth',2);
      end
      Title = ['top ', num2str(rank_i),': ',num2str(type_i),' ',num2str(length(T)),' [',num2str(from), ' ',num2str(to),'] [',num2str(fromNN),' ',num2str(toNN),'] ',num2str(distLength_i,ROUND),' ',num2str(AdInfo(rank_i,AdInfo_names=="distLengthNorm"),ROUND)]; 
      if OVERLAP
        title(Title,'Color', 'b', 'FontWeight', 'bold');  
      else
        title(Title);  
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

function ParamAprx = Arg2SpikeletParam(Arg)
  % The order of the fields in if-statements is the same as the order of operations.
  Operation = struct;
  OS = {};  OSname = [];
  if isfield(Arg,"SuppCutMaTr")
    OS =  [OS,'restrictSupportByMagnitudeRatio'];  
    Operation.restrictSupportByMagnitudeRatio.magnitude_ratio = Arg.SuppCutMaTr;
    OSname = [OSname,'_RSMR'];
  end
  if isfield(Arg,"LegSuppCutMaTrFM")
    OS = [OS,'restrictSupportByLegMagnitudeRatioFirstMinimum'];   
    OSname = [OSname,'_RSLMRFM'];  
    Operation.restrictSupportByLegMagnitudeRatioFirstMinimum.magnitude_ratio = Arg.LegSuppCutMaTrFM;
  end
  Operation.operation_sequence = OS;
  Operation.operation_name = OSname;
  ParamAprx = struct;
  ParamAprx.operation = Operation;
end
