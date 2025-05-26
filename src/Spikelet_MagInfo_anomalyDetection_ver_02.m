% Param = struct; Param.PN = "positive"; Param.clustering_algorithm = "KsdensityGmm";
function Info = Spikelet_MagInfo_anomalyDetection_ver_02(MagInfo,Param,RawData)

    %% Clustering by only magnitude 
    Mag = MagInfo.magnitude; 
    SpikeId_all = find(Mag~=0);
    Mag_all = Mag(SpikeId_all);
   
    InitialCenter_positive = mean(abs(Mag_all)) + 2*std(abs(Mag_all));
    InitialCenters = [(-1)*InitialCenter_positive,0,InitialCenter_positive]';
    [idx_mag_all, C_mag_all] = kmeans(Mag_all, 3, 'Start', InitialCenters);
    std_mag_all = zeros(size(C_mag_all,1),size(C_mag_all,2));
    for cm = 1:size(C_mag_all,1)
      Mag_all_cm = Mag_all(idx_mag_all==cm);  
      std_mag_all(cm,:) = std(Mag_all_cm - mean(Mag_all_cm));
    end
    [min_val,min_pos] = min(abs(C_mag_all));
    SpikeId_filtered = SpikeId_all(idx_mag_all ~= min_pos); 

    if isfield(RawData,"training_region")
      training_region = RawData.training_region;
      Anomaly_region = RawData.anomaly_region;
    end

    %% training test
    posMSF_train_id = []; posMSF_test_id = [];  negMSF_train_id = []; negMSF_test_id = [];
    if exist('training_region','var')
        
      for i = 1:length(SpikeId_filtered)
        SpikeId_i = SpikeId_filtered(i);
        from_i = MagInfo.left(SpikeId_i); 
        to_i = MagInfo.right(SpikeId_i);
         % training=-1, test = 1;
        if from_i > training_region % test
          if MagInfo.magnitude(SpikeId_i) > 0 % positive
             posMSF_test_id = [posMSF_test_id; SpikeId_i]; 
          elseif MagInfo.magnitude(SpikeId_i) < 0 % negative
             negMSF_test_id = [negMSF_test_id; SpikeId_i];  
          end
        elseif to_i < training_region  % trainning
          if MagInfo.magnitude(SpikeId_i) > 0 % positive
             posMSF_train_id = [posMSF_train_id; SpikeId_i]; 
          elseif MagInfo.magnitude(SpikeId_i) < 0 % negative
             negMSF_train_id = [negMSF_train_id; SpikeId_i];  
          end
        end
      end
    end
  
    %% Spikelet_SupportRestriction_BySpikeInSpike
    Train_id_pos = posMSF_train_id; Train_id_neg = negMSF_train_id;
    Train_id = [Train_id_pos; Train_id_neg]; 

    Test_id_pos = posMSF_test_id; Test_id_neg = negMSF_test_id;
    Test_id = [Test_id_pos; Test_id_neg];  

    %% SegInfo test
    % SegDB_names = ["id","type","from","to","length","center","fluctuation"];
    SegInfo = Spikelet_SegInfo_generateFromSpike(MagInfo,Test_id);
    SegInfo = Spikelet_SegInfo_spikeSelect(SegInfo,Param);
    SegInfo.SegDB_initial = SegInfo.SegDB;

    %% Seg statistics test
    SegDB = SegInfo.SegDB;
    SegDB_names = SegInfo.SegDB_names;
    PositiveSeg = SegDB(SegDB(:,SegDB_names == "type") == 2,:);
    NegativeSeg = SegDB(SegDB(:,SegDB_names == "type") == -2,:);
    SegStatistics = struct;
    SegStatistics.positive2.MaT = min(PositiveSeg(:,SegDB_names=="magnitude"));
    SegStatistics.negative2.MaT = max(NegativeSeg(:,SegDB_names=="magnitude"));
    SegStatistics.positive2.MaT_median = median(PositiveSeg(:,SegDB_names=="magnitude"));
    SegStatistics.negative2.MaT_median = median(NegativeSeg(:,SegDB_names=="magnitude"));
    SegStatistics.positive2.MaT_std = std(PositiveSeg(:,SegDB_names=="magnitude"));
    SegStatistics.negative2.MaT_std = std(NegativeSeg(:,SegDB_names=="magnitude"));

    SegStatistics.positive2.Supp_median = median(PositiveSeg(:,SegDB_names=="length"));
    SegStatistics.negative2.Supp_median = median(NegativeSeg(:,SegDB_names=="length"));
    SegStatistics.positive2.Supp_std = std(PositiveSeg(:,SegDB_names=="length"));
    SegStatistics.negative2.Supp_std = std(NegativeSeg(:,SegDB_names=="length"));

    Param.SegStatistics = SegStatistics;
    %% (1) short scale expansion
    if Param.ablation.intra
      SegInfo = Spikelet_SegInfo_extractSpikeInSpike_ver2(SegInfo,Param); 
    end

    %% selected subsequences
    Info_selected = struct;
    Info_selected.RawData = RawData;
    Info_selected.testId_pos = Test_id_pos;
    Info_selected.testId_neg = Test_id_neg;
    Info_selected.testId = Test_id;
    Info_selected.SegInfo = SegInfo;
  
    %% anomaly detection
    if exist('training_region','var')
      %% anomaly detection% 
      DiscordRange = "train";  % DiscordRange = "all";
      if DiscordRange == "train"
        [AdInfo,AdInfo_names] = Spikelet_AD_TrTe_STOMP_ver_04(SegInfo,training_region,Param);
      elseif DiscordRange == "all"
        [AdInfo,AdInfo_names] = Spikelet_AD_TrTe_STOMP_ver_04(SegInfo,training_region,Param);
      end

      [dist_sorted_val, dist_sorted_id] = sort(AdInfo(:,AdInfo_names=="distLengthNorm"),"descend");
      AdInfo = AdInfo(dist_sorted_id,:); 
      AdInfo_1st = AdInfo;
    end

    %% SegInfo-train ( It is better to unify for SegInfo-test)
    if Param.ablation.locMis 
      SegInfo_train = Spikelet_SegInfo_generateFromSpike(MagInfo,Train_id); 
      SegInfo_train = Spikelet_SegInfo_spikeSelect(SegInfo_train,Param);
    end
    
    %% locMis 1
    if Param.ablation.locMis   
      [AdInfo,AdjustInfo_twoSpikeWithSameSign,AdjustInfo_twoSpikeWithSameSign_names,PairwiseInfo_LocMis1] = Spikelet_SegInfo_adjustAnomalyScore_ver_02(AdInfo,AdInfo_names,SegInfo_train,Param);
    end   

    if Param.ablation.inter
     %if (length(PeakValList_1st) >= 2 && PeakValList_1st(2) < PeakDtwVal_Thr) || length(PeakValList_1st) < 2  

        Param.ConsecutiveSpikeExpansion.ConsecutiveSpikeRate_1st = 0.1;
        Param.ConsecutiveSpikeExpansion.ConsecutiveSpikeRate_2nd = 0.2;
        Param.ConsecutiveSpikeExpansion.SharpConvexityByStd = 3;

        % revise SegInfo
        SegInfo_CSE = Spikelet_SegInfo_SegmentExpansionByconsecutiveSpike_ver2(AdInfo,AdInfo_names,SegInfo,Param);
        SegInfo.SegDB = [SegInfo.SegDB; SegInfo_CSE.SegDB];
        SegInfo.last_id = SegInfo_CSE.last_id;

        %% cut noise adjustment (revise AdInfo and SegInfo for SegInfo_CSE)
        % (1) calculate anomaly score for SegInfo_CSE
        [AdInfo_CSE,AdInfo_names] = Spikelet_AD_TrTe_STOMP_ver_04(SegInfo_CSE,training_region,Param);
        if Param.ablation.interNoise 
          [AdInfo_CSE_revised,SegInfo,Info_AS_CSE] = Spikelet_SegInfo_adjustPVAndNoiseCut(AdInfo_CSE,AdInfo_names,SegInfo,training_region,MagInfo.data_org,Param);
        else
          AdInfo_CSE_revised = AdInfo_CSE;
        end
        % revise AdInfo
        AdInfo = [AdInfo; AdInfo_CSE_revised];  
        [ValSorted, IndexSorted] = sort(AdInfo(:,AdInfo_names=="distLengthNorm"),"descend");
        AdInfo = AdInfo(IndexSorted,:);

    end

    %%
    Info = struct;
    Info.MagInfo = MagInfo;
    Info.Param = Param;
    Info.spikeId.all = SpikeId_all;
    Info.mag_all.kmeans.idx = idx_mag_all;
    Info.mag_all.kmeans.mu = C_mag_all; 
    Info.mag_all.kmeans.sigma = std_mag_all;
    Info.mag_all.data = Mag_all;

    Info.Info_selected = Info_selected;
   
    Info.SegInfo = SegInfo;
    if exist('SegInfo_train','var')
      Info.SegInfo_train = SegInfo_train; 
    end
    Info.anomalyInfo.names = AdInfo_names;
    Info.anomalyInfo.positive = AdInfo(AdInfo(:,AdInfo_names=="type")==2,:);
    Info.anomalyInfo.negative = AdInfo(AdInfo(:,AdInfo_names=="type")==-2,:);
    Info.anomalyInfo.expanded = AdInfo(abs(AdInfo(:,AdInfo_names=="type"))~=2,:);
    Info.anomalyInfo.total = AdInfo;
    Info.anomalyInfo.total_1st = AdInfo_1st;
    Info.anomalyInfo.training_region = training_region;
  
    if exist('AdjustInfo_twoSpikeWithSameSign','var')
      Info.adjust.twoSpikeWithSameSign = AdjustInfo_twoSpikeWithSameSign;
      Info.adjust.twoSpikeWithSameSign_names = AdjustInfo_twoSpikeWithSameSign_names;
    end
    if exist('AdjustInfo_ZnormBias','var')
      Info.adjust.AdjustInfo_ZnormBias = AdjustInfo_ZnormBias;
      Info.adjust.AdjustInfo_ZnormBias_names = AdjustInfo_ZnormBias_names;
    end

    if exist('PairwiseInfo_LocMis1','var')
      Info.adjust.PairwiseInfo_LocMis1 = PairwiseInfo_LocMis1;
    end

    if exist("figId_total","var")
      Info.figId_total = figId_total;
    end
end

function SegInfo = Spikelet_SegInfo_generateFromSpike(MagInfo,SpikeId)
  SegDB_names = ["id","type","from","to","length","center","magnitude","fluctuation","component_1","component_2"];
  SegDB = NaN(length(SpikeId), length(SegDB_names));
  D = MagInfo.data_org;
  for i = 1:length(SpikeId)
    spikeId_i = SpikeId(i);
    % if spikeId_i == 29342; fprintf('debug\n'); end
    mag_i = MagInfo.magnitude(spikeId_i);
    from_i = MagInfo.left(spikeId_i);
    to_i = MagInfo.right(spikeId_i);
    % if from_i == 5864; DEBUG = 0; end
    fluctuation_i = max(D(from_i:to_i)) - min(D(from_i:to_i)); 
    length_i = to_i- from_i + 1;
    segDB_i = [spikeId_i,2*sign(mag_i),from_i,to_i,length_i,spikeId_i,mag_i,fluctuation_i];
    SegDB(i,1:length(segDB_i)) = segDB_i;
  end 
  
  SegInfo = struct;
  SegInfo.data_org = MagInfo.data_org;
  SegInfo.MagInfo = MagInfo;
  SegInfo.SegDB_org = SegDB;
  SegInfo.SegDB = SegDB;
  SegInfo.SegDB_names = SegDB_names;
  SegInfo.last_id = max([SegDB(:,SegDB_names=='id');length(MagInfo.data_org)]);
end

function SegInfo = Spikelet_SegInfo_spikeSelect(SegInfo,Param)
      Param_AD = Param.baseline.spikeSelection;
    
      SegInfo = Spikelet_SegInfo_restrictSupportBySpikeInSpike_ver2(SegInfo,Param_AD); 
      SegInfo = Spikelet_SegInfo_deleteSpikesByMedimStisticsWithPV_ver2(SegInfo,Param_AD); 
      SegInfo = Spikelet_SegInfo_deleteSegmentsOfSuppLenThree(SegInfo,Param_AD); 
end


