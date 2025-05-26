function  [result_i,distLengthNorm_i,boundary_NN_i,fromAT_NN_i,toAT_NN_i,fromAT_test_i,toAT_test_i,Info_PCSS_i] ...
    = Spikelet_Pairwise_containTwoSpikesWithSameSign_ver_03(D,type_i,fromTest_i,toTest_i,fromNN_i,toNN_i,SegDB_2,Names,distLengthNorm_top1,Param)
  %{
  Index_Seg2 = abs(SegDB(:,SegNames=="type")) == 2;
  SegDB_2 = SegDB(Index_Seg2,:);
  %}

  FuncName = 'Spikelet_Pairwise_containTwoSpikesWithSameSign_ver_03';
  fig_id = []; 
  DEBUG = true; % DEBUG = false;
  PLOT = true;  PLOT = false;

  Info_PCSS_i = struct;
  Info_PCSS_i.names = ["lengthTest","max_subseqDiff","boudaryRaw_absDiffLr_raio","DtwEuc_ratio"];
  Info_PCSS_i.values = [];

  if isempty(D)
    result_i = []; distLengthNorm_i = []; boundary_NN_i = []; fromAT_NN_i = []; toAT_NN_i=[]; fromAT_test_i=[]; toAT_test_i = []; 
    return;
  end
  %% Parameters
  %{
  ComRange_Thr = Param.LocMis_1.ComRange_Thr;
  BoundaryCloseness_Thr = Param.LocMis_1.BoundaryCloseness_Thr;
  SubseqDiffMin_Thr = Param.LocMis_1.SubseqDiffMin_Thr; 
  DiffFlux_ratio = Param.LocMis_1.DiffFlux_ratio; % Condition for the test subsequence to be relatively flat (using DTW)
  %}
  DistTop1_ratio_Thr = Param.locMis.DistTop1_ratio;
  SubseqDiffMax_Thr = Param.locMis.SubseqDiffMax_Thr;
  SuppMedian_ratio = Param.locMis.SuppMedian_ratio;
  DtwEuc_ratio_Thr = Param.locMis.DtwEuc_ratio;  
  BoudaryRaw_absDiffLr_raio_Thr = Param.locMis.BoudaryRaw_absDiffLr_ratio; 
  MaxDiffFromCenter_ratio = Param.locMis.MaxDiffFromCenter_ratio;
  DistanceTop_Thr  = Param.locMis.DistanceTop_Thr;
  LocMisFluxRatio_Thr  = Param.locMis.LocMisFluxRatio_Thr;

  % Method_Attention = "DtwZeroCross"; 

  %% default return values
  result_i = false; 
  distLengthNorm_i = NaN; boundary_NN_i = NaN; fromAT_NN_i = NaN; toAT_NN_i = NaN; fromAT_test_i = NaN; toAT_test_i = NaN; distLengthNormORG_i = NaN;

  %% debug
  %{
  degre locMis1_new
11  21469	21652	5699	5882
12  21469	21652	12955	13138
119 21469	21652	5699	5882
120 21469	21652	12955	13138
151 14106	14460	1710	2064
222 110521	110708	49441	49628
225 81544	81713	75130	75299
231 55721	55965	2367	2611
232 55721	55965	2367	2611

       53  2765, 2996 (3467, 3586)
       57  15980, 16167
       66  7981 8034
       104 15980 16167
       161 2765 2996
       165 15980, 16167
       247 121735, 122113 (54853	55750)

  %}
  if DEBUG
    DEBUG_range = [21469, 21652];  % [2765,2996]; [2833,2942]; 152[22813, 23160]
    % [fromTest_i,toTest_i]
    if isequal([fromTest_i,toTest_i],DEBUG_range) 
      fprintf('[%s] debug start %d-%d\n',FuncName,fromTest_i,toTest_i);
    end 
  else
    DEBUG_range = [0,0];   
  end

  %% argument
  if type_i > 0
    SuppMedian = Param.SegStatistics.positive2.Supp_median;
  else
    SuppMedian = Param.SegStatistics.negative2.Supp_median;  
  end

  %% center
  lengthNN_org = toNN_i - fromNN_i + 1;
  centerNN_org = (fromNN_i + (lengthNN_org+1)/2 - 1);
  lengthTest = toTest_i - fromTest_i + 1;
  halfLengthTest = floor((lengthTest+1)/2);
  centerTest = (fromTest_i + (lengthTest+1)/2 - 1);

  %% center
  SegCenter_all = SegDB_2(:,Names=="center");
  Index_withSameSign = (SegDB_2(:,Names=="type") == type_i);
  SegCenter_sameSign = SegCenter_all(Index_withSameSign,:);
  SegDB_2_sameSign = SegDB_2(Index_withSameSign,:);

  % z-normalization
  subseqNN_i = D(fromNN_i:toNN_i);
  subseqNN_norm_i = normalize(subseqNN_i);
  subseqTrg_i = D(fromTest_i:toTest_i);
  subseqTrg_norm_i = normalize(subseqTrg_i);
  subseqLen_i = length(subseqTrg_i);
  subseqDiff_norm_i = subseqTrg_norm_i - subseqNN_norm_i; 
  subseqDiff_raw_i =  subseqTrg_i - subseqNN_i;    
  Test_raw_fluc = max(subseqTrg_i ) - min(subseqTrg_i );
  % boudaryRaw_absDiffLr_raio_i = abs(abs(subseqDiff_raw_i(1)) - abs(subseqDiff_raw_i(end)))/Test_raw_fluc; 
  boudaryRaw_absDiffLr_raio_i = max(abs(subseqDiff_raw_i(1)), abs(subseqDiff_raw_i(end)))/Test_raw_fluc;
  
  distLenNorm_i = norm(subseqDiff_norm_i)/sqrt(length(subseqDiff_norm_i));
  distTop1_ratio_i = distLenNorm_i/distLengthNorm_top1;

  %% main
  SegIdList = include([fromNN_i,toNN_i],SegCenter_sameSign);
  if length(SegIdList) >= 2 ...
     && distTop1_ratio_i >= DistTop1_ratio_Thr ...
     && lengthTest > SuppMedian*SuppMedian_ratio ...
     && max(subseqDiff_norm_i) >= SubseqDiffMax_Thr ...
     && boudaryRaw_absDiffLr_raio_i <= BoudaryRaw_absDiffLr_raio_Thr

     dtw_dist_i = sqrt(dtw(subseqTrg_norm_i,subseqNN_norm_i,'squared'));
     euc_dist_i = norm(subseqDiff_norm_i);
     dtwEuc_ratio_i = dtw_dist_i / euc_dist_i;                 
     if dtwEuc_ratio_i >= DtwEuc_ratio_Thr 
  
       FeatureList = ["ratio_insideL2_zeroCross_dtwRaw","ratio_insideL2_zeroCross_dtwZnorm"...
                  "ratio_insideL2Top2DiffSign_zeroCross_eucRaw","ratio_insideL2Top2DiffSign_zeroCross_eucZnorm"...
                  "ratio_10pL2_dtwZnorm","diff_magAbsMax_dtwZnorm"];            
       Param_DSF = struct;
       Param_DSF.FeatureList = FeatureList;
       [Feature_DSF,DSF_names,AdFeature,AdFeature_names] ...
           = Spikelet_SubseqPair_diffSpikeFeature(subseqTrg_i, subseqNN_i, Param_DSF);

       from_top1 = AdFeature(AdFeature_names=="zeroCrossTop1_from_eucZnorm");	
       to_top1 = AdFeature(AdFeature_names=="zeroCrossTop1_to_eucZnorm");
       from_top2 = AdFeature(AdFeature_names=="zeroCrossTop2_from_eucZnorm");	
       to_top2 = AdFeature(AdFeature_names=="zeroCrossTop2_to_eucZnorm");
       to_left = AdFeature(AdFeature_names=="to_left_eucZnorm");
       from_right = AdFeature(AdFeature_names=="from_right_eucZnorm");
       % [max_val,max_pos] = max(abs(subseqDiff_norm_i));
       [max_val,max_pos_rel] = max(abs(subseqDiff_norm_i(from_top1:to_top1)));
       max_pos = from_top1 + max_pos_rel -1;
       diff_center_i = abs(max_pos - subseqLen_i/2);
       topRange1 = (from_top1:to_top1);

       fluxSum_diff_i = sum(abs(diff(subseqDiff_norm_i(topRange1))));
       fluxsum_trg_i = sum(abs(diff(subseqTrg_norm_i(topRange1))));
       locMisFluxRatio_i  = fluxsum_trg_i/fluxSum_diff_i;
  
       if diff_center_i <= subseqLen_i*MaxDiffFromCenter_ratio ...
         && ( distLengthNorm_top1 < DistanceTop_Thr || same_sign(subseqTrg_norm_i,subseqNN_norm_i)) ...
         && distLengthNorm_top1*locMisFluxRatio_i < LocMisFluxRatio_Thr
           
         result_i = true;      
     
         if length(topRange1) >= subseqLen_i*MaxDiffFromCenter_ratio
           distLengthNorm_i  = distLengthNorm_calculate(subseqTrg_norm_i(topRange1),subseqNN_norm_i(topRange1));
         else
           from_12 = min(from_top1,from_top2);  to_12 = max(to_top1,to_top2);
           topRange12 = (from_12:to_12);  
           distLengthNorm_i  = distLengthNorm_calculate(subseqTrg_norm_i(topRange12),subseqNN_norm_i(topRange12));
         end
  
         % distLengthNorm_i = distZnormLengthNorm_calculate_top2(subseqTrg_i, subseqNN_i,from_top1,to_top1,from_top2,to_top2,to_left,from_right);
       else
         result_i = false;  
       end
       %% features
       PCSS_Feature = [lengthTest/SuppMedian, SuppMedian_ratio...
                          ;max(subseqDiff_norm_i),SubseqDiffMax_Thr...
                          ;boudaryRaw_absDiffLr_raio_i, BoudaryRaw_absDiffLr_raio_Thr...
                          ;dtwEuc_ratio_i, DtwEuc_ratio_Thr];
       Info_PCSS_i.values = PCSS_Feature; 
     else
              result_i = false;
     end
 
        % subseqDiffORG_i = (subseqTest_norm_i - subseqNN_norm_i); 
        % distLengthNormORG_i = norm(subseqDiffORG_i)/sqrt(length(subseqDiffORG_i)); 
        if  PLOT
          boundary_NN_i = consecutiveSeg_selected(ConsecutiveSeg_Names=="boundary");
          fromCS2_i = consecutiveSeg_selected(ConsecutiveSeg_Names=="from");
          toCS2_i = consecutiveSeg_selected(ConsecutiveSeg_Names=="to");

          Row = 2; Column = 2; 
          fig_id = figure(); T_test = (fromTest_i:toTest_i); T_train = (fromNN_i:toNN_i);
          
          subplot(Row,Column,1); 
          plot(T_test,D(T_test),'b');  xlim([T_test(1), T_test(end)]);
          title(['test[',num2str(fromTest_i),'-',num2str(toTest_i),'] len=',num2str(length(T_test))]);
          
          subplot(Row,Column,3); plot(T_train,D(T_train),'r'); xlim([min(fromNN_i,fromCS2_i) max(toNN_i,toCS2_i)]);
            T1 = (fromCS2_i:boundary_NN_i); T2 = (boundary_NN_i:toCS2_i);
            hold on; plot(T1,D(T1),'k','LineWidth',2); plot(T2,D(T2),'r','LineWidth',2);
                     xline(fromTest_i,'g'); xline(toTest_i,'g');
          title(['train[',num2str(fromCS2_i),'-',num2str(toCS2_i),'] len=',num2str(length(T_train))]);
         

           subplot(Row,Column,2);      
             plot((1:length(subseqNN_dtw_i)),subseqNN_dtw_i,'r'); hold on;
             plot((1:length(subseqTest_dtw_i)),subseqTest_dtw_i,'b'); 
             title(['len= ',num2str(length(subseqNN_dtw_i))]);
           
           subplot(Row,Column,4);      
             plot((1:length(subseqDiff_dtw_i)),subseqDiff_dtw_i); hold on; 
             xline(OverThrCenter_ref(1)); xline(OverThrCenter_ref(end)); xline(boundary_NN_ref_dtw,':');
             T = (fromAT_NN_dtw_ref:toAT_NN_dtw_ref); plot(T,subseqDiff_dtw_i(T),'r','LineWidth',2);
             title(num2str(distLengthNorm_i));
           % Fig_figtitle(char(Method_Attention)); 
          
          if DEBUG && isequal([fromTest_i,toTest_i],DEBUG_range)
           fprintf('[%s] debug end %d-%d\n',FuncName,fromTest_i,toTest_i);        
          end  
        end
  end
end

function SS = same_sign(D_trg_norm_i, D_NN_norm_i)
   MagInfo_trg = Spikelet_generateInitialSpikelet_ver_02(D_trg_norm_i,Inf); 
   MagInfo_nn = Spikelet_generateInitialSpikelet_ver_02(D_NN_norm_i,Inf); 

   [max_trg,maxPos_trg] = max(abs(MagInfo_trg.magnitude));
   sign_trg = sign(MagInfo_trg.magnitude(maxPos_trg));
   [max_nn,maxPos_nn] = max(abs(MagInfo_trg.magnitude));
   sign_nn = sign(MagInfo_nn.magnitude(maxPos_nn));
   if sign_trg == sign_nn
     SS = true;
   else
     SS = false;
   end
end

function dist_lenNorm = distLengthNorm_calculate(Subseq_1,Subseq_2)
   diff = Subseq_1 - Subseq_2;
   dist_lenNorm = norm(diff,2)/sqrt(length(diff));
end

function distLengthNorm = distZnormLengthNorm_calculate_top2(D1_Raw,D2_Raw,from_top1,to_top1,from_top2,to_top2,to_left,from_right)
   D1_Znorm = normalize(D1_Raw);
   D2_Znorm = normalize(D2_Raw);
   Ddiff_Znorm = D1_Znorm-D2_Znorm;

   Range_at_terminal = [from_top1:to_top1,from_top2:to_top2,1:to_left,from_right:length(Ddiff_Znorm)];
   Range_at_unique = unique(Range_at_terminal);
   Ddiff_Znorm_at_terminal = Ddiff_Znorm(Range_at_unique);
   dist_at = norm(Ddiff_Znorm_at_terminal,2);
   distLengthNorm_terminal = dist_at/sqrt(length(Ddiff_Znorm_at_terminal));

   Range_at_nonterminal = [from_top1:to_top1,from_top2:to_top2];
   Range_at_unique = unique(Range_at_nonterminal);
   Ddiff_Znorm_at_nonterminal = Ddiff_Znorm(Range_at_unique);
   dist_at = norm(Ddiff_Znorm_at_nonterminal,2);
   distLengthNorm_nonterminal = dist_at/sqrt(length(Ddiff_Znorm_at_nonterminal));

   distLengthNorm = max(distLengthNorm_terminal,distLengthNorm_nonterminal); 
end

function  SegIdList = include(Range,SegCenter)
  from = Range(1); to = Range(2);
  
  SegIdList = find(from < SegCenter & SegCenter < to);
end