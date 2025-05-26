function  [result_i,distLengthNorm_i,boundary_NN_i,fromAT_NN_i,toAT_NN_i,fromAT_test_i,toAT_test_i,Info_PCSS_i] ...
    = Spikelet_Pairwise_containTwoSpikesWithSameSign_ver_02(D,type_i,fromTest_i,toTest_i,fromNN_i,toNN_i,SegDB_2,Names,Param)
  %{
  Index_Seg2 = abs(SegDB(:,SegNames=="type")) == 2;
  SegDB_2 = SegDB(Index_Seg2,:);
  %}

  FuncName = 'Spikelet_Pairwise_containTwoSpikesWithSameSign_ver_02';
  fig_id = []; 
  DEBUG = true; % DEBUG = false;
  PLOT = true;  % PLOT = false;

  Info_PCSS_i = struct;
  Info_PCSS_i.names = ["lengthTest","max_subseqDiff","ComRange","BoudaryCloseness","AttentionVariation","boudaryRaw_absDiffLr_raio","DtwEuc_ratio"];
  Info_PCSS_i.values = [];

  if isempty(D)
    result_i = []; distLengthNorm_i = []; boundary_NN_i = []; fromAT_NN_i = []; toAT_NN_i=[]; fromAT_test_i=[]; toAT_test_i = []; 
    return;
  end
  %% Parameters
  ComRange_Thr = Param.LocMis_1.ComRange_Thr;
  BoundaryCloseness_Thr = Param.LocMis_1.BoundaryCloseness_Thr;
  SubseqDiffMax_Thr = Param.LocMis_1.SubseqDiffMax_Thr;
  SubseqDiffMin_Thr = Param.LocMis_1.SubseqDiffMin_Thr; 
  DiffFlux_ratio = Param.LocMis_1.DiffFlux_ratio; % Condition for the test subsequence to be relatively flat (using DTW)
  SuppMedian_ratio = Param.LocMis_1.SuppMedian_ratio;
  DtwEuc_ratio_Thr = Param.LocMis_1.DtwEuc_ratio;  
  BoudaryRaw_absDiffLr_raio_Thr = Param.LocMis_1.BoudaryRaw_absDiffLr_ratio; 

  % Method_Attention = "DtwZeroCross"; 

  %% default return values
  result_i = false; 
  distLengthNorm_i = NaN; boundary_NN_i = NaN; fromAT_NN_i = NaN; toAT_NN_i = NaN; fromAT_test_i = NaN; toAT_test_i = NaN; distLengthNormORG_i = NaN;

  %% debug
  %{
       53  2765, 2996
       57  15980, 16167
       66  7981 8034
       104 15980 16167
       161 2765 2996
       165 15980, 16167
       247 121735, 122113

  %}
  if DEBUG
    DEBUG_range = [5797, 5921];  % [2765,2996]; [2833,2942]; 152[22813, 23160]
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
  subseqTest_i = D(fromTest_i:toTest_i);
  subseqTest_norm_i = normalize(subseqTest_i);
  subseqDiff_norm_i = subseqTest_norm_i - subseqNN_norm_i; 
  subseqDiff_raw_i =  subseqTest_i - subseqNN_i;    
  Test_raw_fluc = max(subseqTest_i ) - min(subseqTest_i );
  
  % boudaryRaw_absDiffLr_raio_i = abs(abs(subseqDiff_raw_i(1)) - abs(subseqDiff_raw_i(end)))/Test_raw_fluc; 
  boudaryRaw_absDiffLr_raio_i = max(abs(subseqDiff_raw_i(1)), abs(subseqDiff_raw_i(end)))/Test_raw_fluc;
  
  %% main
  SegIdList = include([fromNN_i,toNN_i],SegCenter_sameSign);
  if length(SegIdList) >= 2 ...
     && lengthTest > SuppMedian*SuppMedian_ratio ...
     && max(subseqDiff_norm_i) >= SubseqDiffMax_Thr ...
     && boudaryRaw_absDiffLr_raio_i <= BoudaryRaw_absDiffLr_raio_Thr
    From = SegDB_2_sameSign(SegIdList,Names=="from");
    To = SegDB_2_sameSign(SegIdList,Names=="to");
    Center = SegDB_2_sameSign(SegIdList,Names=="center");
    ConsecutiveSeg = [];
    for i = 1:length(From)
      from_i = From(i);  
      to_i = To(i);
      center_i = Center(i);
      for j = 1:length(To)
        from_j = From(j);  
        to_j = To(j);
        center_j = Center(j);
        if to_i == from_j  
          RangeCom = common_range(from_i,to_j,fromNN_i,toNN_i);
          if ~isempty(RangeCom)
            boundary_ij = to_i;

            % comRangeMean
            length_com = RangeCom(2) - RangeCom(1) + 1;
            length_ij = to_j - from_i + 1;
            comRangeMean_ij = (length_com/length_ij + length_com/lengthNN_org)/2;
            comRange_ij = length_com/lengthNN_org;

            % firstSecondDiff_ratio_ij
            length_1st = to_i - from_i + 1;
            length_2nd = to_j - from_j + 1;
            firstSecondDiff_ratio_ij = max(length_1st,length_2nd)/min(length_1st,length_2nd);
            
            % boundaryClosenessToCenter
            boundaryClosenessToCenter_ij =  abs(boundary_ij - centerNN_org)/lengthNN_org;
            
            % ConsecutiveSeg
            ConsecutiveSeg_Names = ["from","to","boundary","center1","center2","comRange","comRangeMean","FirstSecondDiff_ratio","boundaryClosenessToCenter"];
            ConsecutiveSeg = [ConsecutiveSeg...
                             ; from_i, to_j, boundary_ij,center_i,center_j,comRange_ij,comRangeMean_ij,firstSecondDiff_ratio_ij,boundaryClosenessToCenter_ij];
          end
        end
      end
    end
    if ~isempty(ConsecutiveSeg)
      BoundaryCloseness = ConsecutiveSeg(:,ConsecutiveSeg_Names=="boundaryClosenessToCenter");
      ComRange  = ConsecutiveSeg(:,ConsecutiveSeg_Names=="comRange");  
      Candidate_pos = find(ComRange >= ComRange_Thr & BoundaryCloseness <=  BoundaryCloseness_Thr);
      if ~isempty(Candidate_pos)
        result_i = true;  
        
        [max_ComRange, max_posInCandidate_pos] = max(ComRange(Candidate_pos));
        selected_pos = Candidate_pos(max_posInCandidate_pos);
        consecutiveSeg_selected =  ConsecutiveSeg(selected_pos,:);

        min_BoundaryCloseness = min(BoundaryCloseness);
        % feature points
        boundary_NN_i = consecutiveSeg_selected(ConsecutiveSeg_Names=="boundary");
        fromCS2_i = consecutiveSeg_selected(ConsecutiveSeg_Names=="from");
        center1_i = consecutiveSeg_selected(ConsecutiveSeg_Names=="center1");
        toCS2_i = consecutiveSeg_selected(ConsecutiveSeg_Names=="to");
        center2_i = consecutiveSeg_selected(ConsecutiveSeg_Names=="center2");
        lengthCS2_i = toCS2_i - fromCS2_i+ 1;

        %% length threshold
        LeftRange_from = fromNN_i;
        RightRange_to = toNN_i;
        LeftRange_from_ref = LeftRange_from - fromNN_i+1;
        RightRange_to_ref  = RightRange_to -fromNN_i+1; 
        boundary_NN_ref = boundary_NN_i - fromNN_i + 1;

        % decide attention region ("DtwZeroCross")

           [dist_dtw,Path_NN, Path_Test] = dtw(subseqNN_norm_i, subseqTest_norm_i);
           subseqNN_dtw_i = subseqNN_norm_i(Path_NN);
           subseqTest_dtw_i = subseqTest_norm_i(Path_Test);
           subseqDiff_dtw_i = (subseqNN_dtw_i - subseqTest_dtw_i);

           boundary_dtw_index = find(Path_NN == boundary_NN_ref);
           boundary_NN_ref_dtw = round(mean(boundary_dtw_index));

           % find OverThrCenter_ref (over threshold points)
           range_start = max(round(boundary_NN_ref_dtw-0.25*length(subseqDiff_dtw_i)),1);  
           range_end =   min(round(boundary_NN_ref_dtw+0.25*length(subseqDiff_dtw_i)),length(subseqDiff_dtw_i));  
           Range = (range_start:range_end);
           
           OverThr_ref = find(abs(subseqDiff_dtw_i(Range)) >= SubseqDiffMin_Thr);
           if isempty(OverThr_ref)
             OverThr_ref = find(abs(subseqDiff_dtw_i(Range)) >= SubseqDiffMax_Thr);
           end
           OverThrCenter_ref = Range(OverThr_ref);

           if isempty(OverThr_ref)
             [max_val,max_pos_rel] = max(subseqDiff_dtw_i(Range));
             OverThrCenter_ref = range_start + max_pos_rel - 1;
           end
 
           % figure(); plot(subseqDiff_dtw_i); hold on; xline(OverThrCenter_ref(1)); xline(OverThrCenter_ref(end));
           OverThrCenter_Left_ref = OverThrCenter_ref(1);
           range_start_left = max(round(OverThrCenter_Left_ref-0.25*length(subseqDiff_dtw_i)),1); 
           LeftRange_ref = (range_start_left:OverThrCenter_Left_ref);
           zeroCrossings = findZeroCrossings(subseqDiff_dtw_i(LeftRange_ref));
           if isempty(zeroCrossings)
             fromAT_NN_dtw_ref = range_start;  
           else
             fromAT_NN_dtw_ref = LeftRange_ref(1)+ zeroCrossings(end) -1;
           end

           OverThrCenter_Right_ref = OverThrCenter_ref(end);
           range_end_right = min(round(OverThrCenter_Right_ref+0.25*length(subseqDiff_dtw_i)),length(subseqDiff_dtw_i));  
           RightRange_ref = (OverThrCenter_Right_ref:range_end_right);
           zeroCrossings = findZeroCrossings(subseqDiff_dtw_i(RightRange_ref));
           if isempty(zeroCrossings)
             toAT_NN_dtw_ref = range_end;  
           else
             toAT_NN_dtw_ref  = RightRange_ref(1) + zeroCrossings(1);
           end

          % difference at attention
          fromAT_NN_i = fromNN_i + Path_NN(fromAT_NN_dtw_ref) -1;
          fromAT_test_i = fromTest_i + Path_Test(fromAT_NN_dtw_ref);
          toAT_NN_i = fromNN_i + Path_NN(toAT_NN_dtw_ref) -1;
          toAT_test_i = fromTest_i + Path_Test(toAT_NN_dtw_ref);
          lengthAT_NN_i = toAT_NN_i - fromAT_NN_i + 1;
          lengthAT_test_i = toAT_test_i - fromAT_test_i + 1;
          lengthAt_mean_i = (lengthAT_NN_i + lengthAT_test_i)/2;
          lengthAt_max_i  = max(lengthAT_NN_i,lengthAT_test_i);
          RangeAT_dtw_rel = (fromAT_NN_dtw_ref:toAT_NN_dtw_ref); 

          subseqDiffAT_dtw_i = subseqDiff_dtw_i(RangeAT_dtw_rel); 
          subseqTestAT_dtw_i = subseqTest_dtw_i(RangeAT_dtw_rel); 

          variation_subseqDiffAt = sum(abs(diff(subseqDiffAT_dtw_i)));
          variation_subseqTestAt = sum(abs(diff(subseqTestAT_dtw_i)));

          if DEBUG % && isequal([fromTest_i,toTest_i],DEBUG_range) 
            % DEBUG_Feature_Names = ["lengthTest","max_subseqDiff)","ComRange","BoudaryCloseness","AttentionVariation"];
            % lengthTest/SuppMedian

          end
          if variation_subseqTestAt/variation_subseqDiffAt <= DiffFlux_ratio
            dtw_dist_i = sqrt(dtw(subseqTest_norm_i,subseqNN_norm_i,'squared'));
            %{
            [dtw_dist_builtin_i,Path1,Path2] = dtw(subseqTest_norm_i,subseqNN_norm_i,'squared');
            subseqTest_norm_DTW_i = subseqTest_norm_i(Path1);
            subseqNN_norm_DTW_i   = subseqNN_norm_i(Path2);
            subseqDiff_norm_DTW_i = subseqTest_norm_DTW_i - subseqNN_norm_DTW_i;
            dtw_dist_i = norm(subseqDiff_norm_DTW_i);
            %}
            euc_dist_i = norm(subseqDiff_norm_i);
            dtwEuc_ratio_i = dtw_dist_i / euc_dist_i;                 
            if dtwEuc_ratio_i >= DtwEuc_ratio_Thr 
              result_i = true;  
              % At (dtw)
              distLengthNormAt_i_dtw = norm(subseqDiffAT_dtw_i)/sqrt(lengthAt_max_i); 
              % distLengthNorm_i = norm(subseqDiffAT_i)/sqrt(lengthAt_mean_i); 
              % distLengthNorm_i = norm(subseqDiffAT_i)/sqrt(length(subseqDiffAT_i)); 
            
              % total (euc)
              subseqDiffTotal_i_euc = (subseqTest_norm_i - subseqNN_norm_i); 
              distLengthNorm_i_euc = norm(subseqDiffTotal_i_euc)/sqrt(length(subseqDiffTotal_i_euc));   
            
              % final result
              distLengthNorm_i = max(distLengthNormAt_i_dtw,distLengthNorm_i_euc);

              %% features
              PCSS_Feature = [lengthTest/SuppMedian, SuppMedian_ratio...
                          ;max(subseqDiff_norm_i),SubseqDiffMax_Thr...
                          ;max_ComRange,ComRange_Thr...
                          ;min_BoundaryCloseness,BoundaryCloseness_Thr...
                          ;variation_subseqTestAt/variation_subseqDiffAt, DiffFlux_ratio...
                          ;boudaryRaw_absDiffLr_raio_i, BoudaryRaw_absDiffLr_raio_Thr...
                          ;dtwEuc_ratio_i, DtwEuc_ratio_Thr];
              Info_PCSS_i.values = PCSS_Feature; 
            end
          else
            result_i = false;
          end

          DEBUG_local = true; DEBUG_local = false;
           if DEBUG_local
             figure(); plot((1:length(subseqDiff_dtw_i)),subseqDiff_dtw_i); hold on; 
             xline(OverThrCenter_ref(1)); xline(OverThrCenter_ref(end)); xline(boundary_NN_ref_dtw,':');
             T = (fromAT_NN_dtw_ref:toAT_NN_dtw_ref); plot(T,subseqDiff_dtw_i(T),'r','LineWidth',2);
             title(num2str(distLengthNorm_i));
             fprintf('debug local\n');
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
  end
end

function zeroCrossings = findZeroCrossings(Subseq)
    % Get the length of the signal
    n = length(Subseq);
    
    % Detect the indices of zero crossings
    zeroCrossings = [];
    for i = 1:n-1
        % If the product of consecutive elements is negative, a zero crossing occurs
        if Subseq(i) * Subseq(i+1) < 0
            zeroCrossings = [zeroCrossings, i];
        end
    end
end

function  SegIdList = include(Range,SegCenter)
  from = Range(1); to = Range(2);
  
  SegIdList = find(from < SegCenter & SegCenter < to);
end


function Range_common = common_range(from,to,from_NN,to_NN)
  from_max = max(from,from_NN);
  to_min = min(to,to_NN);
  if from_max > to_min
    Range_common = [];
  else
    Range_common =[from_max,to_min];
  end
end