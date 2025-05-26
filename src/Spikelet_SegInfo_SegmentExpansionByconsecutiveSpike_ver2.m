function SegInfo_CSE = Spikelet_SegInfo_SegmentExpansionByconsecutiveSpike_ver2(AdInfo,Names,SegInfo,Param)
  FuncName = 'Spikelet_SegInfo_SegmentExpansionByconsecutiveSpike';
  
  % SCoption = "SCcut_DP"; 
  % SCoption = "SCcut_DPandOP"; % PWcut_2
  % SCoption = "NoRestriction"; % PV4-1
  % If it is PV3, the opposite signs of the vertices are fewer than 2, and if it is PV4, the opposite signs of the vertices are fewer than 3.
  CsExpansion_ExcMethod = "SCcut_OP3"; 
  
  % ExpansionMethod = "first";   % PWcut2 
  % ExpansionMethod = "onlyPV4"; % PWcut2PV4
  ExpansionMethod = "PV3-4";     % PWcut2PV4
  % ExpansionMethod = "PV3-4_NotSpikeInSpike";
  %% argument
  FirstRate =  Param.ConsecutiveSpikeExpansion.ConsecutiveSpikeRate_1st;
  SecondRate = Param.ConsecutiveSpikeExpansion.ConsecutiveSpikeRate_2nd;
  SCstd = Param.ConsecutiveSpikeExpansion.SharpConvexityByStd;
  
  % AdInfo_names = ["id","from","to","type","magnitude","dist","distLengthNorm","from_NN","to_NN","spike_sign"];
  % SegDB_names = ["id","type","from","to","length","center","magnitude","fluctuation"];
  AnDB_org = AdInfo;
  MagInfo = SegInfo.MagInfo;
  D = MagInfo.data_org;
  SegDB_names = SegInfo.SegDB_names;
  last_id = SegInfo.last_id;
  Seg_exist = AnDB_org(:,[find(Names=="from"),find(Names=="to")]);

  %% preprocessing
  Index_spike = abs(AnDB_org(:,Names=="type")) == 2;
  AnDB_spike = AnDB_org(Index_spike,:); 

  [AnDB_val,AnDB_order] = sort(AnDB_spike(:,Names=="distLengthNorm"),'descend');
  AnDB = AnDB_spike(AnDB_order,:);

  %%
  Rank_1st = ceil(size(AnDB,1)*FirstRate);
  Rank_2nd = ceil(size(AnDB,1)*SecondRate);
  SegDB_new = []; id_new = last_id;
  for i = 1:Rank_1st
    from_i = AnDB(i,Names=="from"); 
    to_i = AnDB(i,Names=="to");
    spikeId_i = AnDB(i,Names=="id");
    DEBUG31 = true; % DEBUG31 = false;
    if DEBUG31 && (from_i == 23391 && to_i == 23412)
      fprintf('[%s] from_i= %d\n',FuncName,from_i);
    end
    
    type_i = AnDB(i,Names=="type");
    mag_i = AnDB(i,Names=="magnitude");
    RangeType_2nd = AnDB(i+1:Rank_2nd,[find(Names=="from"),find(Names=="to"),find(Names=="type")]);

    if ExpansionMethod == "first"
      Range_2nd = RangeType_2nd(:,[1,2]);
      [ov_rel,ForB,ov_type] = overlap_first([from_i,to_i],Range_2nd,MagInfo);     
    elseif ExpansionMethod == "onlyPV4"
      if abs(type_i) == 2
        Index_sameType = find(RangeType_2nd(:,3) == type_i);
        Range_2nd_sameType = RangeType_2nd(Index_sameType,[1,2]); 
        [ov_rel_PVfour,ForB,ov_type] = overlap_onlyPVFour([from_i,to_i],Range_2nd_sameType,MagInfo); 
        ov_rel = Index_sameType(ov_rel_PVfour);
      else
        ov_rel = []; ForB = []; ov_type = []; 
      end 
    elseif ExpansionMethod == "PV3-4"   
      if abs(type_i) == 2        
        center_i = spikeId_i; % If spike is abs(type_i == 2), then id_i equals to center_i.
        Index_sameType = find(RangeType_2nd(:,3) == type_i);
        Range_2nd_sameType = RangeType_2nd(Index_sameType,[1,2]); 
        [ov_rel_PV4_relative,ForB_4,ov_type_4] = overlap_point([from_i,to_i],Range_2nd_sameType,MagInfo); 
        ov_rel_PV4 = Index_sameType(ov_rel_PV4_relative);

        Index_diffType = find(RangeType_2nd(:,3) == (-1)*type_i);
        Range_2nd_diffType = RangeType_2nd(Index_diffType,[1,2]); 
        [ov_rel_PV3_relative,ForB_3,ov_type_3] = overlap_leg_overlap([from_i,center_i,to_i],Range_2nd_diffType,MagInfo); 
        % [ov_rel_PV3_relative,ForB_3,ov_type_3] = overlap_leg_point([from_i,center_i,to_i],Range_2nd_diffType,MagInfo); 
        ov_rel_PV3 = Index_diffType(ov_rel_PV3_relative);

        % if ~isempty(ov_rel_PV3);fprintf('debug\n');end

        ov_rel = [ov_rel_PV3;ov_rel_PV4]; ForB = [ForB_3; ForB_4]; ov_type = [ov_type_3;ov_type_4]; 
      else
        ov_rel = []; ForB = []; ov_type = []; 
      end  
    elseif ExpansionMethod == "PV3-4_NotSpikeInSpike"
      if abs(type_i) == 2        
        center_i = spikeId_i; % If spike is abs(type_i == 2), then id_i equals to center_i.
        Range_1st = [from_i,to_i];
        [Num,RangeList] = Spikelet_SegInfo_numInRange(SegInfo,Range_1st,type_i);
        if  Num == 2
          ov_rel = []; ForB = []; ov_type = [];   
        else
          % same type
          Index_sameType = find(RangeType_2nd(:,3) == type_i); 
          Range_2nd_sameType = RangeType_2nd(Index_sameType,[1,2]); 
          Range_2nd_sameType_index = true(size(Range_2nd_sameType,1),1);
          for j = 1:size(Range_2nd_sameType,1)
            Range_2nd_sameType_j = Range_2nd_sameType(j,:);
            Num_j = Spikelet_SegInfo_includeTwoSpike(Range_2nd_sameType_j,type_i,SegInfo);
            if Num_j >= 2
              Range_2nd_sameType_index(j) = false;
            end
          end
          Range_2nd_sameType = Range_2nd_sameType(Range_2nd_sameType_index,:);
          [ov_rel_PV4_relative,ForB_4,ov_type_4] = overlap_point([from_i,to_i],Range_2nd_sameType,MagInfo); 
          ov_rel_PV4 = Index_sameType(ov_rel_PV4_relative);

          % different type
          Index_diffType = find(RangeType_2nd(:,3) == (-1)*type_i);
          Range_2nd_diffType = RangeType_2nd(Index_diffType,[1,2]); 
          Range_2nd_diffType_index = true(size(Range_2nd_diffType,1),1);
          for j = 1:size(Range_2nd_diffType,1)
            Range_2nd_diffType_j = Range_2nd_diffType(j,:);
            Num_j = Spikelet_SegInfo_includeTwoSpike(Range_2nd_sameType_j,type_i,SegInfo);
            if Num_j >= 2
              Range_2nd_diffType_index(j) = false;
            end
          end
          Range_2nd_diffType = Range_2nd_diffType(Range_2nd_diffType_index,:);
          [ov_rel_PV3_relative,ForB_3,ov_type_3] = overlap_leg_overlap([from_i,center_i,to_i],Range_2nd_diffType,MagInfo); 
          % [ov_rel_PV3_relative,ForB_3,ov_type_3] = overlap_leg_point([from_i,center_i,to_i],Range_2nd_diffType,MagInfo); 
          ov_rel_PV3 = Index_diffType(ov_rel_PV3_relative);

          % if ~isempty(ov_rel_PV3);fprintf('debug\n');end

          ov_rel = [ov_rel_PV3;ov_rel_PV4]; ForB = [ForB_3; ForB_4]; ov_type = [ov_type_3;ov_type_4]; 
        end
      else
        ov_rel = []; ForB = []; ov_type = []; 
      end
    end
    for j = 1:length(ov_rel)
        ov_j = i + ov_rel(j);  
        ForB_j = ForB(j);
        ov_type_j = ov_type(j);

        from_ij = AnDB(ov_j,Names == "from");
        to_ij = AnDB(ov_j,Names == "to");
        type_ij = AnDB(ov_j,Names == "type");
        mag_ij = AnDB(ov_j,Names=="magnitude");
        spikeId_ij =  AnDB(ov_j,Names=="id");
        if ForB_j == "front" 
          sign_new = sign(type_ij);  
          from_new = from_ij;
          to_new = to_i;
        elseif ForB_j == "back"
          sign_new = sign(type_i);
          from_new = from_i;
          to_new = to_ij;
        end
        % T = [from_new:to_new]; figure(); plot(T,D(T));
        if ov_type_j == "leg"
          type_new = sign_new * (abs(type_i)+abs(type_ij)-1); 
        elseif ov_type_j == "point"
          type_new = sign_new * (abs(type_i)+abs(type_ij));
        end
        %% new definition: Magnitude refers to the one with the larger absolute value. The signs of type and magnitude do not necessarily match."
        [mag_new_abs, mag_pos] = max([abs(mag_i),abs(mag_ij)]);
        if mag_pos == 1
           mag_new  = mag_i;
        else
           mag_new  = mag_ij;
        end
        center_new = floor((from_i + to_i)/2); 
        fluctuation_new = max(D(from_i:to_i)) - min(D(from_i:to_i)); 
        length_new = to_new- from_new + 1;
        id_new = id_new + 1;
        seg_new = [id_new,type_new,from_new,to_new,length_new,center_new,mag_new,fluctuation_new,spikeId_i,spikeId_ij];
        if length(seg_new) ~= length(SegDB_names)
          fprintf('[%s] Column length is not consistent \n,',FuncName);
          return;
        end
        if sum(from_new == Seg_exist(:,1) & to_new == Seg_exist(:,2)) == 0
          D_new = D(from_new:to_new);
          if CsExpansion_ExcMethod == "SCcut_DP" % SCoption = "SCcut_DPandOP";
            [peak_rate,peak_num] = Spikelet_Subseq_doublePeak(D_new,SCstd);
            if peak_num <= 1
              SegDB_new = insert_SegDB(SegDB_new,SegDB_names,seg_new);
            end
          elseif CsExpansion_ExcMethod == "SCcut_DPandOP"
            [peak_rate,peak_num] = Spikelet_Subseq_doublePeak(D_new,SCstd);
            [peakOP_rate,peakOP_num] = Spikelet_Subseq_peakDiffOpposite(D_new,SCstd);
            if peak_num <= 1 && peakOP_num <= 2
              SegDB_new = insert_SegDB(SegDB_new,SegDB_names,seg_new);
            end
          elseif CsExpansion_ExcMethod == "SCcut_OP3"
            [peakOP_rate,peakOP_num] = Spikelet_Subseq_peakDiffOpposite(D_new,SCstd);
            if    (abs(type_new) == 3 && peakOP_num <= 2) ...
               || (abs(type_new) == 4 && peakOP_num <= 3) 
              SegDB_new = insert_SegDB(SegDB_new,SegDB_names,seg_new);
            end  
          elseif CsExpansion_ExcMethod == "NoRestriction" 
            SegDB_new = insert_SegDB(SegDB_new,SegDB_names,seg_new);  
          end
        else
          DEBUG_duplicate = true; DEBUG_duplicate = false;
          if DEBUG_duplicate
            fprintf('duplicate segment [%d %d]\n',from_new,to_new);
          end
        end
    end
  end

  SegInfo_CSE = SegInfo;
  SegInfo_CSE.SegDB_beforeCSE = SegInfo.SegDB;
  SegInfo_CSE.SegDB = SegDB_new;  
  if isempty(SegDB_new)
    SegInfo_CSE.last_id = SegInfo.last_id;
  else
    SegInfo_CSE.last_id = max(SegDB_new(:,SegDB_names=="id")); 
  end
end

function  [SegDB_new, SegDB_names] = insert_SegDB(SegDB,SegDB_names,seg_i)
  if isempty(SegDB) || size(SegDB,1) == 0
    SegDB_new = seg_i;
    return;
  end

  from_i = seg_i(SegDB_names == "from"); 
  to_i = seg_i(SegDB_names == "to");

  index_i = SegDB(:,SegDB_names == "from") == from_i & SegDB(:,SegDB_names == "to") == to_i;
  if sum(index_i) == 0
    SegDB_new = [SegDB; seg_i];
  else
    SegDB_new = SegDB;
  end
end

function [ov_rel,FrontOrBack,LegOrPoint] = overlap_leg_overlap(Range_1,Range_2nd,MagInfo) 
  from = Range_1(1);
  % center = Range_1(2);
  to = Range_1(3);
  range_1 = [from,to]; 
  mag_i = MagInfo.magnitude;

  ov_rel = []; FrontOrBack = []; LegOrPoint = [];
  for i = 1:size(Range_2nd,1)
   % if range_1(1) == 2833 && Range_2nd(i,1) == 2864; fprintf('debug\n'); end
    range_2 = Range_2nd(i,:);  
    [OVERLAP_i, overlap_length_i, FrontOrBack_i] = Spikelet_overlapSegment(range_1,range_2);
    if OVERLAP_i == true && (FrontOrBack_i == "front" || FrontOrBack_i == "back") 
      ov_rel_i = i;
      if overlap_length_i == 1
        LegOrPoint_i = "point";
      else
        LegOrPoint_i = "leg";
      end
      
      ov_rel = [ov_rel; ov_rel_i];
      LegOrPoint = [LegOrPoint; LegOrPoint_i];
      FrontOrBack = [FrontOrBack; FrontOrBack_i];
    end
  end
end 

function [ov_rel,FrontOrBack,LegOrPoint] = overlap_leg_point(Range_1,Range_2nd,MagInfo) 
  from = Range_1(1);
  center = Range_1(2);
  to = Range_1(3);
  mag_i = MagInfo.magnitude;

  ov_rel = []; FrontOrBack = []; LegOrPoint = [];

  for i = 1:size(Range_2nd,1)
    % if range_1(1) == 2833 && Range_2nd(i,1) == 2864; fprintf('debug\n'); end
    from2_i = Range_2nd(i,1); to2_i = Range_2nd(i,2);

    if center == to2_i || center == from2_i  
      ov_rel_i = i;
      if center == to2_i      
        FrontOrBack_i = "front";
      elseif center == from2_i
        FrontOrBack_i = "back";
      end
      LegOrPoint_i = "leg";
      
      ov_rel = [ov_rel; ov_rel_i];
      LegOrPoint = [LegOrPoint; LegOrPoint_i];
      FrontOrBack = [FrontOrBack; FrontOrBack_i];
    end
  end
end 

function  [ov_rel,FrontOrBack,LegOrPoint] = overlap_point(range_1,Range_2nd,MagInfo)
  ov_rel = []; FrontOrBack = []; LegOrPoint = [];

  %　mag_i = MagInfo.magnitude;
  for i = 1:size(Range_2nd,1)
    % if range_1(1) == 2833 && Range_2nd(i,1) == 2864; fprintf('debug\n'); end
    range_2 = Range_2nd(i,:);  

    if range_1(1) == range_2(2) || range_1(2) == range_2(1)  
      ov_rel_i = i;
      LegOrPoint_i = "point";
      if range_1(1) == range_2(2) 
          FrontOrBack_i = "front";  
      elseif  range_1(2) == range_2(1)  
          FrontOrBack_i = "back";
      end

      ov_rel = [ov_rel; ov_rel_i];
      LegOrPoint = [LegOrPoint; LegOrPoint_i];
      FrontOrBack = [FrontOrBack; FrontOrBack_i];
      
      DEBUG = true; DEBUG = false;
      if DEBUG
        D = MagInfo.data_org; 
        figure(); T = [min([range_1(1),range_2(1)]):max([range_1(2),range_2(2)])]; plot(T,D(T));
      end
    end
  end
end

function  [ov_rel,FrontOrBack,LegOrPoint] = overlap_onlyPVFour(range_1,Range_2nd,MagInfo)
  ov_rel = []; FrontOrBack = []; LegOrPoint = [];

  mag_i = MagInfo.magnitude;
  for i = 1:size(Range_2nd,1)
    % if range_1(1) == 2833 && Range_2nd(i,1) == 2864; fprintf('debug\n'); end
    range_2 = Range_2nd(i,:);  

    [OVERLAP, overlap_length, FrontOrBack_i] = Spikelet_overlapSegment(range_1,range_2);
    if OVERLAP == true && (FrontOrBack_i == "front" || FrontOrBack_i == "back") 
      ov_rel_i = i;
      if overlap_length == 1
        LegOrPoint_i = "point";
      else
        if FrontOrBack_i == "front" && mag_i(range_2(2)) ~= 0  
          LegOrPoint_i = "leg";
        elseif FrontOrBack_i == "back" && mag_i(range_2(1)) ~= 0 
          LegOrPoint_i = "leg";
        else
          LegOrPoint_i = "point";
        end
      end

      DEBUG = true; DEBUG = false;
      if DEBUG
        D = MagInfo.data_org; figure(); T = [min([range_1(1),range_2(1)]):max([range_1(2),range_2(2)])]; plot(T,D(T));
      end
      if LegOrPoint_i == "point"  % continuation only in case of point
        ov_rel = [ov_rel; ov_rel_i];
        LegOrPoint = [LegOrPoint; LegOrPoint_i];
        FrontOrBack = [FrontOrBack; FrontOrBack_i];
      end
    end
  end
end


function [ov_rel,FrontOrBack,LegOrPoint] = overlap_first(range_1,Range_2nd,MagInfo)
  ov_rel = []; FrontOrBack = []; LegOrPoint = [];
  mag_i = MagInfo.magnitude;
  for i = 1:size(Range_2nd,1)
    % if range_1(1) == 2833 && Range_2nd(i,1) == 2864; fprintf('debug\n'); end
    range_2 = Range_2nd(i,:);  
    [OVERLAP, overlap_length, FrontOrBack] = Spikelet_overlapSegment(range_1,range_2);
    if OVERLAP == true && (FrontOrBack == "front" || FrontOrBack == "back") 
      ov_rel = i;
      if overlap_length == 1
        LegOrPoint = "point";
      else
        if FrontOrBack == "front" && mag_i(range_2(2)) ~= 0  
          LegOrPoint = "leg";
        elseif FrontOrBack == "back" && mag_i(range_2(1)) ~= 0 
          LegOrPoint = "leg";
        else
          LegOrPoint = "point";
        end
      end

      DEBUG = true; DEBUG = false;
      if DEBUG
        D = MagInfo.data_org; figure(); T = [min([range_1(1),range_2(1)]):max([range_1(2),range_2(2)])]; plot(T,D(T));
      end
      return;
    end
  end
end
