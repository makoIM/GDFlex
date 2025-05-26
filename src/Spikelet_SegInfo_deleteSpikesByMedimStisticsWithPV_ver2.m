function [SegInfo,RevisedId,DeleteId]  = Spikelet_SegInfo_deleteSpikesByMedimStisticsWithPV_ver2(SegInfo,Param)
  FuncName = 'Spikelet_SegInfo_deleteSpikesByMedimStisticsWithPV_ver2';
  % SR_ByMagRatio = 0.75; SR_BySuppRatio = 0.75; % SR: Support Restriction
  
  SR_ByMagRatio = Param.SR_ByMagRatio;
  SR_BySuppRatio = Param.SR_BySuppRatio;
  MmedianSupp = Param.MmedianSupp;
  MmedianSupp_upperLimit = Param.MmedianSupp_upperLimit;   

  MagInfo = SegInfo.MagInfo;  D = SegInfo.data_org;
  SegDB = SegInfo.SegDB;  SegDB_org = SegDB;
  Names = SegInfo.SegDB_names;

  %% positive spikes
  SegDB_pos = SegDB(SegDB(:,Names=="type") == 2,:);  SegDB_pos_org = SegDB_pos;
  [SegDB_pos,RevisedId_pos,DeleteId_pos] = deleteSpikesByMediamStatisticsWithPV(SegDB_pos,Names,MmedianSupp,MmedianSupp_upperLimit,MagInfo,SR_ByMagRatio, SR_BySuppRatio);

  SegDB_neg = SegDB(SegDB(:,Names=="type") == -2,:);  SegDB_neg_org = SegDB_neg;
  [SegDB_neg,RevisedId_neg,DeleteId_neg] = deleteSpikesByMediamStatisticsWithPV(SegDB_neg,Names,MmedianSupp,MmedianSupp_upperLimit,MagInfo,SR_ByMagRatio, SR_BySuppRatio);

  SegDB = [SegDB_pos;SegDB_neg];
  RevisedId = [RevisedId_pos; RevisedId_neg];
  DeleteId = [DeleteId_pos; DeleteId_neg];
  
  % output
  SegInfo.SegDB = SegDB;
  SegInfo.last_id = max(length(D),max(SegDB(:,Names=="id")));
end
function [SegDB_new,RevisedId, DeleteId] = deleteSpikesByMediamStatisticsWithPV(SegDB,Names,MmedianSupp,MmedianSupp_upperLimit,MagInfo,SR_MagRatio,SR_SuppRatio)
  FuncName = 'Spikelet_SegInfo_deleteSpikesByMedimStisticsWithPV/deleteSpikesByMediamStatisticsWithPV';
  SegDB_new = SegDB;
  Stat = statistics_SegDB(SegDB,Names);
  SuppLen = Stat.suppLen;
  Median = Stat.suppMedian;
  OverSupportOrderId = find(Median * MmedianSupp < SuppLen);
  Mag = MagInfo.magnitude; 
  D = MagInfo.data_org;
  Mag = SegDB(:, Names=="magnitude");
  MagThres_pos = min(Mag(Mag > 0));
  MagThres_neg = max(Mag(Mag < 0));

  RevisedId = []; DeleteId = [];DeleteOrder = [];
  for i = 1:length(OverSupportOrderId)
    order_i = OverSupportOrderId(i);
    id_i = SegDB(order_i,Names=="id");
    mag_i = SegDB(order_i,Names=="magnitude");
    magSign_i = sign(mag_i);
    from_i = SegDB(order_i,Names=="from");
    center_i = SegDB(order_i,Names=="center");
    to_i   = SegDB(order_i,Names=="to");
    length_i = to_i - from_i + 1;
    RangeLeft = (from_i:center_i-1)';
    RangeRight = (center_i+1:to_i)';
    % RangeWithoutCenter = [RangeLeft; RangeRight];
    % MagRange_i = Mag(RangeWithoutCenter);

    DEBUG = true; DEBUG = false;
    if DEBUG && from_i == 14082   % from_i == 14082, 6237
      fprintf('debug [%s]\n',FuncName);
    end
    %%
    overLeft_rel  = find( abs(D(RangeLeft) - D(center_i)) > abs(mag_i)*SR_MagRatio,1,"last");
    overRight_rel = find( abs(D(RangeRight)- D(center_i)) > abs(mag_i)*SR_MagRatio,1,"first");
    if isempty(overLeft_rel) | isempty(overRight_rel)
      absMag_left = max(abs(D(RangeLeft) - D(center_i)));
      absMag_right = max(abs(D(RangeRight)- D(center_i)));
      mag_i_revised = min([absMag_left,absMag_right]);
      overLeft_rel  = find( abs(D(RangeLeft) - D(center_i)) > abs(mag_i_revised)*SR_MagRatio,1,"last");
      overRight_rel = find( abs(D(RangeRight)- D(center_i)) > abs(mag_i_revised)*SR_MagRatio,1,"first");
    end
    overLeft = RangeLeft(1) + overLeft_rel - 1;
    overRight = RangeRight(1) + overRight_rel - 1;

    if ((overRight - overLeft) < length_i*SR_SuppRatio ) ...
       || length_i >= Median * MmedianSupp_upperLimit 
      from_new_i = monotone_extention_left(D,RangeLeft,overLeft);
      to_new_i = monotone_extention_right(D,RangeRight,overRight);

      if DEBUG && from_i == 6237 % from_i == 14082
        fprintf('debug [%s]\n',FuncName);
      end

      % revise Segment
      if from_new_i == from_i && to_new_i == to_i 
        MagInfo_ini = MagInfo.output.generateInitialSpikelet;
        Mag_i_centerExcluded = MagInfo_ini.magnitude;
        Mag_i_centerExcluded(center_i) = 0;
        Mag_i_centerExcluded = Mag_i_centerExcluded(from_i+1:to_i-1);

        type_i = sign(mag_i);
        if type_i > 0
          magThred = MagThres_pos;
        else
          magThred = abs(MagThres_neg);
        end
        if max(abs(Mag_i_centerExcluded)) < magThred
          DeleteId =  [DeleteId; id_i];
          DeleteOrder = [DeleteOrder, order_i];
        end
      else
        SegDB_new(order_i,Names== "from") = from_new_i;
        SegDB_new(order_i,Names== "to")   = to_new_i;
        SegDB_new(order_i,Names== "length") = to_new_i - from_new_i + 1;
        magAbs_left_i = max(abs(D(center_i)-D(RangeLeft)));
        magAbs_right_i = max(abs(D(center_i)-D(RangeRight)));
        SegDB_new(order_i,Names== "magnitude") = magSign_i * min(magAbs_left_i,magAbs_right_i);
        RevisedId = [RevisedId; id_i];
      end
      % debug
      DEBUG = true; DEBUG = false;
      if DEBUG
        figure(); T_i = (from_i:to_i);
        plot(T_i,D(T_i)); hold on;
        xline(overLeft,'--'); xline(from_new_i);
        xline(overRight,'--'); xline(to_new_i);
      end
    end
  end

  if ~isempty(DeleteOrder)
    SegDB_order_remained = setdiff((1:size(SegDB_new,1)),DeleteOrder);
    SegDB_new = SegDB_new(SegDB_order_remained,:);
  end
end

function cur = monotone_extention_left(D,RangeLeft,overLeft)
   D_center = D(RangeLeft(end)); 
   D_overLeft = D(overLeft);
   Sign = sign(D_center - D_overLeft);
   cur = overLeft;
   for t = overLeft-1:-1:RangeLeft(1)
     if sign(D(t+1) - D(t)) == Sign  % Can D(t+1) be extended to D(t)?
       cur = t;  
     else
       break;
     end
   end
end

function cur = monotone_extention_right(D,RangeRight,overRight)
   D_center = D(RangeRight(1)); 
   D_overRight = D(overRight);
   Sign = sign(D_center - D_overRight);
   cur = overRight;
   for t = overRight+1:RangeRight(end)
     if sign(D(t-1) - D(t)) == Sign  % Can D(t-1) be extended to D(t)?
       cur = t;
     else
       break;  
     end
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
  MagPos = Info.magnitude(Info.magnitude > 0);
  MagNeg = Info.magnitude(Info.magnitude < 0);
  Info.magPosMax = max(MagPos); 
  Info.magPosMin = min(MagPos);
  Info.magNegMax = max(MagNeg);
  Info.magNegMin = min(MagNeg);
  Info.magStd = std(Info.magnitude);
  Info.magMedian = median(Info.magnitude);
end

 %{
    RangeLeft = (from_i:center_i-1)';
    RangeRight = (center_i+1:to_i)';

    MagThres = SR_MagRatio*abs(mag_i);
    overPosLeft = find(abs(Mag(RangeLeft)) >= MagThres & sign(Mag(RangeLeft))==magSign_i,1,"last");
    overPosRight = find(abs(Mag(RangeRight)) >= MagThres & sign(Mag(RangeRight))==magSign_i,1,"first");

    if ~isempty(overPosLeft)
      leftBoundary = from_i + overPosLeft -1;  % RangeLeft = (from_i:center_i-1)';
      RangeLeftBoundary = (leftBoundary:center_i);
      [maxVal, maxPos] = max(abs(D(RangeLeftBoundary) - D(center_i)));
      from_new_i = leftBoundary+maxPos-1;
      SegDB_new(order_i,Names=="from") = from_new_i;
    end

    if ~isempty(overPosRight)
      rightBoundary = center_i + overPosRight;   % RangeRight = (center_i+1:to_i)';
      RangeRightBoundary = (center_i:rightBoundary);
      [maxVal, maxPos] = max(abs(D(RangeRightBoundary) - D(center_i)));
      to_new_i = center_i+maxPos-1;
      SegDB_new(order_i,Names=="to") = to_new_i;
    end  
    %}
