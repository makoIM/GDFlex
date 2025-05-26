function [SegInfo,Revised_id] = Spikelet_SegInfo_restrictSupportBySpikeInSpike_ver2(SegInfo,Param)
  %{
   constatint: all the spikes have the same sign. 

   For the spikes included in SpikeIdTarget, 
   exclude regions where the support of spikes with the same sign from SpikeIdAll, 
   including the vertex, results in a continuous area smaller than mediam of the support length of spikes with the same sign
  %}

  %% parameters
  if ~isfield(Param,"MagSumRate")
    Param.MagSumRate = 1;
  end
  if ~isfield(Param,"MagSingleRate")
    Param.MagSingleRate = 0.75;
  end

  %% argument
  MagInfo = SegInfo.MagInfo;
  SegDB = SegInfo.SegDB;
  Names = SegInfo.SegDB_names;
  
 
  SegDB_pos = SegDB(SegDB(:,Names=="type") == 2, :);
  SegDB_neg = SegDB(SegDB(:,Names=="type") == -2,:);
  
  last_id = SegInfo.last_id;
  [SegDB_pos, Revised_id_pos, last_id,SegId_add_pos] = restrictSupportBySpikeInSpike(SegDB_pos,Names,MagInfo,last_id,Param);
  [SegDB_neg, Revised_id_neg, last_id,SegId_add_neg] = restrictSupportBySpikeInSpike(SegDB_neg,Names,MagInfo,last_id,Param);

  % output
  SegInfo.SegDB = [SegDB_pos;SegDB_neg];
  SegInfo.last_id = last_id;
  Revised_id = [Revised_id_pos; Revised_id_neg];
end

function  [SegDB, RevisedSegId, last_id,SegId_add] = restrictSupportBySpikeInSpike(SegDB,Names,MagInfo,last_id,Param)
   MagSumRate = Param.MagSumRate;

   D = MagInfo.data_org;
   SegDB_org = SegDB;
   MagAbs = abs(SegDB(:,Names=="magnitude"));
   [MagAbs_soreted_val,MagAbs_sorted_index] = sort(MagAbs,"descend");
   SegDB = SegDB(MagAbs_sorted_index,:);
   MagAbsMin = min(MagAbs);
   SuppLen = SegDB(:,Names=="length");
   SuppMedian = median(SuppLen);

   RevisedSegId = []; SegId_add = []; SegOrderId_Flag = true(size(SegDB,1),1); 
   for i = 1:length(MagAbs_sorted_index)
     id_i = MagAbs_sorted_index(i); 
     center_i = SegDB(i,Names=="center"); 
     mag_i    = SegDB(i,Names=="magnitude");
     from_i   = SegDB(i,Names=="from"); 
     to_i     = SegDB(i,Names=="to");
     length_i = to_i-from_i+1;
     
     DEBUG_PINPOINT = true; FROM_DEBUG = 97179; DEBUG_PINPOINT = false; % 84012
     if DEBUG_PINPOINT && from_i == FROM_DEBUG
        if exist('figDebug_id','var'); close(figDebug_id); end
        figDebug_id = figure(); T=(from_i:to_i); plot(T,D(T),'b'); 
        Title_debug = [num2str(from_i),'-',num2str(center_i),'-',num2str(to_i),' mag=',num2str(mag_i),' len=',num2str(to_i-from_i+1)]; 
        title(Title_debug); 
        LIM_X = lim_margin(T,0.1); xlim(LIM_X);

        SpikeCenterAll = SegDB(:,Names=="center");
        InSpike = SpikeCenterAll(from_i < SpikeCenterAll & SpikeCenterAll < to_i);
        for ins = 1:length(InSpike)
          if ins == 1; hold on; end;
          center_ins = InSpike(ins);
          if center_ins ~= center_i
            T_ins = (MagInfo.left(center_ins):MagInfo.right(center_ins));
            plot(T_ins,D(T_ins),'r');
          end
        end
        fprintf('InSpike: ');
        disp(setdiff(InSpike,center_i));
        fprintf('\n');
     end

      %% first boundary
      % right first boudary
      DEBUG_69 = true;DEBUG_69 = false;
      if DEBUG_69 && (center_i == 29342 || i == 364)
        fprintf("DEBUG_320\n");  
      end

      [LeftBoudary_i,RightBoudary_i] = find_boundary(SegDB,Names,center_i,D,Param);
      
      if DEBUG_PINPOINT && from_i == FROM_DEBUG && (RightBoudary_i ~= to_i || LeftBoudary_i ~= from_i)
        hold on; xline(LeftBoudary_i); xline(RightBoudary_i);
        Title_debug = [Title_debug,'->',num2str(RightBoudary_i - LeftBoudary_i +1)];
        title(Title_debug);
        fprintf('[left_before-after] %d %d\n',from_i,LeftBoudary_i);
        fprintf('[right_before-after] %d %d\n',to_i,RightBoudary_i);
      end
      
      %% output for each i
      % check SpikeIdTarget_Flag
      if LeftBoudary_i ~= from_i || RightBoudary_i ~= to_i
        T_i = (LeftBoudary_i:RightBoudary_i);
        revised_absMag = max(D(T_i)) - min(D(T_i));
        if revised_absMag < MagAbsMin  
          RevisedSegId = [RevisedSegId;id_i]; 
          SegOrderId_Flag(i) = false; 
        end
      end
      if SegOrderId_Flag(i) && LeftBoudary_i ~= from_i
         SegDB(i,Names=="from") = LeftBoudary_i; 
      end
      if SegOrderId_Flag(i) && RightBoudary_i ~= to_i
         SegDB(i,Names=="to") = RightBoudary_i; 
      end
      if SegOrderId_Flag(id_i) && (LeftBoudary_i ~= from_i || SegOrderId_Flag(id_i) && RightBoudary_i ~= to_i)
         left_mag_abs =  abs(D(center_i) - D(SegDB(id_i,Names=="from")));
         right_mag_abs = abs(D(center_i) - D(SegDB(id_i,Names=="to")));
         SegDB(i,Names=="magnitude") = sign(mag_i)*min(left_mag_abs,right_mag_abs);
         SegDB(i,Names=="length") = SegDB(i,Names=="to") - SegDB(i,Names=="from") + 1;
         RevisedSegId = [RevisedSegId; id_i ];

         DEBUG_RevisedSeg = true; DEBUG_RevisedSeg = false;
         if DEBUG_RevisedSeg
           figDebug_id = figure(); T=(from_i:to_i); plot(T,D(T)); 
           hold on; xline(LeftBoudary_i); xline(RightBoudary_i);
           Title_debug = [num2str(from_i),'-',num2str(center_i),'-',num2str(to_i),' mag=',num2str(mag_i),' len=',num2str(to_i-from_i+1)]; 
           title(Title_debug); 
         end
      end

      DEBUG_117 = false;
      if DEBUG_117 && (LeftBoudary_i ~= from_i || RightBoudary_i ~= to_i)
        if exist('figDebug_id','var'); close(figDebug_id); end  
        figDebug_id = figure(); T=(from_i:to_i); plot(T,D(T)); 
        Title_debug = [num2str(from_i),'-',num2str(center_i),'-',num2str(to_i),' mag=',num2str(mag_i),' len=',num2str(to_i-from_i+1)]; 
        LIM_X = lim_margin(T,0.1); xlim(LIM_X);

        hold on; xline(LeftBoudary_i); xline(RightBoudary_i);
        Title_debug = [Title_debug,'->',num2str(RightBoudary_i - LeftBoudary_i +1)];
        title(Title_debug);
        fprintf('[left_before-after] %d %d\n',from_i,LeftBoudary_i);
        fprintf('[right_before-after] %d %d\n',to_i,RightBoudary_i);

        CenterAll = SegDB(:,Names=="center");
        InSpike = CenterAll(from_i < CenterAll & CenterAll < to_i);
        fprintf('InSpike: ');
        disp(setdiff(InSpike,center_i));
        fprintf('\n');
      end
   end
  SegDB = SegDB(SegOrderId_Flag,:);
end

% return value: The right terminal in the segment segments with a length of TMPS or more, or the boundary point
function [from_i_new, to_i_new] = find_boundary(SegDB,Names,center_i,D,Param)
  %% threshold
  MagSumRate = Param.MagSumRate;
  MagSingleRate = Param.MagSingleRate;

  %% atribute
  index_i = SegDB(:,Names=="center") == center_i;
  from_i = SegDB(index_i,Names=="from");
  to_i = SegDB(index_i,Names=="to");
  mag_i = SegDB(index_i,Names=="magnitude");
  sign_i = sign(mag_i);

  CenterIncluded_order_i = find(from_i < SegDB(:,Names == "center") ...
                     & SegDB(:,Names == "center") < to_i ...
                     & sign(SegDB(:,Names == "magnitude")) == sign_i ...
                     & SegDB(:,Names == "center") ~= center_i );

  magAbsSum = abs(sum(SegDB(CenterIncluded_order_i,Names=="magnitude")));
  if isempty(CenterIncluded_order_i) ...
     || magAbsSum <= abs(mag_i)*MagSumRate ...
     || ( length(CenterIncluded_order_i) == 1 && magAbsSum >= abs(mag_i)*MagSingleRate ) 
    from_i_new = from_i;
    to_i_new = to_i;  
  else
    % find left  
    CenterIncluded_order_i_left_index = SegDB(CenterIncluded_order_i,Names=="center") < center_i;
    if sum(CenterIncluded_order_i_left_index) == 0
      from_i_new = from_i;  
    else
      CenterIncluded_order_i_left = CenterIncluded_order_i(CenterIncluded_order_i_left_index);
      from_i_new = find_left(SegDB,Names,center_i,D,CenterIncluded_order_i_left);
    end
    % find right
    CenterIncluded_order_i_right_index = SegDB(CenterIncluded_order_i,Names=="center") > center_i;
    if sum(CenterIncluded_order_i_right_index ) == 0
      to_i_new = to_i;  
    else
      CenterIncluded_order_i_right = CenterIncluded_order_i(CenterIncluded_order_i_right_index);
      to_i_new = find_right(SegDB,Names,center_i,D,CenterIncluded_order_i_right);
    end
  end
end

function to_i_new = find_right(SegDB,Names,center_i,D,CenterIncluded_order_i_right)
   [MinVal,MinleftMostPos] = min(SegDB(CenterIncluded_order_i_right,Names=="center"));
   boundaryRight_order_i = CenterIncluded_order_i_right(MinleftMostPos);
   boundary_time_i = SegDB(boundaryRight_order_i, Names=="from");
   Range_i = (center_i:boundary_time_i);
   [maxVal, maxPos] = max(abs(D(Range_i)- D(center_i)));
   to_i_new = center_i + maxPos -1;
end

function from_i_new = find_left(SegDB,Names,center_i,D,CenterIncluded_order_i_left)
    [MaxVal,MaxRightMostPos] = max(SegDB(CenterIncluded_order_i_left,Names=="center"));
    boundary_order_i = CenterIncluded_order_i_left(MaxRightMostPos);
    boundary_time_i = SegDB(boundary_order_i, Names=="to");
    Range_i = (boundary_time_i:center_i);
    [maxVal, maxPos] = max(abs(D(Range_i)- D(center_i)));
    from_i_new = boundary_time_i + maxPos -1;  
end


function Supp_independ_i = supp_independent(SegDB,Names,D,order_i)
  % support of order_i
  Supp_i = false(length(D),1);
  from_i = SegDB(order_i,Names=="from");  to_i = SegDB(order_i,Names=="to"); 
  Supp_i(from_i:to_i) = true(to_i-from_i+1,1);

  % support other than order_i
  Supp_OtherThan_i = false(length(D),1);
  for order_j = 1:size(SegDB,1)
    if order_j ~= order_i  
     from_j = SegDB(order_j,Names=="from"); to_j = SegDB(order_j,Names=="to");  
     Supp_OtherThan_i(from_j:to_j) = true(to_j-from_j+1,1);

     DEBUG = true;
     if DEBUG
       Supp_j = false(length(D),1); Supp_i(from_j:to_j) = true(to_j-from_j+1,1);
     end
    end
  end

  % support of order_i which is not included by other segments
  Supp_independ_i = Supp_i&~Supp_OtherThan_i;
  if sum(Supp_independ_i(from_i+1:to_i-1)) == (to_i-1) - (from_i+1) + 1
    Supp_independ_i(from_i) = true;
    Supp_independ_i(to_i) = true;
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
