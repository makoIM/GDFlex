function [PvNumOverCut,SegInfo,AnDB_new_i,SegDB_new_i,Feature,Feature_names,MagInfo_cand] = Spikelet_AnomalyInfoInstance_extractSegmentForPVover5(AnInst,An_names,D,Command,SegInfo,training_region,Param)
   %% parameters     
   ExtractMode = true; ExtractMode = false; 
   MagRatio_Thr = 0.5; %  MagRatio_Thr = 0.6;
   %% output default
   AnDB_new_i = []; SegDB_new_i = [];

   %% Anomaly Instance
   spikeId = AnInst(An_names=="id");
   % if spikeId == 8622;fprintf("debug\n");end
   from = AnInst(An_names=="from");
   to = AnInst(An_names=="to");      
   type = AnInst(An_names=="type");
   dist_original = AnInst(An_names=="distLengthNorm");
   %% spikelet decomposition
   ParamAprx = command2ParamAprx(Command);
   Range_cand = (from:to)';
   D_cand = D(Range_cand);
   MagInfo_cand = Spikelet_aproximation_ver_05(D_cand,ParamAprx);
   Mag_cand = MagInfo_cand.magnitude;
   [Mag_descend,Pos_descend] = sort(Mag_cand,'descend');
   [Mag_ascend,Pos_ascend] = sort(Mag_cand,'ascend');

   %% mag ratio
   PvNumOverCut = false;
   magRatio = 0; type_deemed = type;  pvNum = abs(type);
   if pvNum == 4
        if type == 4 && Mag_descend(3) > 0
          magRatio12 = Mag_descend(2)/Mag_descend(1);
          magRatio23 = Mag_descend(3)/Mag_descend(2);
          MagPos = [Pos_descend(1:3);Pos_ascend(1:2)];
          MagVal = [Mag_descend(1:3);Mag_ascend(1:2)];
          magRatio = min([magRatio12,magRatio23]);
        elseif type == - 4 && Mag_ascend(3) < 0
          magRatio12 = Mag_ascend(2)/Mag_ascend(1);  
          magRatio23 = Mag_ascend(3)/Mag_ascend(2);
          MagPos = [Pos_ascend(1:3);Pos_descend(1:2)];
          MagVal = [Mag_ascend(1:3);Mag_descend(1:3)];
          magRatio = min([magRatio12,magRatio23]);
        end
        
        if magRatio > MagRatio_Thr
          PvNumOverCut = true;  
          type_deemed = sign(type)*6;
          if ExtractMode
            [SegDB_new_i, segLast_id] = extractLRsegment([from,to,spikeId],MagPos,MagVal,type,D,SegInfo,pvNum);
          end
          % MagPoint = [1,,length(Range_cand_i)];
        end
      elseif pvNum == 3 
        if Mag_descend(2) > 0 && Mag_ascend(2) < 0
           mag_12_n = Mag_ascend(2)/Mag_ascend(1);
           mag_12_p = Mag_descend(2)/Mag_descend(1);
           magRatio = min(mag_12_n,mag_12_p);
           if magRatio > MagRatio_Thr
             PvNumOverCut = true;   
             type_deemed = sign(type)*5;
             if ExtractMode
               MagPos = [Pos_ascend(1:2);Pos_descend(1:2)];
               MagVal = [Mag_ascend(1:2);Mag_descend(1:2)];
               [SegDB_new_i, segLast_id] = extractLRsegment([from,to,spikeId],MagPos,MagVal,type,D,SegInfo,pvNum);
             end
           end
        end
   end

   Feature = [dist_original,magRatio,type_deemed];
   Feature_names = ["dist_org","magRatio","type_deemed"];

   if ~isempty(SegDB_new_i)
     SegInfo.SegDB = [SegInfo.SegDB; SegDB_new_i];
     SegInfo.last_id = max(SegDB_new_i(:,SegInfo.SegDB_names=="id"));

     SegInfo_new = SegInfo;
     SegInfo_new.SegDB = SegDB_new_i;
     [AnDB_new_i,AdDB_names] = Spikelet_AD_TrTe_STOMP_ver_02(SegInfo_new,training_region,Param);
   end
end

function [SegDB_new, segLast_id] = extractLRsegment(range_origin,MagPos,MagVal,type_org,D,SegInfo,pvNum)
   %% output
   SegDB_new = [];
   %% SegInfo
   SegDB = SegInfo.SegDB;
   SegDB_names = SegInfo.SegDB_names;
   segLast_id = SegInfo.last_id;

   from_origin = range_origin(1);
   to_origin = range_origin(2);
   spikeId_origin = range_origin(3);
   index_spikeId = SegInfo.SegDB(:,SegInfo.SegDB_names == "id") == spikeId_origin;
   component_1 = SegInfo.SegDB(index_spikeId,SegInfo.SegDB_names == "component_1");
   component_2 = SegInfo.SegDB(index_spikeId,SegInfo.SegDB_names == "component_2");
   length_cand = to_origin - from_origin + 1;

   %% extract new segment
   [MagPos_sorted,MagPos_soretedOrder] = sort(MagPos,"ascend");
   MagSeq_val = MagVal(MagPos_soretedOrder);
   MagSeq_time = [1;MagPos_sorted;length_cand];
   
   %% left 
   from_left = MagSeq_time(1);
   to_left = MagSeq_time(pvNum+1);
   Mag_left =  MagSeq_val(1:pvNum-1);

   % SegDB_names = ["id","type","from","to","length","center","magnitude","fluctuation"];
   Exist = checkIfExist(from_left,to_left,SegDB,SegDB_names);
   if ~Exist
      [seg_i,segLast_id] = create_segment(D,from_origin,from_left,to_left,type_org,MagVal,MagPos,Mag_left,segLast_id,component_1,component_2);  
      SegDB_new = [SegDB_new; seg_i];
   end

   if pvNum == 3
     from_left = MagSeq_time(1);
     to_left = MagSeq_time(pvNum+1+1);
     Mag_left =  MagSeq_val(1:pvNum-1+1);

     % SegDB_names = ["id","type","from","to","length","center","magnitude","fluctuation"];
     Exist = checkIfExist(from_left,to_left,SegDB,SegDB_names);
     if ~Exist
      [seg_i,segLast_id] = create_segment(D,from_origin,from_left,to_left,type_org,MagVal,MagPos,Mag_left,segLast_id,component_1,component_2);  
      SegDB_new = [SegDB_new; seg_i];
     end
   end
   %% right
   from_right = MagSeq_time(end-pvNum);
   to_right = MagSeq_time(end);
   Mag_right = MagSeq_val(end-pvNum+1:end);
   Exist = checkIfExist(from_left,to_left,SegDB,SegDB_names);
   if ~Exist
     [seg_i,segLast_id] = create_segment(D,from_origin,from_right,to_right,type_org,MagVal,MagPos,Mag_right,segLast_id,component_1,component_2);  
     SegDB_new = [SegDB_new; seg_i];
   end 
   if pvNum == 3
     from_right = MagSeq_time(end-pvNum-1);
     to_right = MagSeq_time(end);
     Mag_right = MagSeq_val(end-pvNum+1-1:end);
     Exist = checkIfExist(from_left,to_left,SegDB,SegDB_names);
     if ~Exist
      [seg_i,segLast_id] = create_segment(D,from_origin,from_right,to_right,type_org,MagVal,MagPos,Mag_right,segLast_id,component_1,component_2);  
      SegDB_new = [SegDB_new; seg_i];
     end 
   end
end

function [seg_i,segLast_id] = create_segment(D,from_origin,from_left,to_left,type_org,MagVal,MagPos,Mag_left,segLast_id,component_1,component_2)
     id_i = segLast_id + 1; segLast_id = segLast_id + 1;
     type_i = sign(type_org)*3;
     from_i = from_left + from_origin - 1;
     to_i = to_left + from_origin - 1;
     length_i = to_i - from_i + 1;
     [maxVal, argMaxVal_pos] = max(abs(MagVal));
     center_i = MagPos(argMaxVal_pos) + from_origin -1;
     magnitude_i = sign(type_i)*min(abs(Mag_left));
     fluctuation_i = max(D(from_i:to_i)) - min(D(from_i:to_i));

     seg_i = [id_i, type_i, from_i, to_i, length_i, center_i, magnitude_i, fluctuation_i,component_1,component_2];
end

function Exist = checkIfExist(from_i,to_i,SegDB,SegDB_names)
  Index_i = SegDB(:,SegDB_names == "from") == from_i & SegDB(:,SegDB_names == "to") == to_i; 
  if sum(Index_i) >= 1
    Exist = true;
  else
    Exist = false;
  end
end

function ParamAprx = command2ParamAprx(Command)
  Arg = command2Arg(Command);
  ParamAprx = Arg2SpikeletParam(Arg);
end

function Arg = command2Arg(Command)
  Arg = struct;
  Split = regexp(Command,'_','split');
  for i = 1:length(Split)
    Arg_i = Split{i};
    if length(Arg_i) >= 10 && strcmp(Arg_i(1:9),'ScutMaTr-')  
      Arg.SuppCutMaTr = str2num(strrep(Arg_i(10:end),'p','.'));   
    elseif length(Arg_i) >= 13 && strcmp(Arg_i(1:12),'LScutMaTrFM-')  
      Arg.LegSuppCutMaTrFM = str2num(strrep(Arg_i(13:end),'p','.'));
    elseif length(Arg_i) >= 7 && strcmp(Arg_i(1:6),'Spike-') 
      Arg.SpikeType= Arg_i(7:end);    
    elseif length(Arg_i) >= 8 && strcmp(Arg_i(1:7),'Select-') 
      Hyphen_i = Arg_i(8:end);   
      Split_i = regexp(Hyphen_i,'-','split');
      Arg.SelectMethod = Split_i{1};
      Arg.SelectParam = str2num(strrep(Split_i{2},'p','.'));
    elseif length(Arg_i) >= strcmp(Arg_i(1:5),'dist-')
       Hyphen_i = Arg_i(6:end);   
       Split_i = regexp(Hyphen_i,'-','split');
       Arg.dist = Split_i{1};
       Arg.subseq_norm = Split_i{2};     
    end
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
