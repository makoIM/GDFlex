function Info = Spikelet_SpikeInSpike_constantSegmentation(MagInfo,Param)

  %% parameters
  ParamIntra = Param.intra;
  % SegStatistics = Param.SegStatistics;

  MagBandWidthRate = ParamIntra.Const_MagBandRate; 
  LengthRate = ParamIntra.Const_LengthRate; 

  %% extract spike using constant segmentation 
  %[ConstRange,SegmentList,Names] = Spikelet_MagInfo_PvAbs2_constantExtraction(MagInfo,Const_MagBandRate,Const_LengthRate);

  % MagBandRate = 0.25 LengthRate = 0.5
  % Rate_PeakInCenter = 0.5;
  if ~exist('Method','var')
    MethodSearch = "mean_first_cross";
    MethodCenter = "median"; % MethodCenter = "mean";
  end

  SegmentList_Names = ["from","to","type","center","parentSpikeId"];  %% magnitude can be known for z-normalized subsequences; 
  %% subsequence extraction
  D = MagInfo.data_org;
  Mag = MagInfo.magnitude;
  [MagAbsMax, peak] = max(abs(Mag)); 
  magVal = Mag(peak);
  magSign = sign(magVal);
  %% constant MagBand extraction
  MagBand = MagAbsMax*MagBandWidthRate;
  [Clen, Names] = Spikelet_Subseq_ConstantExtractionbyMagBand(D,MagBand);

  if isempty(Clen) 
    ConstRange = [];
    LocalSegmentList = [];
  else
      from_org = Clen(1,Names=="from");
      to_org   = Clen(1,Names=="to");
      ConstRange = [from_org,to_org];
      if Clen(1,Names=="length") < length(D)*LengthRate
        LocalSegmentList = []; 
      else    
        if MethodCenter == "median"
          const_center_val = median(D(from_org:to_org));
        elseif MethodCenter == "mean"
          const_center_val = mean(D(from_org:to_org));  
        end

        if to_org < peak 
          FindOption = "monotone_extension";  
          from_const_ext = find_constant_left(peak,D,magSign,const_center_val,FindOption);
          FindOption = "default";  
          from_const_def = find_constant_left(peak,D,magSign,const_center_val,FindOption);
          type_const = magSign;
          % from_const = find_constant_right(1,D,(-1)*magSign,const_center_val);
          %                  ["from",        "to",    "type","center"];
          LocalSegmentList = [from_const_ext,length(D),(-1)*magSign,peak,peak...
                             ;from_const_def,length(D),(-1)*magSign,peak,peak...
                             ;1,peak,magSign,peak,peak];
        elseif peak < from_org
          FindOption = "monotone_extension";    
          to_const_ext = find_constant_right(peak,D,magSign,const_center_val,FindOption);
          FindOption = "default";    
          to_const_def = find_constant_right(peak,D,magSign,const_center_val,FindOption);
          % to_const = find_constant_left(length(D),D,(-1)*magSign,const_center_val);
          LocalSegmentList = [1,to_const_ext,magSign,peak,peak...
                             ;1,to_const_def,magSign,peak,peak...
                             ;peak,length(D),(-1)*magSign,peak,peak];
        else
          % [from_const,to_const] = search_lr(peak, D, MagBand*Rate_PeakInCenter); 
          LocalSegmentList = [1,   peak,     magSign,peak,peak...      % left leg
                             ;peak,length(D),(-1)*magSign,peak,peak];  % right leg   
        end
      end
  end

  %%
  Info = struct;
  Info.constantRange = ConstRange; 
  Info.segmentList = LocalSegmentList;
  Info.segmentList_names = SegmentList_Names;
end

function to_constant = find_constant_left(mag_pos,D,mag_sign,constant_center,Option)
   if ~exist('Option','var')
     Option = "default";
   end

   F2P = (1:mag_pos);
   if mag_sign > 0
       to_cross = find( D(F2P) <= constant_center, 1,'last');
   else
       to_cross = find( D(F2P) >= constant_center, 1,'last');
   end

   if Option == "default"
     to_cross_pre = to_cross + 1;
     [min_val, min_pos] = min([abs(D(to_cross)-constant_center),abs(D(to_cross_pre)-constant_center)]);
     if min_pos == 1
       to_constant = to_cross;  
     elseif min_pos == 2
       to_constant = to_cross_pre;
     end
   elseif Option == "monotone_extension"
     to_constant = to_cross;  
     if mag_sign > 0   
       for i = to_cross-1:-1:1
         if (D(i) - D(i+1)) < 0
           to_constant = i;
         else
           break;  
         end
       end
     elseif mag_sign < 0   
       for i = to_cross-1:-1:1
         if (D(i) - D(i+1)) > 0
           to_constant = i;
         else
           break;  
         end
       end
     end
   end
end

function to_constant =  find_constant_right(mag_pos,D,magSign,constant_center,Option)
   if ~exist('Option','var')
     Option = "default";
   end

   P2T = (mag_pos:length(D));
   if magSign > 0
       to_cross_rel = find( D(P2T) <= constant_center, 1,'first');
   else
       to_cross_rel = find( D(P2T) >= constant_center, 1,'first');
   end
   to_cross = mag_pos + to_cross_rel -1;

   if Option == "default"
     to_cross_pre = to_cross -1;
     [min_val, min_pos] = min([abs(D(to_cross)-constant_center),abs(D(to_cross_pre)-constant_center)]);
     if min_pos == 1
       to_constant = to_cross;  
     elseif min_pos == 2
       to_constant = to_cross_pre;
     end

   elseif Option == "monotone_extension"
     to_constant = to_cross;  
     if magSign > 0   
       for i = to_cross+1:length(D)-1
         if (D(i+1) - D(i)) < 0
           to_constant = i;
         else
           break;  
         end
       end
     elseif magSign < 0   
       for i = to_cross+1:length(D)-1
         if (D(i+1) - D(i)) > 0
           to_constant = i;
         else
           break;  
         end
       end
     end
   end
end

function [from,to] = search_lr(i, D, MagThr)
  count_right = search_right(i, D, MagThr);    
  count_left = search_left(i, D, MagThr);    
  from = i - (count_left -1);
  to = i + count_right -1;
end

function count = search_right(pos_i,D,MagThr)
  MagThr_half = MagThr/2;
  count = 0;
  for i = pos_i:length(D)
      if abs(D(i) - D(pos_i)) <= MagThr_half
        count = count + 1;
      else
        break;  
      end
  end
end

function count = search_left(pos_i,D,MagThr)
  MagThr_half = MagThr/2;
  count = 0;
  for i = pos_i:-1:1
      if abs(D(i) - D(pos_i)) <= MagThr_half
        count = count + 1;
      else
        break;  
      end
  end
end
