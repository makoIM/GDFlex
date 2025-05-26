function [Feature,Feature_names, AdFeature, AdFeature_names] = Spikelet_SubseqPair_diffSpikeFeature(D1_Raw, D2_Raw, Param)
     if exist('Param','var') && ~isempty(Param) && isfield(Param,"Command")
       Command = Param.Command;
     else
       Command = 'dist-Euc-Znorm_ScutMaTr-0p5';
     end
     if exist('Param','var') && ~isempty(Param) && isfield(Param,"FeatureList")
       FeatureList = Param.FeatureList;  
     else
       FeatureList = ["ratio_EucDtwL2_raw","ratio_EucDtwL2_znorm"...
                    ,"ratio_insideL2_zeroCross_eucRaw", "ratio_insideL2_spike_eucRaw"...
                    ,"ratio_insideL2_zeroCross_eucZnorm","ratio_insideL2_spike_eucZnorm"...
                    ,"ratio_insideL2_zeroCross_dtwRaw","ratio_insideL2_spike_dtwRaw"...
                    ,"ratio_insideL2_zeroCross_dtwZnorm","ratio_insideL2_spike_dtwZnorm"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_eucRaw"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_eucZnorm"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_dtwRaw"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_dtwZnorm"...
                    ,"ratio_10pL2_dtwZnorm","ratio_leftL2_dtwZnorm","ratio_rightL2_dtwZnorm","diff_magAbsMax_dtwZnorm"];
       FeatureList = ["ratio_EucDtwL2_raw","ratio_EucDtwL2_znorm"...
                    ,"ratio_insideL2_spike_eucRaw"...
                    ,"ratio_insideL2_spike_eucZnorm"...
                    ,"ratio_insideL2_spike_dtwRaw"...
                    ,"ratio_insideL2_spike_dtwZnorm"...
                    ,"ratio_insideL2Top2CS_spike_eucRaw"...
                    ,"ratio_insideL2Top2CS_spike_eucZnorm"...
                    ,"ratio_insideL2Top2CS_spike_dtwRaw"...
                    ,"ratio_insideL2Top2CS_spike_dtwZnorm"...
                    ,"ratio_insideL2_zeroCross_eucRaw"...
                    ,"ratio_insideL2_zeroCross_eucZnorm"...
                    ,"ratio_insideL2_zeroCross_dtwRaw"...
                    ,"ratio_insideL2_zeroCross_dtwZnorm"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_eucRaw"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_eucZnorm"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_dtwRaw"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_dtwZnorm"...
                    ,"ratio_10pL2_dtwZnorm","diff_magAbsMax_dtwZnorm"...
                    ,"ratio_leftL2_dtwZnorm","ratio_rightL2_dtwZnorm"];    
        FeatureList = ["ratio_EucDtwL2_raw","ratio_EucDtwL2_znorm"...
                    ,"ratio_insideL2_spike_eucRaw"...
                    ,"ratio_insideL2_spike_eucZnorm"...
                    ,"ratio_insideL2_spike_dtwRaw"...
                    ,"ratio_insideL2_spike_dtwZnorm"...
                    ,"ratio_insideL2Top2CS_spike_eucRaw"...
                    ,"ratio_insideL2Top2CS_spike_eucZnorm"...
                    ,"ratio_insideL2Top2CS_spike_dtwRaw"...
                    ,"ratio_insideL2Top2CS_spike_dtwZnorm"...
                    ,"ratio_insideL2_zeroCross_eucRaw"...
                    ,"ratio_insideL2_zeroCross_eucZnorm"...
                    ,"ratio_insideL2_zeroCross_dtwRaw"...
                    ,"ratio_insideL2_zeroCross_dtwZnorm"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_eucRaw"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_eucZnorm"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_dtwRaw"...
                    ,"ratio_insideL2Top2DiffSign_zeroCross_dtwZnorm"...
                    ,"ratio_10pL2_dtwZnorm"...
                    ,"ratio_leftL2_dtwZnorm","ratio_rightL2_dtwZnorm"...
                    ,"diff_magAbsMax_eucZnorm","diff_magAbsMax_dtwZnorm"]; 
     end

     %{
      "ratio_EucDtw_L2"  => "ratio_EucDtwL2_raw"
                            "ratio_EucDtwL2_znorm"
       ratio_10p_distL2_dtw => ratio_10pL2_dtwZnorm
       ratioL2_left_dtw  =>   ratio_leftL2_dtwZnorm
       ratioL2_right_dtw =>    ratio_leftL2_dtwZnorm
       diffDtw_magMax    =>   diffMagMax_dtwZnorm
     %}

   %% initialize output
   Feature = [];
   Feature_names = [];

   %% znormalization / dtw alignment path

   if contain_inAnyOf(FeatureList,"eucZnorm") 
     D1_Znorm = normalize(D1_Raw);
     D2_Znorm = normalize(D2_Raw);
   end

   if contain_inAnyOf(FeatureList,"dtwRaw") 
     [distDtw,PathD1, PathD2] = dtw(D1_Raw, D2_Raw);
     D1_RawDtw = D1_Raw(PathD1);
     D2_RawDtw = D2_Raw(PathD2);
   end

   if contain_inAnyOf(FeatureList,"dtwZnorm") 
     D1_Znorm = normalize(D1_Raw);
     D2_Znorm = normalize(D2_Raw);
     [distDtw,PathD1, PathD2] = dtw(D1_Znorm, D2_Znorm);
     D1_ZnormDtw = D1_Znorm(PathD1);
     D2_ZnormDtw = D2_Znorm(PathD2);
   end

   %% calculate feature
   [TF,Pos_org] = ismember(FeatureList,["ratio_EucDtwL2_raw"]);
   if sum(TF) ~= 0
     ratio_EucDtwL2 = ratio_EucDTWL2_calcurate(D1_Raw,D2_Raw,D1_RawDtw,D2_RawDtw);  
     
     % revise (Feature,Feature_names,FeatureList)
     Feature = [Feature,ratio_EucDtwL2];
     Feature_names = [Feature_names,"ratio_EucDtwL2_raw"];   
     FeatureList = FeatureList(~TF);
   end

   [TF,Pos_org] =ismember(FeatureList,["ratio_EucDtwL2_znorm"]);
   if sum(TF) ~= 0
     ratio_EucDtwL2 = ratio_EucDTWL2_calcurate(D1_Znorm,D2_Znorm,D1_ZnormDtw,D2_ZnormDtw);
     
     % revise (Feature,Feature_names,FeatureList)
     Feature = [Feature,ratio_EucDtwL2];
     Feature_names = [Feature_names,"ratio_EucDtwL2_znorm"];   
     FeatureList = FeatureList(~TF);
   end

   %% EUC
   TargetFeatureNames_stem = ["ratio_insideL2_spike","ratio_insideL2Top2CS_spike","ratio_insideL2_zeroCross","ratio_10pL2","ratio_leftL2","ratio_rightL2","diff_magAbsMax","ratio_insideL2Top2DiffSign_zeroCross"];
   AdFeature = []; AdFeature_names = [];
   % ratio_insideL2 eucRaw
   if contain_inAnyOf(FeatureList,"eucRaw") 
     TargetFeatureNames = strcat(TargetFeatureNames_stem,"_eucRaw");
     [Feature,Feature_names,FeatureList,AdFeature_eucRaw,AdFeature_names_stem] = call_diffSpikeFeature(D1_Raw,D2_Raw,TargetFeatureNames,FeatureList,Feature,Feature_names,Command); 
     AdFeature = [AdFeature,AdFeature_eucRaw];  AdFeature_names = [AdFeature_names, AdFeature_names_stem + "_eucRaw"];
   end
   % ratio_insideL2 eucZnorm
   if contain_inAnyOf(FeatureList,"eucZnorm") 
     TargetFeatureNames = strcat(TargetFeatureNames_stem,"_eucZnorm");
     [Feature,Feature_names,FeatureList,AdFeature_eucZnorm,AdFeature_names_stem] = call_diffSpikeFeature(D1_Znorm,D2_Znorm,TargetFeatureNames,FeatureList,Feature,Feature_names,Command);
     AdFeature = [AdFeature,AdFeature_eucZnorm];  AdFeature_names = [AdFeature_names, AdFeature_names_stem + "_eucZnorm"];
   end
   %% DTW
   % ratio_insideL2 dtwRaw
   if contain_inAnyOf(FeatureList,"dtwRaw") 
     TargetFeatureNames = strcat(TargetFeatureNames_stem,"_dtwRaw");
     [Feature,Feature_names,FeatureList,AdFeature_dtwRaw,AdFeature_names_stem] = call_diffSpikeFeature(D1_RawDtw,D2_RawDtw,TargetFeatureNames,FeatureList,Feature,Feature_names,Command);  
     AdFeature = [AdFeature,AdFeature_dtwRaw];  AdFeature_names = [AdFeature_names, AdFeature_names_stem + "_dtwRaw"];
   end
   % ratio_insideL2 dtwZnorm
   if contain_inAnyOf(FeatureList,"dtwZnorm") 
     TargetFeatureNames = strcat(TargetFeatureNames_stem,"_dtwZnorm");  
     [Feature,Feature_names,FeatureList,AdFeature_dtwZnorm,AdFeature_names_stem] = call_diffSpikeFeature(D1_ZnormDtw,D2_ZnormDtw,TargetFeatureNames,FeatureList,Feature,Feature_names,Command); 
     AdFeature = [AdFeature,AdFeature_dtwZnorm];  AdFeature_names = [AdFeature_names, AdFeature_names_stem + "_dtwZnorm"];
   end
end

function [Feature,Feature_names, FeatureList,adfeature,adfeatureName] = call_diffSpikeFeature(D1,D2,TargetFeatureNames,FeatureList,Feature,Feature_names,Command); 
   [TF,Pos_org] = ismember(FeatureList,TargetFeatureNames);
   if sum(TF) ~= 0
     [feature, featureName_calculated,adfeature,adfeatureName] = diffSpikeFeature(D1,D2,Command);
     FeatureList_in = FeatureList(TF);
     for i = 1:length(FeatureList_in)
       featureName_i = FeatureList_in(i);  
       [Bool,fnc_pos_i] = contained_inAnyOf(featureName_calculated,featureName_i);
       Feature = [Feature,feature(fnc_pos_i)];
       Feature_names = [Feature_names,featureName_i];    
     end
     FeatureList = FeatureList(~TF);
   end
end

function ratio_EucDtwL2 = ratio_EucDTWL2_calcurate(D1,D2,D1_dtw,D2_dtw)    
     LP= 2;
     Ddiff_raw = D1 - D2;  
     Ddiff_dtw = D1_dtw - D2_dtw;
     dist_lenNorm_euc = norm(Ddiff_raw,LP)/sqrt(length(Ddiff_raw));
     dist_lenNorm_dtw = norm(Ddiff_dtw,LP)/sqrt(length(Ddiff_dtw));
     ratio_EucDtwL2 =  dist_lenNorm_dtw/dist_lenNorm_euc;
end  


function [Feature, Feature_names, AdFeature,AdFeature_names] = diffSpikeFeature(D1,D2,Command)
   if ~exist("LP","var")
     LP = 2;
   end
   %% 
   Ddiff = D1 - D2;
   len_whole = length(Ddiff);
   distL2_whole = norm(Ddiff,2);
   %% magnitude statistics
   ParamAprx = command2ParamAprx(Command);
   MagInfo = Spikelet_aproximation_ver_05(Ddiff,ParamAprx);
   Mag = MagInfo.magnitude;
   MagAbs_ZC = abs(Mag);
   [MagAbsVal_sorted,MagAbsPos_sorted] = sort(MagAbs_ZC,'descend');
   
   %% zero cross statistics
   % MagInfo_ZC = Spikelet_MagInfo_zeroCross_ver_1(MagInfo);
   MagInfoZC = Spikelet_MagInfoZC_calculate(Ddiff);
   Center = MagInfoZC.center;

   % diff_magAbsMax
   Mag_ZC = MagInfoZC.magnitude;
   MagAbs_ZC = abs(Mag_ZC);
   diff_magAbsMax = max(MagAbs_ZC);
   
   Debug_MagInfo = true; Debug_MagInfo = false;
   if Debug_MagInfo
     figure(); 
     subplot(3,1,1); plot(D1,'b','LineWidth',2); hold on; plot(D2,'r');
     subplot(3,1,2); plot(Ddiff); 
     subplot(3,1,3); plot(MagInfo.magnitude);
   end
   %% Zerocross
   topK = 2;
   [distL2_topK,Len_topK,Mag_topK,Center_topK,From_topK,To_topK] = MagInfoZC_dist_topK(MagInfoZC,topK,LP);
   % [distL2_topK,Len_topK,Mag_topK,From_topK,To_topK] = dist_topK_byMag(MagInfo_ZC,Ddiff,MagAbsPos_sorted,topK,LP);
   ratio_insideL2_zeroCross = 0;  ratio_insideL2Top2DiffSign_zeroCross = 0;
   zeroCrossTop1_from = NaN;  zeroCrossTop1_to = NaN;
   zeroCrossTop2_from = NaN; zeroCrossTop2_to = NaN; 
   if length(Center) >= 1
     ratio_insideL2_zeroCross = distL2_topK(1)/distL2_whole;   
     zeroCrossTop1_from = From_topK(1);  zeroCrossTop1_to = To_topK(1);
      % ratio_insideL2_zeroCross = (distL2_top1*sqrt(len_top1))/(distL2_whole*sqrt(len_whole)); 
   end
   if length(Center) >= 2
     distL2_top12 = 0;  
     zeroCrossTop2_from = From_topK(2);  zeroCrossTop2_to = To_topK(2); 
     if sign(Mag_topK(1)*Mag_topK(2)) == -1
       Range_1 = [From_topK(1),To_topK(1)];  Range_2 = [From_topK(2),To_topK(2)];     
       [OVERLAP, overlap_length,ov_type] = Spikelet_overlapSegment(Range_1,Range_2);
       if OVERLAP
          from_12 = min([From_topK(1),From_topK(2)]); to_12 = max([To_topK(1),To_topK(2)]);
          distL2_top12 = norm(Ddiff(from_12:to_12),LP);
         
       else
          % distL2_1 = norm(Ddiff(Range_1(1):Range_1(2)),LP);distL2_2 = norm(Ddiff(Range_2(1):Range_2(2)),LP);
          distL2_1 = distL2_topK(1);   distL2_2 = distL2_topK(2);  
          distL2_top12_sqare = distL2_1^2 + distL2_2^2;
          distL2_top12 = sqrt(distL2_top12_sqare);
          
          %{ 
          Range_between = [min(To_topK):max(From_topK)];
          distL2_between = norm(Ddiff(Range_between),LP);
          distL2_top12_sqare = distL2_topK(1)^2 + distL2_topK(2)^2 - distL2_between^2;
          if distL2_top12_sqare > 0
            distL2_top12 = sqrt(distL2_top12_sqare);
          end
          %}
       end    
     end
     ratio_insideL2Top2DiffSign_zeroCross = distL2_top12/distL2_whole;
   end

   %% Spike
   LP = 2; topK = 2;
   [distL2_topK,len_top1,Mag_topK,Center_topK,From_topK,To_topK] = MagInfo_dist_topK(MagInfo,Ddiff,MagAbsPos_sorted,topK,LP);
   ratio_insideL2_spike = 0;  ratio_insideL2Top2CS_spike = 0;  
   spikeTop1_from = NaN;  spikeTop1_to = NaN;
   spikeTop12_from = NaN; spikeTop12_to = NaN; 
   if length(Center_topK) >= 1
     ratio_insideL2_spike = distL2_topK(1)/distL2_whole; 
     spikeTop1_from = From_topK(1); spikeTop1_to = To_topK(1); 
   end
   if length(Center_topK) >= 2
     spikeTop2_from = From_topK(2); spikeTop2_to = To_topK(2);   
     from_1 = From_topK(1); to_1 = To_topK(1); Range_1 = [from_1, to_1];
     from_2 = From_topK(2); to_2 = To_topK(2); Range_2 = [from_2, to_2];
     [OVERLAP, overlap_length,ov_type] = Spikelet_overlapSegment(Range_1,Range_2);
     if OVERLAP
       from_12 = min([From_topK(1),From_topK(2)]); to_12 = max([To_topK(1),To_topK(2)]);
       distL2_top12 = norm(Ddiff(from_12:to_12),LP); 
       ratio_insideL2Top2CS_spike = distL2_top12/distL2_whole;

       spikeTop12_from = from_12; spikeTop12_to = to_12; 
     end
   end


   %% ratio_leftL2
   RangeOutside_left = Spikelet_findNearestZeroCross(1,Ddiff(1),Ddiff);
   from_left = 1; to_left = RangeOutside_left(2);
   dist_leftL2 = norm(Ddiff(from_left:to_left),LP);
   ratio_leftL2 = dist_leftL2/distL2_whole; 
   
   %% ratio_rightL2
   RangeOutside_right = Spikelet_findNearestZeroCross(length(Ddiff),Ddiff(end),Ddiff);
   from_right = RangeOutside_right(1); to_right = length(Ddiff);
   dist_rightL2 = norm(Ddiff(from_right:to_right),LP);
   ratio_rightL2 = dist_rightL2/distL2_whole; 

   %% ratio_10pL2
   [sort_val, sort_order] = sort(abs(Ddiff),'descend');
   sort_order_10p = sort_order(1:round(length(sort_order)*0.1));
   D_diff_10p = Ddiff(sort_order_10p);
   dist_10pL2 = norm(D_diff_10p,LP);
   ratio_10pL2 = dist_10pL2/distL2_whole;

   %% Feature
   Feature_names = ["ratio_insideL2_zeroCross","ratio_insideL2Top2DiffSign_zeroCross" ...
                    ,"ratio_insideL2_spike",   "ratio_insideL2Top2CS_spike"...
                    ,"ratio_10pL2",            "diff_magAbsMax"...
                    ,"ratio_leftL2",           "ratio_rightL2"];
   Feature = [ratio_insideL2_zeroCross, ratio_insideL2Top2DiffSign_zeroCross...
             ,ratio_insideL2_spike,     ratio_insideL2Top2CS_spike...
             ,ratio_10pL2,              diff_magAbsMax...
             ,ratio_leftL2,             ratio_rightL2];

   AdFeature_names = ["spikeTop1_from","spikeTop1_to","spikeTop2_from","spikeTop2_to"...
                     ,"zeroCrossTop1_from","zeroCrossTop1_to","zeroCrossTop2_from","zeroCrossTop2_to"...
                     ,"to_left","from_right"];
   AdFeature = [spikeTop1_from,spikeTop1_to,spikeTop12_from,spikeTop12_to...
               ,zeroCrossTop1_from,zeroCrossTop1_to,zeroCrossTop2_from,zeroCrossTop2_to...
               ,to_left,from_right];
end

function [Bool,pos] = contained_inAnyOf(StrList,Str)
   Bool = false;

   for i = 1:length(StrList)
     Bool = contains(Str,StrList(i));
     if Bool == true
       pos = i;  
       return;
     end
   end
end

function [Bool,pos] = contain_inAnyOf(StrList,Str)
   Bool = false;

   for i = 1:length(StrList)
     Bool = contains(StrList(i),Str);
     if Bool == true
       pos = i;  
       return;
     end
   end

end

function  [Dist_sorted,Len,Mag,Center_sorted,From,To] = MagInfoZC_dist_topK(MagInfo_ZC,topK,LP)  
  Ddiff = MagInfo_ZC.data_org;
  Center = MagInfo_ZC.center;
  if isempty(Center)
    Dist_sorted = []; Len =[]; Mag = [];Center_sorted = []; From = []; To = [];
  end

  Dist = NaN(length(Center),1);
  for i = 1:length(Center)
    center_i = Center(i);  
    from = MagInfo_ZC.left(i);
    to =  MagInfo_ZC.right(i);
    Dist(i) = norm(Ddiff(from:to),LP);
  end

  if length(Center) >= 1
    [Dist_sorted,Order] = sort(Dist,"descend");
    Center_sorted = Center(Order);
    Mag = MagInfo_ZC.magnitude(Order);
    From = MagInfo_ZC.left(Order);
    To = MagInfo_ZC.right(Order);
    Len = To - From + 1;
  end
end


function  [distLP_topK,Len,Mag,Center,From,To] = MagInfo_dist_topK(MagInfo,Ddiff,MagAbsPos_sorted,topK,LP)  
  magTopPos = MagAbsPos_sorted(1:topK);
  Center = MagAbsPos_sorted;
  Mag = MagInfo.magnitude(magTopPos);
  From = MagInfo.left(magTopPos);
  To =  MagInfo.right(magTopPos);
  Len = To - From + 1;

  distLP_topK = NaN(length(From),topK);
  for i = 1:length(From)
    from_i = From(i);
    to_i = To(i);
    distLP_topK(i) = norm(Ddiff(from_i:to_i),LP);
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
