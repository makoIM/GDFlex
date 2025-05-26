function [AnDB_CSE_rev,SegInfo_revised,Info] = Spikelet_SegInfo_adjustPVAndNoiseCut(AnDB_CSE,AnDB_names,SegInfo,training_region,D,Param)
   %% algorithm parameters
   % ZnormDiffThr = 2.25;
   Command = 'dist-Euc-Znorm_ScutMaTr-0p5';
   %% ID constant
   ADJUST = 32;  % adjust id for cutNoisePattern
   
   %% output
   AnDB_CSE_rev = AnDB_CSE; 
   Feature_cutNoise = []; Feature_cutNoise_names = [];
   Feature_shortenPVover5 = []; Feature_shortenPVover5_names = [];
   SegInfo_revised = SegInfo;
   %% DTW distance for CSE
   AnDB_CSE_add = []; SegDB_CSE_add = [];
   for i = 1:size(AnDB_CSE,1)
     AnDbCse_i = AnDB_CSE(i,:); 
     fromCse_i = AnDbCse_i(AnDB_names=="from");
     toCse_i = AnDbCse_i(AnDB_names=="to");        
     fromCseNN_i = AnDbCse_i(AnDB_names=="from_NN");
     toCseNN_i  = AnDbCse_i(AnDB_names=="to_NN");

     % if fromCse_i == 22728 && toCse_i == 22787; fprintf("debug\n");end

     %% cut noise pattern
     % cutResult = cutNoisePattern(D,fromCse_i,toCse_i,fromCseNN_i,toCseNN_i,ZnormDiffThr);
    
    %% extractSegmentForPVover5
    [PvNumOverCut,SegInfo_revised,AnDB_CSE_i,SegDB_new_i,Feature_shortenPVover5_i,Feature_shortenPVover5_names] = Spikelet_AnomalyInfoInstance_extractSegmentForPVover5(AnDbCse_i,AnDB_names,D,Command,SegInfo_revised,training_region,Param);
     if PvNumOverCut % ~isempty(SegDB_new_i)
       AnDB_CSE_rev(i, AnDB_names  == "distLengthNorm") = 0;   
       AnDB_CSE_rev(i, AnDB_names  == "distLengthNorm_pre") = AnDB_CSE(i, AnDB_names  == "distLengthNorm");
       AnDB_CSE_rev(i, AnDB_names  == "adjust") = ADJUST ;
       Feature_shortenPVover5 = [Feature_shortenPVover5; Feature_shortenPVover5_i];
     end    
     AnDB_CSE_rev = [AnDB_CSE_rev;AnDB_CSE_i];
   end

   for i = 1:size(AnDB_CSE_rev,1)
     AnDbCse_i = AnDB_CSE_rev(i,:);   
     [cutResult,Feature_CutNoise_i,Feature_cutNoise_names] = Spikelet_AnomalyInfoInstance_cutNoise_usingDtw(AnDbCse_i,AnDB_names,D,Param);
     % AD6
     % [cutResult,Feature_CutNoise_i,Feature_cutNoise_names] = Spikelet_AnomalyInfoInstance_cutNoise(AnDbCse_i,AnDB_names,D,Param);
     if cutResult
       AnDB_CSE_rev(i, AnDB_names  == "distLengthNorm") = 0;   
       AnDB_CSE_rev(i, AnDB_names  == "distLengthNorm_pre") = AnDB_CSE_rev(i, AnDB_names  == "distLengthNorm");
       AnDB_CSE_rev(i, AnDB_names  == "adjust") = ADJUST ;

       Feature_cutNoise = [Feature_cutNoise; Feature_CutNoise_i];
     end
   end
   %% output
   Info = struct;
   Info.Names = AnDB_names;
   Info.revised = AnDB_CSE_rev;
   Info.input = AnDB_CSE; 
   Info.Feature_cutNoise       = Feature_cutNoise;
   Info.Feature_cutNoise_names = Feature_cutNoise_names;
   Info.Feature_shortenPVover5       = Feature_shortenPVover5;
   Info.Feature_shortenPVover5_names = Feature_shortenPVover5_names;
end

function cutResult = cutNoisePattern(D,from,to,fromNN,toNN,ZnormDiffThr)
     cutResult = false;     

     TimeCand = (from:to);
     TimeNN = (fromNN:toNN);
   
     D_cand   = D(TimeCand);
     D_nn = D(TimeNN);
   
     D_cand_norm = normalize(D_cand);
     D_nn_norm = normalize(D_nn);

     Ddiff = D_cand_norm - D_nn_norm;
     if max(abs(Ddiff)) <= ZnormDiffThr
       cutResult = true;
     end

end


