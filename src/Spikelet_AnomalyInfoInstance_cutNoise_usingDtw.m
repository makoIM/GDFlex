function   [cutResult,DsFeature,DsFeature_names] = Spikelet_AnomalyInfoInstance_cutNoise_usingDtw(AnInst,Names,D,Param)
   %% deafult output
   % cutResult = false;

   %% Feature list
   FeatureList = ["ratio_10pL2_dtwZnorm","diff_magAbsMax_dtwZnorm"];               
   Param_DSF = struct;
   Param_DSF.FeatureList = FeatureList;

   %% Spikelet_SubseqPair_diffSpikeFeature
   if ismember("spikelet_from",Names)
     from_i = AnInst(Names=="spikelet_from");
   elseif ismember("from",Names)
     from_i = AnInst(Names=="from");
   end
   if ismember("spikelet_to",Names)
     to_i = AnInst(Names=="spikelet_to");  
   elseif ismember("to",Names)
     to_i = AnInst(Names=="to");        
   end 
   fromNN_i = AnInst(Names=="from_NN");
   toNN_i  = AnInst(Names=="to_NN");

   % if from_i == 22728 && to_i == 22787; fprintf('debug\n');end

   D_cand = D(from_i:to_i);
   D_nn = D(fromNN_i:toNN_i);
   [DsFeature,DsFeature_names,AdFeature, AdFeature_names] = Spikelet_SubseqPair_diffSpikeFeature(D_cand, D_nn, Param_DSF);
 
   EscEval_i = [];
   for i = 1:length(FeatureList)
         featureName_i = FeatureList(i);
         feature_pos_i = (featureName_i == DsFeature_names);
         feature_i = DsFeature(feature_pos_i);
         EscEval_i = [EscEval_i,feature_i];

         featureName_char_i = char(featureName_i);
         eval([featureName_char_i '= feature_i;']);
   end

   %% discriminate whether a subsequence pair is cut or not
   survive =   ratio_10pL2_dtwZnorm >= 0.95 ...
                | (ratio_10pL2_dtwZnorm >= 0.9 & diff_magAbsMax_dtwZnorm > 2);  
   cutResult = ~survive;
end
