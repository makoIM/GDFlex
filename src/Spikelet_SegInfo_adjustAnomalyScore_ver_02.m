function [AnDB,AnDB_rev,AnDB_rev_names,PairwiseInfo] = Spikelet_SegInfo_adjustAnomalyScore_ver_02(AnDB_org,AnNames,SegInfo_train,Param)
  FuncName = 'Spikelet_SegInfo_adjustAnomalyScore';
  
  %% argument  
  % AdInfo_names = ["id","from","to","type","magnitude","dist","distLengthNorm","from_NN","to_NN","spike_sign"];
  % SegDB_names = ["id","type","from","to","length","center","magnitude","fluctuation"];

  %% output (PairwiseInfo)
  [~,~,~,~,~,~,~,Info_PCSS] ...
    = Spikelet_Pairwise_containTwoSpikesWithSameSign_ver_03([]);
  PairwiseInfo_Stem_Names = ["from","to","type"];
  Info_PCSS_Names = Info_PCSS.names;
  PairwiseInfo_names = [PairwiseInfo_Stem_Names, "distNorm_new", "distNorm_pre",Info_PCSS_Names];
  PairwiseInfo_values = NaN(size(AnDB_org,1),length(PairwiseInfo_names));
  
  %% output
  AnDB = AnDB_org;
  AnDB_rev = []; 
  AnDB_rev_names = ["id","from","to"...
                         ,"distLengthNorm_org","distLengthNorm_new","distLengthNorm_ratio"...
                         ,"boundary_NN","fromAT_NN","toAT_NN"...
                         ,"fromAT_test","toAT_test"];

  %% argument
  MagInfo = SegInfo_train.MagInfo;
  D = MagInfo.data_org;
  SegDB = SegInfo_train.SegDB;
  SegNames = SegInfo_train.SegDB_names;
  distLengthNorm_top1 = AnDB_org(1,AnNames=="distLengthNorm");
  % last_id = SegInfo.last_id;

  %% main loop for AnDB to SegDB_2
  Index_Seg2 = abs(SegDB(:,SegNames=="type")) == 2;
  SegDB_2 = SegDB(Index_Seg2,:);
  info_PCSS_i = [];
  for i = 1:size(AnDB,1)
    AnDB_i = AnDB(i,:);
    fromTest_i = AnDB_i(AnNames=="from");
    toTest_i = AnDB_i(AnNames=="to");  
    type_i = AnDB_i(AnNames=="type");

    PairwiseInfo_values(i,1:3) = [fromTest_i,toTest_i,type_i];
    result_i = false;
    if abs(type_i) == 2
 
      DEBUG = true; % DEBUG = false; 
      if DEBUG && isequal([fromTest_i,toTest_i],[2864, 2977])
        fprintf('[%s] debug\n',FuncName);
      end
      fromNN_i = AnDB_i(AnNames=="from_NN");
      toNN_i = AnDB_i(AnNames=="to_NN");
      
      % = Spikelet_Pairwise_containTwoSpikesWithSameSign(D,type_i,fromTest_i,toTest_i,fromNN_i,toNN_i,SegDB_2,SegNames,Param);
      ContainTwoSpikesWithSameSign_VERSION = 3;
      if ContainTwoSpikesWithSameSign_VERSION == 3
        [result_i,distLengthNorm_i,boundary_NN_i,fromAT_NN_i,toAT_NN_i,fromAT_test_i,toAT_test_i,Info_PCSS_i] ...
        = Spikelet_Pairwise_containTwoSpikesWithSameSign_ver_03(D,type_i,fromTest_i,toTest_i,fromNN_i,toNN_i,SegDB_2,SegNames,distLengthNorm_top1,Param);   
      elseif ContainTwoSpikesWithSameSign_VERSION == 2
        [result_i,distLengthNorm_i,boundary_NN_i,fromAT_NN_i,toAT_NN_i,fromAT_test_i,toAT_test_i,Info_PCSS_i] ...
        = Spikelet_Pairwise_containTwoSpikesWithSameSign_ver_02(D,type_i,fromTest_i,toTest_i,fromNN_i,toNN_i,SegDB_2,SegNames,Param);
      end

      if result_i && distLengthNorm_i > AnDB_org(i,AnNames=="distLengthNorm")  

        info_PCSS_i = Info_PCSS_i.values;    
        %% revise AnDB
        AnDB(i,AnNames=="distLengthNorm_pre") = AnDB_org(i,AnNames=="distLengthNorm");
        AnDB(i,AnNames=="adjust") = 1;
        AnDB(i,AnNames=="distLengthNorm") = distLengthNorm_i;
        distLengthNorm_org_i = AnDB_org(i,AnNames=="distLengthNorm");
        AnDB_rev = [AnDB_rev...
                   ;AnDB_i(AnNames=="id"),AnDB_i(AnNames=="from"),AnDB_i(AnNames=="to")...
                   ,distLengthNorm_org_i,distLengthNorm_i,distLengthNorm_i/distLengthNorm_org_i...
                   ,boundary_NN_i,fromAT_NN_i,toAT_NN_i...
                   ,fromAT_test_i,toAT_test_i];

        %% revise PairwiseInfo
        PairwiseInfo_values(i,length(PairwiseInfo_Stem_Names)+1:end) = [distLengthNorm_i,distLengthNorm_org_i,info_PCSS_i(:,1)'];
        PairwiseInfo_Thres = info_PCSS_i(:,2)';
      end   
    end
  end

  % AnDB
  [sort_val, sort_pos] = sort(AnDB(:,AnNames=="distLengthNorm"),'descend');
  AnDB = AnDB(sort_pos,:);

  % PairwiseInfo
  PairwiseInfo_values = PairwiseInfo_values(sort_pos,:);
  PairwiseInfo = struct;
  PairwiseInfo.names = PairwiseInfo_names; 
  PairwiseInfo.values = PairwiseInfo_values;
  if exist('PairwiseInfo_Thres','var')
    PairwiseInfo.threshold = PairwiseInfo_Thres;
  end
end


