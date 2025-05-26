function [AdInfo,AdInfo_names] = Spikelet_AD_TrTe_STOMP_ver_02(SegInfo,training_region,Param)                                 
  FuncName = 'Spikelet_AD_TrTe_STOM_ver_02';
  AdInfo_names = ["id","from","to","type","magnitude","dist","distLengthNorm","from_NN","to_NN","spike_sign","adjust","distLengthNorm_pre"];
  
  %% argument
  DIST = Param.dist;
  % FeatureType = "C22"; FeatureType = "Euclidean";
  
  if DIST == "Euc"
    EXCLUSION_WINDOW_RATE = 0.2;

    EXCLUSION_FlacToFlac_0_RATE = [0.4, 2];     % 2024-12-03 AD4
    EXCLUSION_FlacToFlac_0_RATE = [0.35, 2.5];  % 2024-12-07 AD6
    EXCLUSION_FlacToFlac_1_RATE = [0.4, 2];  % [0.5, 2] [0,inf]; 
    % EXCLUSION_FlacToFlac_1_RATE = [0.3,2.5]; % 2024-11-02
    
    EXCLUSION_FlacToFlac_2_RATE = [0.75,1.5]; % [0.75,2];
    % EXCLUSION_FlacToFlac_2_RATE = [0.8,1.5];  % 2024-12-02
      
    EXCLUSION_FlacToFlac_3_RATE = [0.75,2];  % 2024-11-07 CS0
    EXCLUSION_FlacToFlac_3_RATE = [0.5,2];  % 2024-11-30 DA2
    EXCLUSION_FlacToFlac_3_RATE = [0.75,2];  % 2024-11-07 DA3
    EXCLUSION_FlacToFlac_3_RATE = [0.75,2.5];  % 2024-12-02 AD3

    EXCLUSION_FlacToFlac_5_RATE = [0.5,2];  
    % EXCLUSION_FlacToFlac_5_RATE = [0.6,2]; 

    % EXCLUSION_FlacToFlac_3_RATE = [0.4,2]; % 2024-11-10 226 test fail

    EXCLUSION_MatToMat_RATE = [0.5,2];    % [0.5,2]; [0.6,inf]; 
    % EXCLUSION_MatToMat_RATE = [0.6,2];  % 2024-12-02
    % EXCLUSION_MatToMat_RATE = [0,2];  % 2024-11-03
  elseif DIST == "C22"
    EXCLUSION_WINDOW_RATE = 0.2;
    EXCLUSION_FlacToFlac_1_RATE = [0.4, 2]; % [0.5, 2] [0,inf]; 
    EXCLUSION_FlacToFlac_2_RATE = [0.75,1.5]; % [0.75,2];
    EXCLUSION_FlacToFlac_3_RATE = [0.75,2];
    EXCLUSION_MatToMat_RATE     = [0.5,2]; % [0.5,2]; [0.6,inf]; 
  end

  D_all = SegInfo.data_org; 
  D_train = SegInfo.data_org(1:training_region);
  SegDB = SegInfo.SegDB;
  SegDB_names = SegInfo.SegDB_names;
  MagInfo_ini = SegInfo.MagInfo.output.generateInitialSpikelet;

  AdInfo = NaN(size(SegDB,1),length(AdInfo_names));
  for i = 1:size(SegDB,1)  % i:test, j:train
    if i == 1; fprintf('dist-calc(%d %d) ',size(SegDB,1),length(D_train)); UNIT = 50; end;   
    %% test subsequence 
    spikeId_i = SegDB(i,SegDB_names=="id");
    from_i = SegDB(i,SegDB_names=="from");
    to_i = SegDB(i,SegDB_names=="to");
    mag_i = SegDB(i,SegDB_names=="magnitude");
    fluctuation_i = abs(SegDB(i,SegDB_names=="fluctuation"));
    type_i =  SegDB(i,SegDB_names=="type");
    sign_i = sign(type_i);
    subseq_len_i = to_i - from_i + 1;
    subseq_i = D_all(from_i:to_i);

    DEBUG_i = true; DEBUG_i = false; 
    DEBUG_Range_i = false;
    if DEBUG_i && isequal([from_i, to_i],[4019, 4097]) 
      fprintf('[%s] debug spikeId %d\n',FuncName, spikeId_i);
      figure(); T_i = (from_i:to_i); plot(T_i,D_all(T_i));
      DEBUG_Range_i = true;
    end

    EXCLUSION_ZONE_i = false(length(D_train)-subseq_len_i+1,1);
    if length(D_train) == length(D_all)
      [OVERLAP_i, overlap_length_i] = Spikelet_overlapSegment([1,length(D_train)],[from_i,to_i]);
      if overlap_length_i/subseq_len_i <= EXCLUSION_WINDOW_RATE
        EXCUSION_LENGHT_i = ceil(subseq_len_i*EXCLUSION_WINDOW_RATE);  
        Ex_from_i = max(1,from_i-EXCUSION_LENGHT_i); Ex_to_i = min(length(D_train),to_i+EXCUSION_LENGHT_i); 
        EXCLUSION_ZONE_i(Ex_from_i :Ex_to_i) = true(Ex_to_i-Ex_from_i+1,1); 
      end
    end
   
    if type_i == 0
      % EXCLUSION_ZONE_Mag_i = false(length(D_train)-subseq_len_i+1,1);

      EXCLUSION_FlacToFlac_RANGE = [fluctuation_i*EXCLUSION_FlacToFlac_0_RATE(1),fluctuation_i*EXCLUSION_FlacToFlac_0_RATE(2)];         
      EXCLUSION_ZONE_FlacToFlac_i = exclusion_byFlacToFlac(D_train,subseq_len_i,EXCLUSION_FlacToFlac_RANGE);
      EXCLUSION_ZONE_Mag_i = EXCLUSION_ZONE_FlacToFlac_i;
    end

    if abs(type_i) == 1  
      EXCLUSION_FlacToFlac_RANGE = [fluctuation_i*EXCLUSION_FlacToFlac_1_RATE(1),fluctuation_i*EXCLUSION_FlacToFlac_1_RATE(2)];         
      EXCLUSION_ZONE_FlacToFlac_i = exclusion_byFlacToFlac(D_train,subseq_len_i,EXCLUSION_FlacToFlac_RANGE);
      EXCLUSION_ZONE_Mag_i = EXCLUSION_ZONE_FlacToFlac_i;
    end
    if abs(type_i) == 2 
      if DEBUG_Range_i 
        fprintf('debug\n');
      end
      EXCLUSION_FlacToFlac_RANGE = [fluctuation_i*EXCLUSION_FlacToFlac_2_RATE(1),fluctuation_i*EXCLUSION_FlacToFlac_2_RATE(2)];         
      EXCLUSION_ZONE_FlacToFlac_i = exclusion_byFlacToFlac(D_train,subseq_len_i,EXCLUSION_FlacToFlac_RANGE);
      
      Magnitude_ini = MagInfo_ini.magnitude(1:length(D_train));
      if mag_i > 0
        EXCLUSION_MatToMat_RANGE = [mag_i*EXCLUSION_MatToMat_RATE(1),mag_i*EXCLUSION_MatToMat_RATE(2)];
      elseif mag_i < 0
        EXCLUSION_MatToMat_RANGE = [mag_i*EXCLUSION_MatToMat_RATE(2),mag_i*EXCLUSION_MatToMat_RATE(1)];
      end
      EXCLUSION_ZONE_MatToMat_i = exclusion_byMatToMat(D_train,subseq_len_i,Magnitude_ini,EXCLUSION_MatToMat_RANGE,from_i);

      % EXCLUSION_ZONE_Mag_i = EXCLUSION_ZONE_FlacToFlac_i | EXCLUSION_ZONE_MatToMat_i; 
      EXCLUSION_ZONE_Mag_i = EXCLUSION_ZONE_FlacToFlac_i & EXCLUSION_ZONE_MatToMat_i;

      if DEBUG_Range_i 
        fprintf('[%s] [FL,Mat]= %d %d\n',FuncName,sum(EXCLUSION_ZONE_Mag_i),sum(EXCLUSION_ZONE_FlacToFlac_i));
      end
    end

    if abs(type_i) >= 3
      EXCLUSION_FlacToFlac_RANGE = [fluctuation_i*EXCLUSION_FlacToFlac_3_RATE(1),fluctuation_i*EXCLUSION_FlacToFlac_3_RATE(2)];         
      EXCLUSION_ZONE_FlacToFlac_i = exclusion_byFlacToFlac(D_train,subseq_len_i,EXCLUSION_FlacToFlac_RANGE);
      EXCLUSION_ZONE_Mag_i = EXCLUSION_ZONE_FlacToFlac_i;

      %  CS ver.2 
      %{
      Magnitude_ini = MagInfo_ini.magnitude(1:length(D_train));
      if mag_i > 0
        EXCLUSION_MatToMat_RANGE = [mag_i*EXCLUSION_MatToMat_RATE(1),mag_i*EXCLUSION_MatToMat_RATE(2)];
      elseif mag_i < 0
        EXCLUSION_MatToMat_RANGE = [mag_i*EXCLUSION_MatToMat_RATE(2),mag_i*EXCLUSION_MatToMat_RATE(1)];
      end
      EXCLUSION_ZONE_MatToMat_i = exclusion_byMatToMat(D_train,subseq_len_i,Magnitude_ini,EXCLUSION_MatToMat_RANGE,from_i);
      EXCLUSION_ZONE_Mag_i =  EXCLUSION_ZONE_FlacToFlac_i & EXCLUSION_ZONE_MatToMat_i;
      %}
    end

     if abs(type_i) >= 5
      EXCLUSION_FlacToFlac_RANGE = [fluctuation_i*EXCLUSION_FlacToFlac_5_RATE(1),fluctuation_i*EXCLUSION_FlacToFlac_5_RATE(2)];         
      EXCLUSION_ZONE_FlacToFlac_i = exclusion_byFlacToFlac(D_train,subseq_len_i,EXCLUSION_FlacToFlac_RANGE);
      EXCLUSION_ZONE_Mag_i = EXCLUSION_ZONE_FlacToFlac_i;
    end
   
    EXCLUSION_ZONE_i = EXCLUSION_ZONE_i | EXCLUSION_ZONE_Mag_i;

    %% 
    %{
    if abs(type_i) == 1 && exist("EXCLUSION_FlacToFlac_RATE","var") 
      EXCLUSION_FlacToFlac_RANGE = [fluctuation_i*EXCLUSION_FlacToFlac_RATE(1),fluctuation_i*EXCLUSION_FlacToFlac_RATE(2)];         
      EXCLUSION_ZONE_FlacToFlac_i = exclusion_byFlacToFlac(D_train,subseq_len_i,EXCLUSION_FlacToFlac_RANGE);
      EXCLUSION_ZONE_i = EXCLUSION_ZONE_i | EXCLUSION_ZONE_FlacToFlac_i;
    end

    if abs(type_i) == 2 && exist("EXCLUSION_MatToMat_RATE","var")
      Magnitude_ini = MagInfo_ini.magnitude(1:length(D_train));

      if mag_i > 0
        EXCLUSION_MatToMat_RANGE = [mag_i*EXCLUSION_MatToMat_RATE(1),mag_i*EXCLUSION_MatToMat_RATE(2)];
      elseif mag_i < 0
        EXCLUSION_MatToMat_RANGE = [mag_i*EXCLUSION_MatToMat_RATE(2),mag_i*EXCLUSION_MatToMat_RATE(1)];
      end
      EXCLUSION_ZONE_MatToMat_i = exclusion_byMatToMat(D_train,subseq_len_i,Magnitude_ini,EXCLUSION_MatToMat_RANGE,from_i);
      EXCLUSION_ZONE_i = EXCLUSION_ZONE_i | EXCLUSION_ZONE_MatToMat_i;
    end
    %}
    %% find nearest neighbor
    if length(subseq_i) >= 4
      if DIST == "Euc"
        DistanceProfile_i = Spikelet_STOMP_DistanceProfile(D_train, subseq_i, Param);
      elseif DIST == "C22"
        DistanceProfile_i = Spikelet_C22_DistanceProfile(D_train, subseq_i, Param,EXCLUSION_ZONE_i);
        %       [out,Name] = C22_call(y) 
      end
    elseif length(subseq_i) == 3
      fprintf('[%s] [%d %d] subseq length is shorter than 4\n',FuncName,from_i,to_i);  
      if DIST == "Euc"   
        DistanceProfile_i = Spikelet_EUC_DistanceProfile(D_train, subseq_i, Param,EXCLUSION_ZONE_i);
      elseif DIST == "C22"
        DistanceProfile_i = Spikelet_C22_DistanceProfile(D_train, subseq_i, Param,EXCLUSION_ZONE_i);
      end
    else
      fprintf('[%s] [%d %d] subseq length is shorter than 3\n',FuncName,from_i,to_i);  
      DistanceProfile_i = zeros(length(D_train)-length(subseq_i)+1,1);
    end

    if DEBUG_Range_i 
      fprintf('debug\n');
    end
    if exist('EXCLUSION_ZONE_i','var') % && ~isempty(EXCLUSION_ZONE_i)
      DistanceProfile_i(EXCLUSION_ZONE_i) = inf(sum(EXCLUSION_ZONE_i),1); 
    end
    [min_dist_i,min_pos_i] = min(DistanceProfile_i);
    fromNN_i = min_pos_i;
    toNN_i = min_pos_i+subseq_len_i-1; 
    
    if DIST == "Euc"
      min_dist_lengthNorm_i = min_dist_i/sqrt(subseq_len_i);
    elseif DIST == "C22"
     min_dist_lengthNorm_i = min_dist_i;
    end 
    adInfo_i = [spikeId_i,from_i,to_i,type_i,mag_i,min_dist_i,min_dist_lengthNorm_i,fromNN_i,toNN_i,sign_i];
    AdInfo(i,1:length(adInfo_i)) = adInfo_i;
  
    if rem(i, (10*UNIT)) == 0
     fprintf('%d\n',i);
    elseif rem(i, UNIT) == 0
     fprintf('%d ',i);  
    end
  end
  [dist_sorted_val, dist_sorted_id] = sort(AdInfo(:,AdInfo_names=="distLengthNorm"),"descend");
  AdInfo = AdInfo(dist_sorted_id,:);

  %% DEBUG
  DEBUG = true; DEBUG = false;
  if DEBUG
    fig_id = figure(); 
    ROW = 2; COLUMN = 3;  RANK = ROW*COLUMN; distLengthNorm = zeros(RANK,1);
    for i = 1:RANK
      subplot(ROW,COLUMN,i);
      from = AdInfo(i,AdInfo_names=="from"); to = AdInfo(i,AdInfo_names=="to");
      fromNN = AdInfo(i,AdInfo_names=="from_NN"); toNN = AdInfo(i,AdInfo_names=="to_NN");
      T = (fromNN:toNN); plot(D_all(T)); Train_normNN = normalize(D_all(T)); hold on;
      T = (from:to); plot(D_all(T));     Test_norm = normalize(D_all(T)); 
      distLengthNorm(i) = norm(Test_norm-Train_normNN)/sqrt(length(T));
      title(['top ', num2str(i),': ',num2str(length(T)),' [',num2str(from), ' ',num2str(to),'] [',num2str(fromNN),' ',num2str(toNN),'] ',num2str(AdInfo(i,AdInfo_names=="distLengthNorm"))]);  
    end
    disp(distLengthNorm);
  else
    fig_id = [];    
  end
end

function  EXCLUSION_ZONE = exclusion_byFlacToFlac(D_train,subseq_len,EXCLUSION_RANGE_MAT)
  EXCLUSION_ZONE = false(length(D_train) - subseq_len + 1,1);

  % DEBUG = false; 
  for i = 1:length(D_train) - subseq_len + 1  
    % if DEBUG && i == 213; fprintf('debug');end
    subseq_i = D_train(i:i+subseq_len-1);  
    max_i = max(subseq_i);
    min_i = min(subseq_i);
    valueRange_i = max_i - min_i; 
    if valueRange_i < EXCLUSION_RANGE_MAT(1) || EXCLUSION_RANGE_MAT(2) < valueRange_i 
      EXCLUSION_ZONE(i) = true; 
    end
  end
end

function EXCLUSION_ZONE = exclusion_byMatToMat(D_train,subseq_len,Magnitude,EXCLUSION_MatToMat_RANGE,from_cmp)
  EXCLUSION_ZONE = false(length(D_train) - subseq_len + 1,1); 
  Time_train = (1:length(D_train));
  candCenterTimeIndex = (  (EXCLUSION_MatToMat_RANGE(1) <= Magnitude(Time_train)) ...
                         &  (Magnitude(Time_train) <= EXCLUSION_MatToMat_RANGE(2)) );
  for i = 1:length(D_train) - subseq_len + 1
    from_i = i;
    to_i = i + subseq_len - 1;
    candCenterTime_rel = find(candCenterTimeIndex(from_i:to_i));
    candCenterTime = from_i + candCenterTime_rel -1;
    if isempty(candCenterTime) 
      EXCLUSION_ZONE(i) = true;
    else
      for j = 1:length(candCenterTime)  
        center_ij = candCenterTime(j);
        magCenter_ij = Magnitude(center_ij);
        % 
        DEBUG_ij = true; prev_center_id = 0;DEBUG_ij = false;
        if DEBUG_ij && from_cmp == 107364 && center_id ~= prev_center_id
           fprintf('%d ',center_ij);   
           prev_center_id = center_ij;
           figure(); 
             N = 2;
             from_debug_j = max(1,center_ij-N*subseq_len);
             to_debug_j = min(center_ij+N*subseq_len,length(D_train));
             T = (from_debug_j:to_debug_j); 
             plot(T, D_train(T)); hold on; xline(center_ij); title(num2str(center_ij));        
        end
        if magCenter_ij > 0
          leftMag_ij =  max(D_train(center_ij) - D_train(from_i:center_ij));
          rightMag_ij = max(D_train(center_ij) - D_train(center_ij:to_i));
          mag_ij = min(rightMag_ij,leftMag_ij);
        elseif magCenter_ij < 0
          leftMag_ij =  min(D_train(center_ij) - D_train(from_i:center_ij));
          rightMag_ij = min(D_train(center_ij) - D_train(center_ij:to_i));
          mag_ij = max(rightMag_ij,leftMag_ij);
        else
          mag_ij = 0;  
        end
        if mag_ij < EXCLUSION_MatToMat_RANGE(1) || EXCLUSION_MatToMat_RANGE(2) < mag_ij 
          EXCLUSION_ZONE(i) = true;      
        else
          DEBUG_ij2 = false; % if DEBUG_ij2; fprintf('debug_ij2'); end 
        end
      end
    end
  end
end

%{
  
%}
