function SegInfo = Spikelet_SegInfo_extractSpikeInSpike_ver2(SegInfo,Param)
  FuncName = 'Spikelet_SegInfo_extractSpikeInSpike_ver2';
  %  ExpandOption = "magnitude_half"; ExpandOption = "LegSupport_half";
  SegDB_names =  ["id","type","from","to","length","center","magnitude","fluctuation"];
  SegDB_exp = []; SegDB_expPlus = []; SegDB_expPlus_names = ["spikeParentId"];

  %% parameter 
  SegStatistics = Param.SegStatistics;
  MaT_pos2 = SegStatistics.positive2.MaT;
  MaT_neg2 = SegStatistics.negative2.MaT;
  % Opp_fluctuationRate = Param.SpikeInSpike.Opp_fluctuationRate;
  %{
  ParamSis = Param.SpikeInSpike;
  Const_MagBandRate = ParamSis.Const_MagBandRate; 
  Const_LengthRate = ParamSis.Const_LengthRate; 
  Opp_MagRate = ParamSis.Opp_MagRate; 
  Opp_BandExpansionRate = ParamSis.Opp_BandExpansionRate; 
  Opp_method = ParamSis.Opp_method;
  %}
  %% MagInfo
  MagInfo = SegInfo.MagInfo;
  Mag = MagInfo.magnitude;
  SegDB = SegInfo.SegDB;
  Names = SegInfo.SegDB_names;
  D = SegInfo.data_org;

  %% main loop
  SpikeOrderList = find(abs(SegDB(:,Names=="type")) == 2);
  id_i = SegInfo.last_id;
  for i = 1:length(SpikeOrderList)
    order_i = SpikeOrderList(i);
    SegDb_i = SegDB(order_i,:); 
    segId_i = SegDB(order_i,Names=="id");  
    from_i = SegDB(order_i,Names=="from"); 
    to_i =  SegDB(order_i,Names=="to"); 
    center_i = SegDB(order_i,Names=="center"); 
    mag_i = SegDB(order_i,Names=="magnitude"); 
    if sign(mag_i) > 0
      MaT2abs_i = MaT_pos2;
    else
      MaT2abs_i = abs(MaT_neg2);
    end
    %{
    67[7014 7021][7000,7050]  76[87549 87682][87009 87025][169422,169460][92239,92340]
    103[7291 7302][8484 8517][3564, 3590][9264 9276]
    98[4965 4994][1634,1655][4190 4199][6782 6805]
    %}
    DEBUG_PV2 = true; Range_debug = [87549 87682]; DEBUG_PV2 = false; 
    if DEBUG_PV2 ...
       && from_i <= Range_debug(1) && Range_debug(2) <= to_i     
       %% && (Range_debug(1) <= from_i && from_i <= Range_debug(2) ...
       %%   || Range_debug(1) <= to_i && to_i <= Range_debug(1)) 
      fprintf('[%s]debug PV2 \n',FuncName);
      figure();
      subplot(2,1,1)
        T = (from_i:to_i); plot(T,D(T)); hold on;
        T_deb = (Range_debug(1):Range_debug(2)); plot(T_deb,D(T_deb),'r','LineWidth',2);
      title([num2str(from_i),' ',num2str(to_i)]);
      subplot(2,1,2)
        normDt = normalize(D(T));
        T = (1:length(T));  plot(T,normDt); 
    end

    %% spike decomposition for z-normalized subsequence
    D_i = D(from_i:to_i);
    Dnorm_i = normalize(D_i);
    MagInfo_i = Spikelet_generateInitialSpikelet_ver_02(Dnorm_i,Inf);

    %% (1) Spike in Spikes: opposite spikes
    % Info_Opp_i = Spikelet_SpikeInSpike_oppositeSpikeSegmentation(MagInfo_i,Param); % GDFlex_ver1
    Info_Opp_i = Spikelet_SpikeInSpike_oppositeSpikeSegmentation_ver2(MagInfo_i,Param,SegDb_i,Names,D);  % GDFlex_ver1_beta
    SegmentList_opp_local = Info_Opp_i.segmentList;
    Names_opp = Info_Opp_i.segmentList_names;
    SegmentList_opp = SegmentList_opp_local;     
    FTCTPpos = [find(Names_opp=="from")...
                ,find(Names_opp=="to")...
                ,find(Names_opp=="center")...
                ,find(Names_opp=="oppTerminal")...
                ,find(Names_opp=="parentSpikeId")];
    for j = 1:size(SegmentList_opp_local,1)
        SegmentList_opp(j,FTCTPpos) = SegmentList_opp_local(j,FTCTPpos) + (from_i -1)*ones(1,length(FTCTPpos));
    end

    for j = 1:size(SegmentList_opp,1)
      id_i = id_i + 1;  
      DEBUG_Opp = true; DEBUG_Opp = false; 
      if DEBUG_Opp && from_i == 4164
        fprintf('[%s] debug opposite spike\n',FuncName);
      end
      from_ij = SegmentList_opp(j,Names_opp=="from");
      to_ij   = SegmentList_opp(j,Names_opp=="to");
      type_ij = SegmentList_opp(j,Names_opp=="type");
      center_ij =  SegmentList_opp(j,Names_opp=="center");
      magnitude_ij = MagInfo.magnitude(center_ij);
      parentSpikeId_ij =  SegmentList_opp(j,Names_opp=="parentSpikeId");
      length_ij = to_ij - from_ij + 1;
      fluctuation_ij = max(D(from_ij:to_ij)) - min(D(from_ij:to_ij)); 
      signTerminalDiff_ij = sign(D(to_ij) - D(from_ij));

      %% SegDB_names =  ["id","type","from","to","length","center","magnitude","fluctuation"];
      % if fluctuation_ij >= MaT2abs_i*Opp_fluctuationRate
        if sign(magnitude_ij) > 0 && magnitude_ij >= MaT_pos2
           type_ij = 2;    
        elseif sign(magnitude_ij) < 0 && magnitude_ij <= MaT_neg2
           type_ij = -2;
        elseif sign(magnitude_ij) > 0 && fluctuation_ij >= MaT_pos2
           type_ij = signTerminalDiff_ij;
        elseif sign(magnitude_ij) > 0 && fluctuation_ij >= abs(MaT_neg2)
           type_ij = signTerminalDiff_ij;
        else
           type_ij = 0; 
        end
        seg_ij =  [id_i,type_ij,from_ij,to_ij,length_ij,center_ij,magnitude_ij,fluctuation_ij];
        [SegDB_exp,SegDB_names,Inserted] = insert_SegDB(SegDB_exp,SegDB_names,seg_ij);
        if Inserted
          SegDB_expPlus = [SegDB_expPlus; seg_ij,parentSpikeId_ij];
        end
      % end
    end

    %% (2) Spike in Spikes: using constant spikes
    Info_Const_i = Spikelet_SpikeInSpike_constantSegmentation(MagInfo_i,Param);
    SegmentList_constant_local = Info_Const_i.segmentList;
    Names_const = Info_Const_i.segmentList_names;
    SegmentList_const = SegmentList_constant_local;
    for j = 1:size(SegmentList_constant_local,1)
        FTCPpos = [find(Names_const=="from"),find(Names_const=="to"),find(Names_const=="center"),find(Names_const=="parentSpikeId")];
        SegmentList_const(j,FTCPpos) = SegmentList_constant_local(j,FTCPpos) + (from_i -1)*ones(1,length(FTCPpos));
    end
    
    for j = 1:size(SegmentList_const,1)
      id_i = id_i + 1;  

      DEBUG_Const = true; DEBUG_Const = false; 
      if DEBUG_Const && id_i == 6685
        fprintf('[%s] debug const\n',FuncName);
      end
      from_ij = SegmentList_const(j,Names_const=="from");
      to_ij   = SegmentList_const(j,Names_const=="to");
      length_ij = to_ij - from_ij + 1;
      center_ij =  SegmentList_const(j,Names_const=="center");
      type_ij = SegmentList_const(j,Names_const=="type");
      parentSpikeId_ij =  SegmentList_const(j,Names_const=="parentSpikeId");
      fluctuation_ij = max(D(from_ij:to_ij)) - min(D(from_ij:to_ij)); 
      magnitude_ij = sign(type_ij)*fluctuation_ij;
     
      mag_left_ij = D(center_ij) - D(from_ij); 
      mag_right_ij = D(center_ij) - D(to_ij);
      magnitudeReal_ij = sign(mag_left_ij)*min(abs(mag_left_ij), abs(mag_right_ij));
      if mag_left_ij > 0 && mag_right_ij > 0 % sign(type_ij) > 0 
        if magnitudeReal_ij > MaT_pos2
          type_ij = 2;  
          magnitude_ij = magnitudeReal_ij;
        end
      elseif  mag_left_ij < 0 && mag_right_ij < 0 % sign(type_ij) < 0
       if magnitude_ij > MaT_neg2
          type_ij = -2;  
          magnitude_ij = magnitudeReal_ij;     
       end
      end
   
      %% SegDB_names =  ["id","type","from","to","length","center","magnitude","fluctuation"];
      seg_ij =  [id_i,type_ij,from_ij,to_ij,length_ij,center_ij,magnitude_ij,fluctuation_ij];
      [SegDB_exp,SegDB_names,Inserted] = insert_SegDB(SegDB_exp,SegDB_names,seg_ij);
      if Inserted
        SegDB_expPlus = [SegDB_expPlus; seg_ij,parentSpikeId_ij];
      end
    end
  end

  %% add SegDB_exp to SegInfo.SegDB 
  SegDB_exp_final = NaN(size(SegDB_exp,1),length(SegInfo.SegDB_names));
  SegDB_exp_final(:,1:size(SegDB_exp,2)) = SegDB_exp; 
  SegInfo.SegDB = [SegInfo.SegDB; SegDB_exp_final];
  SegInfo.last_id = max([SegInfo.SegDB(:,SegDB_names=='id');length(MagInfo.data_org)]);
  SegInfo.SegDB_SpikeInSpike = SegDB_expPlus;
  SegInfo.SegDB_SpikeInSpike_names = [SegDB_names,SegDB_expPlus_names];
end
%% center(peak)
    %{
    The ones with small magnitudes are not suitable as candidates 
    because random bias becomes larger.
    %}
    %{
    id_i = id_i + 1; 
    if ExpandOption == "magnitude_half"
      from_peak_rel_i = find(abs(D(from_i:spikeId_i) - D(spikeId_i)) <= abs(mag_i)/2, 1,'first');
      to_peak_rel__i  = find(abs(D(spikeId_i:to_i) - D(spikeId_i)) <= abs(mag_i)/2, 1,'last');
              
      from_peak_i = from_i + from_peak_rel_i - 1; 
      to_peak_i  =  spikeId_i + to_peak_rel__i -1;
      
    elseif ExpandOption == "LegSupport_half"
      from_peak_i = ceil( (from_i + spikeId_i)/2);
      to_peak_i   = floor( (spikeId_i + to_i)/2);
    end
    type_peak_i = 0;
    length_peak_i = to_peak_i - from_peak_i + 1;
    center_peak_i = floor((from_peak_i + to_peak_i)/2);
    fluctuation_peak_i = max(D(from_peak_i:to_peak_i)) - min(D(from_peak_i:to_peak_i));
    seg_peak_i =  [id_i,type_peak_i,from_peak_i,to_peak_i,length_peak_i,center_peak_i,fluctuation_peak_i];
    if length_peak_i >= 3
      [SegDB,SegDB_names] = insert_SpikeDB(SegDB,SegDB_names,seg_peak_i);
    end
    %}

function  [SegDB_new, SegDB_names,New] = insert_SegDB(SegDB,SegDB_names,seg_i)
  if isempty(SegDB) || size(SegDB,1) == 0
    SegDB_new = seg_i;
    New = true;
    return;
  end

  from_i = seg_i(SegDB_names == "from"); 
  to_i = seg_i(SegDB_names == "to");

  index_i = SegDB(:,SegDB_names == "from") == from_i & SegDB(:,SegDB_names == "to") == to_i;
  if sum(index_i) == 0
    SegDB_new = [SegDB; seg_i];
    New = true;
  else
    SegDB_new = SegDB;
    New = false;
  end
end