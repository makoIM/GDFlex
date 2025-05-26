function [SpikeDb, Names] = Spikelet_MagInfo_extractSpikeDB(MagInfo)
  Names = ["id","center","magnitude","left","right","support"];

  SpikeDb = NaN(length(MagInfo.center),length(Names));

  for i = 1:length(Names)
    if Names(i) == "id"
     SpikeDb(:,Names==Names(i)) = (1:length(MagInfo.center));        
    elseif Names(i) == "center"
     SpikeDb(:,Names==Names(i)) = MagInfo.center;
    elseif Names(i) == "support"
     SpikeDb(:,Names=="support") ...
       = SpikeDb(:,Names=="right") - SpikeDb(:,Names=="left") + 1;   
    else
     Feature_i = MagInfo.(Names(i));
     SpikeDb(:,Names==Names(i)) = Feature_i(MagInfo.center);
    end
  end

end