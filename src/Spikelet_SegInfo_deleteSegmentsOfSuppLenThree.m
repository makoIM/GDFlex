function [SegInfo,RemovedId] = Spikelet_SegInfo_deleteSegmentsOfSuppLenThree(SegInfo,Param)
  % argument
  % MmedianSupp = Param.SuppCut_Mmedian;   
  %MagInfo = SegInfo.MagInfo;  
  D = SegInfo.data_org;
  SegDB = SegInfo.SegDB;  SegDB_org = SegDB;
  Names = SegInfo.SegDB_names;

  Order_pos = find(SegDB(:,Names=="length") == 3 & SegDB(:,Names=="type") == 2);
  RemovedId_pos = removeExceptForMax(SegDB,Names,Order_pos);

  Id_neg = find(SegDB(:,Names=="length") == 3 & SegDB(:,Names=="type") == -2);
  RemovedId_neg = removeExceptForMax(SegDB,Names,Id_neg);

  % output
  RemovedId = [RemovedId_pos; RemovedId_neg]; 
  if isempty(RemovedId)
  else
    OrderId_new = setdiff((1:size(SegDB,1)), RemovedId);
    SegDB = SegDB(OrderId_new,:);
  end
  SegInfo.SegDB = SegDB;
  SegInfo.last_id = max(length(D),max(SegDB(:,Names=="id")));

end

function RemovedId = removeExceptForMax(SegDB,Names,Order_pos)
  if isempty(Order_pos)
    RemovedId = [];
  else
    [maxVal,maxOrderPosRef] = max(abs(SegDB(Order_pos, Names=="magnitude")));
    order_max_pos = Order_pos(maxOrderPosRef);
    RemovedId = setdiff(Order_pos, order_max_pos);
  end
end