function Data = Spikelet_get_TimeSeriesForm_from_SpikeletDecomposition(MagInfo)
 
  Data_org = MagInfo.data_org;
  Center = find(~isnan(MagInfo.type));
  Left = MagInfo.left(Center);
  Right = MagInfo.right(Center);
  
  Apx_time = unique([1;Left;Center;Right;length(Data_org)]); 
  if length(Apx_time) >= 2
    Data = interp1(Apx_time, Data_org(Apx_time),(1:length(Data_org))');
  else
    Data =  Data_org; 
  end
  Center_const =  find(MagInfo.type == 0);
  for i = 1:length(Center_const)
    from = MagInfo.left(Center_const(i));
    to = MagInfo.right(Center_const(i));
    value = MagInfo.value(Center_const(i));
    Data(from:to) = value*ones(to-from+1,1);
  end
end