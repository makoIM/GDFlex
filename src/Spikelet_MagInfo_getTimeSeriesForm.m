function Data = Spikelet_MagInfo_getTimeSeriesForm(MagInfo)
 
  Data_org = MagInfo.data_org;
  Center = find(MagInfo.magnitude ~= 0);
  Left = MagInfo.left(Center);
  Right = MagInfo.right(Center);
  
  Apx_time = unique([1;Left;Center;Right;length(Data_org)]); 
  if length(Apx_time) >= 2
    Data = interp1(Apx_time, Data_org(Apx_time),(1:length(Data_org))');
  else
    Data =  Data_org;
  end

  if isfield(MagInfo,'constant')
    CSeg = MagInfo.constant;
    Names = MagInfo.constant_names;
    for i = 1:size(CSeg,1)
      from_i = CSeg(i,Names=="left");
      to_i = CSeg(i,Names=="right");
      val_i = CSeg(i,Names=="mean"); 
      Data(from_i:to_i) = val_i *ones(to_i-from_i+1,1); 
    end
  end

  %{
  for i = 1:length(Center)
    from = MagInfo.left(Center(i));
    to = MagInfo.right(Center(i));
    Data(from:to) = Data_org(from:to);
  end
  %}
end