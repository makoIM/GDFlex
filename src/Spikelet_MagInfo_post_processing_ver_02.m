function MagInfo = Spikelet_MagInfo_post_processing_ver_02(MagInfo, OpName, Param)
   % Spikelet_get_TimeSeriesForm_from_SpikeletDecomposition(MagInfo);  oldversion
   MagInfo.data = Spikelet_MagInfo_getTimeSeriesForm(MagInfo);
   MagInfo.center = find(MagInfo.magnitude ~= 0);
   
   MagInfo.output.(OpName).magnitude = MagInfo.magnitude; 
   MagInfo.output.(OpName).left = MagInfo.left;
   MagInfo.output.(OpName).right = MagInfo.right;
   MagInfo.output.(OpName).data  = MagInfo.data;
   MagInfo.output.(OpName).param_str  = MagInfo.param_str;
   MagInfo.output.(OpName).param  = Param;
   if isfield(MagInfo,'leg_magnitude')
     MagInfo.output.(OpName).leg_magnitude = MagInfo.leg_magnitude; 
   end

end
