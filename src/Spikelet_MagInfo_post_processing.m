function MagInfo = Spikelet_MagInfo_post_processing(MagInfo, OpName)
   if isfield(MagInfo,'output')
      % latest dependent fields 
      MagInfo.center = find(~isnan(MagInfo.type));  
      MagInfo.data = Spikelet_get_TimeSeriesForm_from_SpikeletDecomposition(MagInfo); 
       
      % output fileds (independent)
      MagInfo.output.(OpName).type      = MagInfo.type;
      MagInfo.output.(OpName).value     = MagInfo.value;
      MagInfo.output.(OpName).magnitude = MagInfo.magnitude; 
      MagInfo.output.(OpName).leg_magnitude = MagInfo.leg_magnitude; 
      MagInfo.output.(OpName).left = MagInfo.left;
      MagInfo.output.(OpName).right = MagInfo.right;
      
      if strcmp(OpName,'reduceSpikeByQueryPSR')
         MagInfo.output.(OpName).PSRqueryInfo = MagInfo.PSRqueryInfo;
      else
         if isfield(MagInfo, 'PSRqueryInfo')
            MagInfo = rmfield(MagInfo, 'PSRqueryInfo');
         end
      end

      % dependent
      MagInfo.output.(OpName).data      = MagInfo.data;
      MagInfo.output.(OpName).center    = MagInfo.center;
   else
      MagInfo.operation.(OpName).data      = MagInfo.data;
      MagInfo.operation.(OpName).magnitude = MagInfo.magnitude; 
      MagInfo.operation.(OpName).left = MagInfo.left;
      MagInfo.operation.(OpName).right = MagInfo.right; 
   end
   
end
