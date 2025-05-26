function MagInfo = Spikelet_Op_reduceSpikeByMagnitude_ver_02(MagInfo) 
  PLOT = true; PLOT = false;

  %%
  FuncName = 'Spikelet_reduceSpikeByMagnitude_ver_02';
  OpName = 'reduceSpikeByMagnitude';
  
  %% parameters
  Param = MagInfo.param.operation.(OpName);
  Query =  Param.query; 
  Mag = MagInfo.magnitude;
  reduced_time = find(eval(Query));
  
  if ~isempty(reduced_time)
    MagInfo.magnitude(reduced_time) = zeros(length(reduced_time),1);
    MagInfo.left(reduced_time) = reduced_time;
    MagInfo.right(reduced_time) = reduced_time;
    MagInfo.center = find(MagInfo.magnitude ~= 0);
  end

 %% post processing
 MagInfo.param_str = [Param.query];
 MagInfo = Spikelet_MagInfo_post_processing_ver_02(MagInfo, OpName, Param);
 
 %% PLOT
 if exist('PLOT','var') && PLOT
   fig_id = Spikelet_MagInfo_plot_current(MagInfo);
 end
end


