function fig_id = Spikelet_MagInfo_plot_current(MagInfo)

     Current = MagInfo.current_operation;
     Current_operatin = MagInfo.param.operation.operation_sequence{Current};
     if Current == 1
       Prev_operation = 'generateInitialSpikelet';
     else
       Prev_operation = MagInfo.param.operation.operation_sequence{Current-1};
     end

     D = MagInfo.data;
     Left = MagInfo.left;
     Right = MagInfo.right;
     XLIM = lim_margin([1 length(D)], 0.025);

     fig_id = figure();  Row = 5;
     subplot(Row,1,1); plot(MagInfo.data_org);        xlim(XLIM); title('original');
     subplot(Row,1,2); plot(MagInfo.output.(Prev_operation).data); xlim(XLIM); 
                       Title = ['prev: ',  Prev_operation]; title(strrep(Title,'_','\_'));
     subplot(Row,1,3); plot(MagInfo.data);            xlim(XLIM); 
       if contains(Current_operatin,'Constant') && isfield(MagInfo,'constant')
         Seg = MagInfo.constant; Names = MagInfo.constant_names;
         hold on; 
         for i = 1:size(Seg,1)
           left = Seg(i,Names=="left");  right = Seg(i,Names=="right");
           range = (left:right);
           plot(range,MagInfo.data(range),'r','LineWidth',1.5);
         end
         Center = find(MagInfo.magnitude ~= 0);
         count_center = 0;
         for i = 1:length(Center)
           for j = 1:size(Seg,1)   
             if Seg(j,Names=="left") <= Center(i) && Center(i) <= Seg(j,Names=="right") 
               count_center = count_center + 1;  
             end
           end
         end
         title(['time series form center-in/const-num= ',num2str(count_center),'/',num2str(size(Seg,1))]);
       else
         title('time series form');
       end
     subplot(Row,1,4); plot(MagInfo.magnitude);       xlim(XLIM); title('magnitude');
     subplot(Row,1,5); 
       D_left = Left - (1:length(Left))';
       plot(D_left);
   
       hold on;
       D_right = Right - (1:length(Right))';
       plot(D_right); xlim(XLIM); title('support');
   
     Title = MagInfo.param_str;
     Fig_figtitle(strrep(Title,'_','\_'));
end

function LIM = lim_margin(D, Rate) 
  Margin = Rate*(max(D) - min(D));
  if Margin ~= 0
    LIM = [min(D)-Margin, max(D)+Margin];
  else
    LIM = [];
  end
end
