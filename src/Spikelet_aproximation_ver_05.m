function MagInfo = Spikelet_aproximation_ver_05(X, Param) 
  FigOn_seg = true; FigOn_seg = false;
  FigOn_leg = true; FigOn_leg = false;
  Debug_On = true; Debug_On = false;
  FuncName = 'Spikelet_aproximation_ver_05';
   
  if ~exist('Param','var')
    Param = [];
  end
  %% Initial spikelet
  start_SAt = tic;
  MagInfo = Spikelet_generateInitialSpikelet_ver_02(X,Inf);
  if isfield(Param,'dataname')
    MagInfo.dataname = Param.dataname;
  else
    MagInfo.dataname = 'unknown';      
  end
  %% size
  InitialOpName = 'generateInitialSpikelet';
  MagInfo.size.original = length(X);
  MagInfo.time.(InitialOpName) = toc(start_SAt);
  MagInfo.size.(InitialOpName) = sum(~isnan(MagInfo.output.(InitialOpName).type));   
  MagInfo.param = Param; 
  fprintf([InitialOpName,': time = ',       num2str(MagInfo.time.(InitialOpName)), ', size= ',num2str(MagInfo.size.(InitialOpName)),'\n']);
  
  if isfield(Param, 'operation')
    MagInfo = Spikelet_Op_call_approximation_ver_03(MagInfo);
  end
  % end_SAt = toc(start_SAt);
end

%% Word_all:  continuation region, Word_maximal = Maximal continuation region
function  [Word_maximal, Word_all] = leg_continuation(NonzeroTime, SupLeft, SupRight)
 Convex_time = [SupLeft, NonzeroTime, SupRight];
 
 Cont_matrix = zeros(length(NonzeroTime),2); 
 Cont_Max = zeros(length(NonzeroTime),2);
 %% find continuation for each convex pattern
 for i = 1:length(SupRight)
   Cont_right_i = find((Convex_time(:,1) == Convex_time(i,2)) & (Convex_time(:,2) == Convex_time(i,3))); 
   if ~isempty(Cont_right_i)
     Cont_matrix(i,Cont_right_i) =  Convex_time(Cont_right_i,3);
   end
 end
 
 %% concatinate continuation pair
 Word_apx_order = zeros(length(NonzeroTime),2); Word_apx_order(:,1) = [1:length(NonzeroTime)];
 Word_apx_time = zeros(length(NonzeroTime),2);
 Check = zeros(length(NonzeroTime),1);
 
 for i = 1:length(NonzeroTime)
   from = i;  
   if Check(i) ~= 1
     Check(i) = 1;    
     [val,pos] = max(Cont_matrix(i,:));
     while val ~= 0 
       Check(pos) = 1;  
       Word_apx_order(from,:) = [from,pos];  
       Cont_Max(from,pos) = SupRight(pos);
       [val,pos] = max(Cont_matrix(pos,:)); 
     end
   end
 end

 Cont_Max_Row = max(Cont_Max,[],2); % maximum for each row
 cont_index = (Cont_Max_Row > 0);
 Cont_Max = Cont_Max(cont_index,:);
 word_num = size(Cont_Max,1);
 
 if word_num >= 1
   Word_all = cell(word_num,1);
   for i = 1:length(Word_all)
     Word_i = struct;
     order_list = find(Cont_Max(i,:) ~= 0);
     Word_i.apx_order = order_list; 
     common_leg = NaN(length(order_list) - 1, 2);
     for j = 1:length(order_list)-1
       common_leg(j,:) = [NonzeroTime(order_list(j)),NonzeroTime(order_list(j+1))];        
     end
     Word_i.common_leg = common_leg;
     Word_i.time_peak_range = [NonzeroTime(order_list(1)),NonzeroTime(order_list(end))];
     Word_i.time_support_range = [SupLeft(order_list(1)),SupRight(order_list(end))];   
   
     Word_all{i} = Word_i;
   end
 else
    word_num = length(NonzeroTime);
    Word_all = cell(word_num,1);
     
    for i = 1:length(Word_all)
     Word_i = struct;
     order_list = [1:length(NonzeroTime)];
     Word_i.apx_order = order_list; 
     Word_i.time_peak_range = [NonzeroTime(order_list(1)),NonzeroTime(order_list(end))];
     Word_i.time_support_range = [SupLeft(order_list(1)),SupRight(order_list(end))];   
     Word_i.common_leg = NaN(0, 2);
     Word_all{i} = Word_i;
    end
 end
 
 maximal_index = true(word_num,1);
 for i = 1:word_num
  tpr_i = Word_all{i}.time_peak_range;   
  for j = 1:word_num
   if j ~= i   
     Inc_ij = range_included(tpr_i, Word_all{j}.time_peak_range);
     if Inc_ij
       maximal_index(i) = false;  
       break; 
     end
   end
  end
 end
 
 Word_maximal = Word_all(maximal_index);
 
end

function Inc = range_included(range1, range2)
  if (range2(1) < range1(1)) & (range1(2) < range2(2))
    Inc = true;
  else
    Inc = false;
  end
end