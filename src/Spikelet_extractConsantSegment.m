function [Cmax_sorted,Names] = Spikelet_extractConsantSegment(D,CLT,CMT,LimitLen)
  % extract all constant segment
  [Clen, Names] = Spikelet_extractConsantSegment_all(D,CLT,CMT,LimitLen);

  %{
  DEBUG = true; DEBUG = false; 
  if DEBUG;figure(); subplot(2,1,1); plot_Clen(D,Clen,Names); title(['Spikelet-extractConsantSegment-all size=',num2str(size(Clen,1)),]); end
  %}

  % select maximal segment
  Clen = maximalClen(Clen,Names);
  Clen(:,Names=='id') = (1:size(Clen,1))';

  % select miminum std segment
  Len_unique = unique(Clen(:,Names=="support"));
  Len_unique = sort(Len_unique,"descend");
  
  Cmax = []; Clen_RemainIndex = true(size(Clen,1),1);
  for i = 1:length(Len_unique)
    len_i = Len_unique(i);
    %{
    if len_i == 313; 
        fprintf('debug\n');
    end
    %}
    Clen_i_index = (Clen(:,Names=="support")==len_i & Clen_RemainIndex);
    Clen_i = Clen(Clen_i_index,:);
    for j = 1:size(Clen_i,1)
     Clen_ij = Clen_i(j,:); 
     if Clen_RemainIndex(Clen_ij(Names=="id"))
       %{  
       if ismember(Clen_ij(Names=="id"),[1,2])
          fprintf('id %d\n',Clen_i(j,1));
       end
       %}
       Clen_ij_interact = interacting_list(Clen_ij,Clen(Clen_RemainIndex,:),Names);
       [Cmax_i, Clen_RemainIndex] = longest_and_minimumStd(Clen_ij_interact,Names,Clen_RemainIndex);
       Cmax = [Cmax; Cmax_i];
     end
    end
  end

  % sort Cmax
  [Left_sorted_val, Left_sorted_index]  = sort(Cmax(:,Names=="left"), "ascend");
  Cmax_sorted = Cmax(Left_sorted_index,:);
  Cmax_sorted(:,Names=="id") = (1:size(Cmax_sorted,1))';

  %{
  if DEBUG;subplot(2,1,2); plot_Clen(D,Cmax,Names); title(['Spikelet-extractConsantSegment size=',num2str(size(Cmax,1))]); end
  %}
end

function Clen_new = maximalClen(Clen,Names)
  Clen_1 = maximalClenFromLeft(Clen,Names);
  Clen_new = maximalClenFromRight(Clen_1,Names);
end

function Clen_new = maximalClenFromLeft(Clen,Names)
  [Left_val,Left_pos] = sort(Clen(:,Names=="left"));
  Left_uniq = unique(Left_val);
  Flag = true(size(Clen,1),1);
  for i = 1:length(Left_uniq)
    left_i = Left_uniq(i);  
    % if left_i == 930; fprintf('debug [%d]\n',left_i);end
    Index_i = find(Clen(:,Names=="left") == left_i);
    right_ij = Clen(Index_i,Names=="right");
    [right_max, max_pos_rel] = max(right_ij);  % The right_max is unique even though multiple max_pos_rel 
    maximal_index = Index_i(max_pos_rel);
    del_i = setdiff(Index_i,maximal_index);
    if ~isempty(del_i)
     Flag(del_i) = false(length(del_i),1);
    end
  end

  Clen_new = Clen(Flag,:);
end

function Clen_new = maximalClenFromRight(Clen,Names)
  [Right_val,Right_pos] = sort(Clen(:,Names=="right"));
  Right_uniq = unique(Right_val);
  Flag = true(size(Clen,1),1);
  for i = 1:length(Right_uniq)
    right_i = Right_uniq(i);  
    Index_i = find(Clen(:,Names=="right") == right_i);
    left_ij = Clen(Index_i,Names=="left");
    [left_min, min_pos_rel] = min(left_ij);  % The right_max is unique even though multiple max_pos_rel 
    minmal_index = Index_i(min_pos_rel);
    del_i = setdiff(Index_i,minmal_index);
    if ~isempty(del_i)
     Flag(del_i) = false(length(del_i),1);
    end
  end

  Clen_new = Clen(Flag,:);
end

function plot_Clen(D,Clen,Names)
  plot(D); hold on;
  for i = 1:size(Clen,1)
    val = Clen(i,Names=="mean");
    left = Clen(i,Names=="left");
    right = Clen(i,Names=="right");
    plot(left:right,val*ones(right-left+1,1),'r');
  end
end

function [Clen_max, Cmax_RemainIndex] = longest_and_minimumStd(Clen,Names,Cmax_RemainIndex)
  % if ~isempty(Clen) && sum(Clen(:,Names=="left") == 930) ~= 0; fprintf('debug\n'); end
  [max_supp,max_supp_pos] = max(Clen(:,Names=='support'));
  max_index = find(Clen(:,Names=="support") == max_supp);

  if sum(max_index) >= 2
    Clen_max_cand = Clen(max_index,:);
    [min_std, min_std_pos_rel] = min(Clen_max_cand(:,Names=="std"));
    maxAndminStd_id = max_index(min_std_pos_rel);
    Clen_max = Clen(maxAndminStd_id,:);
  else
    Clen_max = Clen(max_index,:);
  end
  Cmax_checked = Clen(:,Names=="id");
  Cmax_RemainIndex(Cmax_checked) = false(length(Cmax_checked),1);
  %{
  Cmax_checked = setdiff(Clen(:,Names=="id"), Clen_max(Names=="id"));
  if ~isempty(Cmax_checked)
    Cmax_RemainIndex(Cmax_checked) = false(length(Cmax_checked),1);
  end
  %}
end

function   [min_id,int_id] = min_std_in_intersectingRegion(Clen,Names,i)
   left_i = Clen(i,Names=="left");
   right_i = Clen(i,Names=="right"); 
   std_i = Clen(i,Names=="std");

   min_id = i; int_id = [i]; min_std = std_i;
   for j = 1:size(Clen,1)
     if j ~= i  
      left_j = Clen(j,Names=="left");
      right_j = Clen(j,Names=="right"); 
      if (left_j <= left_i  && left_i <= right_j) || (left_i <= left_j  && left_j <= right_i)
        std_j = Clen(j,Names=="std");
        if std_j < min_std
          min_id = j;
          min_std = std_j;
        end
      end
     end
   end
end

function Clen_cand = interacting_list(Clen_i,Clen,Names)
  Clen_cand = [];
  left_i = Clen_i(Names=="left"); 
  right_i = Clen_i(Names=="right");
  for j = 1:size(Clen,1)
    left_j = Clen(j,Names=="left");
    right_j = Clen(j,Names=="right");
    if (left_j <= left_i && left_i <= right_j) || (left_j <= right_i && right_i <= right_j) ...
      || (left_i <= left_j && right_j <= right_i) 
      % The following case should not occur since it is processed from a long subsequence  
      % || (left_j <= left_i && right_i <= right_j) 
      Clen_cand = [Clen_cand; Clen(j,:)];
    end
  end
end


function [Clen, Names] = Spikelet_extractConsantSegment_all(D,CLT,CMT,LimitLen)
  Names = ["id","left","right","magnitude","support","mean","std"];
  Clen = []; count = 0;
  for i = 1:length(D)
    % if i == 2662; fprintf('[%s] debug\n',FuncName);  end
    right = i + search_right(i, D, CMT, LimitLen) - 1;
    left = i - search_left(i, D, CMT, LimitLen);
    support  = right - left + 1;
    magnitude = 0;
    if support >= CLT ...
            && ( isempty(Clen) || (sum(Clen(:,Names=="left") == left & Clen(:,Names=="right") == right) == 0) )
       count = count + 1; 
       Clen = [Clen; count, left, right,magnitude,support,mean(D(left:right)),std(D(left:right))];
    end
  end
end

function count = search_right(pos_i,D,MagThr, LimitLen)
  MagThr_half = MagThr/2;
  count = 0;
  for i = pos_i:length(D)
      if abs(D(i) - D(pos_i)) <= MagThr_half
        count = count + 1;
        if count >= LimitLen
          break;
        end
      else
        break;  
      end
  end
end

function count = search_left(pos_i,D,MagThr, LimitLen)
  MagThr_half = MagThr/2;
  count = 0;
  for i = pos_i:-1:1
      if abs(D(i) - D(pos_i)) <= MagThr_half
        count = count + 1;
        if count >= LimitLen
          break;
        end
      else
        break;  
      end
  end
end
