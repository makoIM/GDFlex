function MagInfo_zeroCross = Spikelet_MagInfoZC_calculate(D_org)
  Pos_nonzero = find(D_org ~= 0);
  Pos_zero = find(D_org == 0);
  D = D_org(Pos_nonzero);
  D_org_time = (1:length(D_org));

  From = [1:length(D)]';
  To   = D_org_time(Pos_nonzero)';
  Map_center = [From, To];
 
  % find zeroCross
  ZCtime = find(diff(sign(D)) ~= 0);
  
  Center = []; Mag = []; Left = []; Right = [];
  for i = 1:length(ZCtime)-1
    cur = ZCtime(i)+1;  
    next = ZCtime(i+1);

    Region_sameSign_i = D(cur:next);
    % sign_i = sign(D(ceil((cur+next)/2)));
    [max_val_i,max_pos_i] = max(abs(Region_sameSign_i));
    center_i = cur + (max_pos_i - 1);

    Center = [Center; center_i];
    Mag = [Mag; D(center_i)];
    Left = [Left; cur];
    Right = [Right; next];
  end

  Center_final = map_cener(Center,Map_center);
  Mag_final  = Mag;
  Left_final =  map_left(Left,Map_center,Pos_zero);
  Right_final = map_right(Right,Map_center,Pos_zero);

  MagInfo_zeroCross = struct;
  MagInfo_zeroCross.data_org = D_org;
  MagInfo_zeroCross.center = Center_final;
  MagInfo_zeroCross.magnitude = Mag_final;
  MagInfo_zeroCross.left = Left_final;
  MagInfo_zeroCross.right = Right_final;

  %%
  DEBUG = true; DEBUG = false;
  if DEBUG
    LCR = [Left, Center, Right];

    LCR_final = [Left_final, Center_final, Right_final];
    disp([LCR,LCR_final]);
    
    figure();
    plot(D_org,'k'); hold on;
    for i = 1:length(Left_final)
      if rem(i,2) == 0
        Cleft = 'b:'; Cright = 'b--';  Ccenter = 'b:';
      else
        Cleft = 'r:'; Cright = 'r--';  Ccenter = 'r:';
      end
      yline(0,'g--');  
      xline(Left_final(i),Cleft,'LineWidth',1.5);
      xline(Right_final(i),Cright,'LineWidth',1.5);
      xline(Center_final(i),Ccenter);
    end
  end
end

function Center_new = map_cener(Center,Map_center)
   Center_new = NaN(length(Center),1);
   for i = 1:length(Center)
     Center_new(i) = Map_center(Center(i) == Map_center(:,1),2);
   end
end

function Left_new = map_left(Left, Map_center,Pos_zero)
  Left_new = Map_center(Left,2);
  for i = 1:length(Left)
   leftZC_i = Left_new(i)-1;
   if sum(leftZC_i == Pos_zero) >= 1
     Left_new(i) = leftZC_i;   
   end
  end
end

function Right_new = map_right(Right,Map_center,Pos_zero)
  Right_new = Map_center(Right,2);
  for i = 1:length(Right)
    rightZC_i = Right_new(i)+1;
    if sum(rightZC_i  == Pos_zero) >= 1
     Right_new(i) = rightZC_i;   
   end
  end
end