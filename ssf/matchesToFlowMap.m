function F = matchesToFlowMap(p_matched,w,h)

  % compute flow map
  Fu = zeros(h,w);
  Fv = zeros(h,w);
  Fi = zeros(h,w);

  iD = sub2ind([h,w],p_matched(2,:),p_matched(1,:));

  Fu(iD) = p_matched(3,:)-p_matched(1,:);
  Fv(iD) = p_matched(4,:)-p_matched(2,:);
  Fi(iD) = 1;

  F = cat(3,Fu,Fv,Fi);

end