% 第一维维lat，第二维为lon，第三维为时间或者深度,（***把第三维看成深度，对时间循环）
% unit of V is m/s
% u(depth,lat,lon,time)为了函数能够通用以后m函数同一用这种格式
function  [ux,uy]=zh_grad2(u,lon,lat)
%-----------------------------------------------------
% [u,v]=zh_grad2(u,lon,lat)
% horizontal gradient [m,n] --> [m,n]
% assuming [1:m] is northward and positive
% similar to sun_grad2.m, but use different differential method.
%-----------------------------------------------------
% u=squeeze(permute(u,[2,3,1,4]));
[m,n,k]=size(u); dx=ones(m,n)*NaN; dy=dx; ux=u*NaN; uy=ux;
ux(:,1:n-1,:)=u(:,2:n,:)-u(:,1:n-1,:); 
uy(1:m-1,:,:)=u(2:m,:,:)-u(1:m-1,:,:);
dlon(1:n-1)=lon(2:n)-lon(1:n-1);
dlat(1:m-1)=lat(2:m)-lat(1:m-1);

for j=1:m-1,
 for i=1:n-1, 
  dx(j,i)=dlon(3)/abs(dlon(3))*sw_dist([lat(j) lat(j)],[0,dlon(i)],'km')*1e3;%*abs(fg(j));
 end
dy(j,:)=dlat(3)/abs(dlat(3))*sw_dist([0, dlat(j)],[0,0],'km')*1e3;%*abs(fg(j));
end

for i=1:k,
ux(:,:,i)=ux(:,:,i)./dx; uy(:,:,i)=uy(:,:,i)./dy;
end
ux(:,n,:)=ux(:,n-1,:);
uy(m,:,:)=uy(m-1,:,:);

% ux=permute(ux,[3,1,2]);
% uy=permute(uy,[3,1,2]);
