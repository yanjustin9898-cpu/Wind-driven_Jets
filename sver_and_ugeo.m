load('IfremerdailyWindstressAllPacific.mat');
% input
lonp=lons;clear lons
latp=lats;clear lats
rho=1025;
[I,J,K]=size(taux);

load etopo_globe.mat
landmaskp = ones(size(taux(:,:,1)));
for i=1:I
    for j=1:J
        lonind = find(lon>=lonp(i));  
        latind = find(lat<=latp(j));
        if topo(latind(1),lonind(1))<=0
            landmaskp(i,j)=topo(latind(1),lonind(1));
        else
            landmaskp(i,j)=NaN;
        end
    end
end    

A = find(lonp>260);        
B = find(latp>8);
for j=1:length(B)
    nanloc = find(taux(A,B(j),1)~=taux(A,B(j),1));
    taux(A(nanloc(1)):I,B(j),:)=NaN;
    tauy(A(nanloc(1)):I,B(j),:)=NaN;
end


for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    taux(i,nanlocy,:)=NaN;
    tauy(i,nanlocy,:)=NaN;
end

taux = permute(taux,[2,1,3]);
tauy = permute(tauy,[2,1,3]);

taux1=movmean(taux,4,1,'omitnan');
taux=movmean(taux1,4,2,'omitnan');
tauy1=movmean(tauy,4,1,'omitnan');
tauy=movmean(tauy1,4,2,'omitnan');
clear taux1 tauy1 sst1

% partial derivatives of taux and tauy
[txx,txy]=zh_grad2(taux,lonp,latp);
clear txx taux
[tyx,tyy]=zh_grad2(tauy,lonp,latp);
clear tyy tauy
% curl of tau
wscs=tyx-txy;
clear txy tyx
wscs_mean=nanmean(wscs,3);clear wscs

% gradient of tau curl
[wscsx,wscsy]=zh_grad2(wscs_mean,lonp,latp);
clear wscsx
fs = zeros(size(wscsy));
[I,J,K]=size(wscsy);
for i=1:I
    for j=1:J
        fs(i,j,:)=sw_f(latp(i));
    end
end
      
[fsx betas] = zh_grad2(fs,lonp,latp);
clear fs fsx


% contourf(lonp,latp,nanmean(wscs,3),-2e-7:1e-8:2e-7,'LineStyle','None','edgecolor','None');hold on
% cm = redblue(101);
% colorbar();
% colormap(cm);
% worldmap3(2);
% xlim([140 290]);
% ylim([-30 0]);

% Integration to calculatse Sverdrup u
% rho,my_beta,latp,lonp

for i=1:I
    for j=1:J
        dist=sw_dist([latp(i),latp(i)],[lonp(1),lonp(2)],'km')*1000;  % meter
        U_s(i,j)=nansum(wscsy(i,j:J))*dist/(betas(i,j)*rho); %Integration
    end
end
clear betas

U_s(isnan(landmaskp'))=NaN;
Us_mean=nanmean(U_s,3);
%Us_mean=nanmean(U_s,3);
clear U_s;
contourf(lonp,latp,Us_mean,-200:1:200,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
worldmap3(2);
xlim([120 290]);
ylim([-30 30]);
title(strcat('Sverdrup Transport mean from Nov 2006 to Nov 2009'));

Us_mean_lowpassed1=movmean(Us_mean,16,1,'omitnan');
Us_mean_lowpassed2=movmean(Us_mean_lowpassed1,16,2,'omitnan');clear Us_mean_lowpassed1
Us_mean_lowpassed3=movmean(Us_mean-Us_mean_lowpassed2,16,1,'omitnan');
Us_mean_lowpassed4=movmean(Us_mean_lowpassed3,16,2,'omitnan');clear Us_mean_lowpassed3
Us_mean_highpassed=Us_mean-(Us_mean_lowpassed4+Us_mean_lowpassed2);

% wscs_mean_lowpassed1=movmean(wscs_mean,16,1,'omitnan');
% wscs_mean_lowpassed2=movmean(wscs_mean_lowpassed1,16,2,'omitnan');clear wscs_mean_lowpassed1
% wscs_mean_lowpassed3=movmean(wscs_mean-wscs_mean_lowpassed2,16,1,'omitnan');
% wscs_mean_lowpassed4=movmean(wscs_mean_lowpassed3,16,2,'omitnan');clear wscs_mean_lowpassed3
% wscs_mean_highpassed=wscs_mean-(wscs_mean_lowpassed4+wscs_mean_lowpassed2);


Us_mean_highpassed=smoother_2d(Us_mean_highpassed,4);
Us_mean_highpassed(isnan(landmaskp'))=NaN;
Us_mean_lowpassed2(isnan(landmaskp'))=NaN;

plot(-1*nanmean(Us_mean_lowpassed2(:,441:445),2),latp,'LineWidth',1,'color','blue');
hold on;
plot(-1*nanmean(Us_mean_highpassed(:,441:445),2),latp,'LineWidth',1);
grid on;
legend('4^o lowpassed','4^o highpassed')
ylim([-30 30]);
xlim([-100 100]);
xlabel('Zonal Transport m^2/s');
ylabel('Latitude (^oN)');


%zonal average of U_ssverdrup
U_smeridian=nanmean(Us_mean(:,401:521),2);
wscsy_smoothed = smoother_2d(wscsy,4);
wscsy_smeridian=nanmean(wscsy_smoothed(:,401:521),2);
%win = hann(16);
%win=win/sum(win);
%U_slowlow=conv(U_smeridian,win,'same');
U_slowlow=movmean(U_smeridian,16,'omitnan');
U_sbandpassed=nanmean(Us_mean_highpassed(:,401:521),2);


%%
longitude_p=ncread('CMEMS_GeostrophicU_Pacific.nc','longitude');
latitude_p=ncread('CMEMS_GeostrophicU_Pacific.nc','latitude');
u_p=ncread('CMEMS_GeostrophicU_Pacific.nc','ugos');
time=ncread('CMEMS_GeostrophicU_Pacific.nc','time');

longitude_p=longitude_p+360;

[I,J,K]=size(u_p);
u_p1=movmean(u_p,4,1,'omitnan');
u_p=movmean(u_p1,4,2,'omitnan');
clear u_p1


load etopo_globe.mat
landmaskp =size(u_p(:,:,1));
for i=1:I
    for j=1:J
        lonind = find(lon>=longitude_p(i));  
        latind = find(lat<=latitude_p(j));
        if topo(latind(1),lonind(1))<=0
            landmaskp(i,j)=topo(latind(1),lonind(1));
        else
            landmaskp(i,j)=NaN;
        end
    end
end    

for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    u_p(i,nanlocy,:)=NaN;
end

A = find(longitude_p>260);
B = find(latitude_p>0);
for j=1:length(B)
    nanloc = find(u_p(A,B(j),1)~=u_p(A,B(j),1));
    u_p(A(nanloc(1)):I,B(j),:)=0;
end


for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    u_p(i,nanlocy,:)=NaN;
end

u_p=permute(u_p,[2 1 3]);

u_mean=nanmean(u_p,3);clear u_p

u_lowpassed1=movmean(u_mean,16,1,'omitnan');
u_lowpassed2=movmean(u_lowpassed1,16,2,'omitnan');clear u_lowpassed1
u_lowpassed3=movmean(u_mean-u_lowpassed2,16,1,'omitnan');
u_lowpassed4=movmean(u_lowpassed3,16,2,'omitnan');clear u_lowpassed3
up_highpassed=u_mean-(u_lowpassed4+u_lowpassed2);

up_highpassed=smoother_2d(up_highpassed,4);

[J,I,K]=size(up_highpassed);
A = find(longitude_p>260);
B = find(latitude_p>0);
for j=1:length(B)
    nanloc = find(up_highpassed(B(j),A)~=up_highpassed(B(j),A));
    up_highpassed(B(j),A(nanloc(1)):I)=0;
end

for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    up_highpassed(nanlocy,i)=NaN;
end


%%
plot(U_sbandpassed,latp);
hold on;
plot(nanmean(up_highpassed(:,41:81),2)*1000,latitude_p);
ylim([-30 30]);
xlim([-80 80]);grid on


for i=1:80
    plot(up_highpassed(:,i*4),latitude_p);
    hold on;
end

plot(U_sbandpassed,latp);
hold on;
plot(nanmean(up_highpassed(:,101:141),2)*1000,latitude_p);
ylim([-30 30]);
xlim([-80 80]);grid on





subplot(141)
h1=plot(nanmean(Us_mean_highpassed(:,361:401),2),latp);
xlabel('Zonal Sverdrup Transport');
ylabel('Latitude');
hold on;
h2=plot(nanmean(up_highpassed(:,1:41),2)*1000,latitude_p);
ax1 = gca;                          % Get current axes
ax2 = axes('Position', ax1.Position, ...  % Overlay axes
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ...     % Transparent background
           'XColor', 'r', ...       % Red top axis
           'YTick', []);            % Remove y-ticks
ylim([-30 30]);
xlim(ax2,[-80 80]);
xlim(ax1,[-80 80]);grid on
xlabel(ax2,'u_g');
legend([h1, h2], {'U sverdrup', 'u geostrophic'},'Location','northeast');
title('150W-140W');

interp1(latp,nanmean(Us_mean_highpassed(:,361:401),2),latitude_p);
[r,p]=corrcoef(interp1(latp,nanmean(Us_mean_highpassed(:,361:401),2),latitude_p),nanmean(up_highpassed(:,1:41),2))

subplot(142)
h1=plot(nanmean(Us_mean_highpassed(:,441:481),2),latp);
xlabel('Zonal Sverdrup Transport');
ylabel('Latitude');
hold on;
h2=plot(nanmean(up_highpassed(:,81:121),2)*1000,latitude_p);
ax1 = gca;                          % Get current axes
ax2 = axes('Position', ax1.Position, ...  % Overlay axes
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ...     % Transparent background
           'XColor', 'r', ...       % Red top axis
           'YTick', []);            % Remove y-ticks
ylim([-30 30]);
xlim(ax2,[-80 80]);
xlim(ax1,[-80 80]);grid on
xlabel(ax2,'u_g');
legend([h1, h2], {'U sverdrup', 'u geostrophic'},'Location','northeast');
title('130W-120W');

subplot(143)
h1=plot(1.5*nanmean(Us_mean_highpassed(:,521:561),2),latp);
xlabel('Zonal Sverdrup Transport');
ylabel('Latitude');
hold on;
h2=plot(nanmean(up_highpassed(:,161:201),2)*1000,latitude_p);
ax1 = gca;                          % Get current axes
ax2 = axes('Position', ax1.Position, ...  % Overlay axes
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ...     % Transparent background
           'XColor', 'r', ...       % Red top axis
           'YTick', []);            % Remove y-ticks
ylim([-30 30]);
xlim(ax1,[-120 120]);
xlim(ax2,[-120 120]);grid on
xlabel(ax2,'u_g');
legend([h1, h2], {'U sverdrup', 'u geostrophic'},'Location','northeast');
title('110W-100W');

subplot(144)
h1=plot(nanmean(Us_mean_highpassed(:,601:641),2),latp);
xlabel('Zonal Sverdrup Transport');
ylabel('Latitude');
hold on;
h2=plot(nanmean(up_highpassed(:,241:281),2)*500,latitude_p);
ax1 = gca;                          % Get current axes
ax2 = axes('Position', ax1.Position, ...  % Overlay axes
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ...     % Transparent background
           'XColor', 'r', ...       % Red top axis
           'YTick', []);            % Remove y-ticks
ylim([-30 30]);
xlim(ax2,[-60 60]);
xlim(ax1,[-30 30]);grid on
xlabel(ax2,'u_g');
legend([h1, h2], {'U sverdrup', 'u geostrophic'},'Location','northeast');
title('90W-80W');
