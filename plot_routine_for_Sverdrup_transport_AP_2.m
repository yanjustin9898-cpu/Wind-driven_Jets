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

tauxmean=nanmean(taux,3);
tauymean=nanmean(tauy,3);

% partial derivatives of taux and tauy
[txx,txy]=zh_grad2(tauxmean,lonp,latp);
clear txx taux
[tyx,tyy]=zh_grad2(tauymean,lonp,latp);
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
U_smeridian=nanmean(Us_mean(:,361:521),2);
wscsy_smoothed = smoother_2d(wscsy,4);
wscsy_smeridian=nanmean(wscsy_smoothed(:,361:521),2);
%win = hann(16);
%win=win/sum(win);
%U_slowlow=conv(U_smeridian,win,'same');
U_slowlow=zh_filter(U_smeridian,16);
U_sbandpassed=nanmean(Us_mean_highpassed(:,361:521),2);

rms_low=rms(U_slowlow,'omitnan');
rms_high=rms(U_sbandpassed(9:234),'omitnan');

%win = hann(8);
%win=win/sum(win);
%Ulowpassed=conv(Umeridian,win,'same');

transport_bandpassed=zeros(size(U_sbandpassed));
dist=zeros(size(U_sbandpassed));
for i=1:I
        dist(i)=sw_dist([latp(i),latp(i)],[lonp(1),lonp(2)],'km')*1000;  % meter
end

nansum(U_sbandpassed(214:222).*dist(214:222))

ax1 = axes;

h1 = plot(ax1,U_slowlow, latp, 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_sbandpassed, latp, 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ax1.XColor='blue';
ylim([-30 30]);
xlim([-150 150]);
xlabel('Zonal Transport m^2/s');
ylabel('Latitude (^oN)');
legend([h1, h2, h3], {'4^o lowpassed', '4^o highpassed', 'WSC_y'});


plot(U_sbandpassed,latp);
hold on;
plot(transport_band,lat);
ylim([-30 30]);
xlim([-150 150]);

% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match WSC_y line color
hold(ax2,'on');
plot(ax2, 1e12*wscy_smeridian, latp, 'LineWidth', 1, 'color', [0 0.5 0]);
xlabel(ax2, 'WSC_y (m^2/s)');
xlim(ax2, [-1 1]);

%%
%Atlantic
load('IfremerdailyWindstressAtlantic.mat');
% input
lona=lons;clear lons
lata=lats;clear lats
rho=1025;
[I,J,K]=size(taux);

landmaska = ones(size(taux(:,:,1)));
for i=1:I
    for j=1:J
        if lona(i)<0
            lonind = find(lon>=lona(i)+360);  
        else
            lonind = find(lon>=lona(i));  
        end
        latind = find(lat<=lata(j));
        if topo(latind(1),lonind(1))<=0
            landmaska(i,j)=topo(latind(1),lonind(1));
        else
            landmaska(i,j)=NaN;
        end
    end
end    
clear topo

for i =1:I
    nanlocy =  find(landmaska(i,:)~=landmaska(i,:));
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

tauxmean=nanmean(taux,3);
tauymean=nanmean(tauy,3);

% partial derivatives of taux and tauy
[txx,txy]=zh_grad2(taux,lona,lata);
clear txx taux
[tyx,tyy]=zh_grad2(tauy,lona,lata);
clear tyy tauy
% curl of tau
wscs=tyx-txy;
clear txy tyx
wscs_mean=nanmean(wscs,3);clear wscs

% gradient of tau curl
[wscsx,wscsy]=zh_grad2(wscs_mean,lona,lata);
clear wscsx
fs = zeros(size(wscsy));
[I,J,K]=size(wscsy);
for i=1:I
    for j=1:J
        fs(i,j,:)=sw_f(lata(i));
    end
end
      
[fsx betas] = zh_grad2(fs,lona,lata);
clear fs fsx


% contourf(lona,lata,nanmean(wscs,3),-2e-7:1e-8:2e-7,'LineStyle','None','edgecolor','None');hold on
% cm = redblue(101);
% colorbar();
% colormap(cm);
% worldmap3(2);
% xlim([140 290]);
% ylim([-30 0]);

% Integration to calculatse Sverdrup u
% rho,my_beta,lata,lona

for i=1:I
    for j=1:J
        dist=sw_dist([lata(i),lata(i)],[lona(1),lona(2)],'km')*1000;  % meter
        U_a(i,j)=(-1)*nansum(wscsy(i,j:J))*dist/(betas(i,j)*rho); %Integration
    end
end
clear betas

U_a(isnan(landmaska'))=NaN;
Ua_mean=nanmean(U_a,3);

%Us_mean=nanmean(U_s,3);
clear U_a;
contourf(lona,lata,-1*Ua_mean,-200:1:200,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
worldmap3Atl(2);
xlim([-80 10]);
ylim([-30 30]);
title(strcat('Sverdrup Transport mean from Nov 2005 to Nov 2009'));

Ua_mean_lowpassed1=movmean(Ua_mean,16,1,'omitnan');
Ua_mean_lowpassed2=movmean(Ua_mean_lowpassed1,16,2,'omitnan');clear Ua_mean_lowpassed1
Ua_mean_lowpassed3=movmean(Ua_mean-Ua_mean_lowpassed2,16,1,'omitnan');
Ua_mean_lowpassed4=movmean(Ua_mean_lowpassed3,16,2,'omitnan');clear Ua_mean_lowpassed3
Ua_mean_highpassed=Ua_mean-(Ua_mean_lowpassed4+Ua_mean_lowpassed2);

Ua_mean_highpassed=smoother_2d(Ua_mean_highpassed,4);
Ua_mean_highpassed(isnan(landmaska'))=NaN;
Ua_mean_lowpassed2(isnan(landmaska'))=NaN;

A = find(lona<-65);        
B = find(lata<10);
for j=1:length(B)
    nanloc = find(Ua_mean(B(j),A)==Ua_mean(B(j),A));
    test = size(nanloc);
    if test(2)~=0
        Ua_mean(B(j),1:A(nanloc(end)))=0;
        Ua_mean_highpassed(B(j),1:A(nanloc(end)))=0;
        Ua_mean_lowpassed2(B(j),1:A(nanloc(end)))=0;
    end
end

%zonal average of U_ssverdrup
U_ameridian=-1*nanmean(Ua_mean(:,161:321),2);
wscay_smoothed = smoother_2d(wscsy,4);
wscay_smeridian=nanmean(wscay_smoothed(:,161:321),2);
%win = hann(16);
%win=win/sum(win);
%U_slowlow=conv(U_smeridian,win,'same');
U_alowlow=zh_filter(U_ameridian,16);
U_abandpassed=-1*nanmean(Ua_mean_highpassed(:,161:321),2);

rms_low=rms(U_alowlow,'omitnan');
rms_high=rms(U_abandpassed(9:234),'omitnan');

ax1 = axes;

h1 = plot(ax1,U_alowlow, lata, 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_abandpassed, lata, 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ax1.XColor='blue';
ylim([-30 30]);
xlim([-100 100]);
xlabel('Zonal Transport m^2/s');
ylabel('Latitude (^oN)');
legend([h1, h2, h3], {'4^o lowpassed', '4^o highpassed', 'WSC_y'});

% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match WSC_y line color
hold(ax2,'on');
plot(ax2, 1e12*wscay_smeridian, latp, 'LineWidth', 1, 'color', [0 0.5 0]);
xlabel(ax2, 'WSC_y (m^2/s)');
xlim(ax2, [-0.7 0.7]);


%%
%plotting part of 8 subplots
marg_h =[0.06 0.06]; 
marg_w =[0.05 0.05];
Nh = 2;
Nw = 5;
gap=[0.08 0.03];

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-2*gap(2))/2.5;

axh1 = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw1 = (1-sum(marg_w)-2*gap(2))/5;

px = [marg_w(1) marg_w(1)+axw+gap(2) marg_w(1)+axw*2+2*gap(2)];
py = [1-marg_h(1)-axh 1-marg_h(1)-axh*2-gap(1)]; 


ss(1)=subplot('Position',[px(2) py(1) axw axh]);
contourf(lonp,latp,-1*Us_mean_highpassed,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([120 290]);
ylim([-30 30]);
caxis([-80,80]);
xticks([120 150 180 210 240 270]);
xticklabels({'120¡ãE','150¡ãE','180¡ã','150¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('b)','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 12,'TickDir','out'); 
worldmap3(1);

ss(2)=subplot('Position',[px(1) py(1) axw axh]);
contourf(lonp,latp,-1*Us_mean,-200:1:200,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([120 290]);
ylim([-30 30]);
xticks([120 150 180 210 240 270]);
xticklabels({'120¡ãE','150¡ãE','180¡ã','150¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('a)','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 12,'TickDir','out'); 
worldmap3(1);

ss(3)=subplot('Position',[px(3) py(1) axw1 axh1]);
ax1 = gca;

h1 = plot(ax1,U_slowlow,latp, 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_sbandpassed,latp, 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ylim([-30 30]);
xlim([-60 60]);
xticks([-60,-40,-20,0,20,40,60]);
xticklabels({'-60','-40','-20','0','20','40','60'});
xlabel('Zonal Transport m^2/s');
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
set(gca,'FontSize', 12,'TickDir','out'); 
ax1.XLabel.Position(2) = ax1.XLabel.Position(2) + 8.5; 
ax1.XLabel.Position(1) = ax1.XLabel.Position(1) + 0.05; 
ax1.XAxis.TickLabelGapOffset = -6;  % For bottom x-axis (blue)
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match dWSC/dy line color
hold(ax2,'on');
plot(ax2, 1e12*wscsy_smeridian, latp, 'LineWidth', 1, 'color', [0 0.5 0]);
xlabel(ax2, 'dWSC/dy (10^{-12}m^2/s)');
xlim(ax2, [-0.6 0.6]);
ax2.XLabel.Position(2) = ax2.XLabel.Position(2) - 3.85; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1) + 0.09; 
ax2.XAxis.TickLabelGapOffset = -6;  % For top x-axis (orange)
xticks([-0.6,-0.40,-0.20,0,0.20,0.40,0.60]);
xticklabels({'-0.6','-0.4','-0.2','0','0.2','0.4','0.6'});
legend([h1, h2, h3], {'4^o lowpassed', '4^o highpassed', 'dWSC/dy'},'Location','northwest');
title('c) 140¡ã-110¡ãW','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'FontSize', 12,'TickDir','out'); 


ss(4)=subplot('Position',[px(2) py(2) axw axh]);
contourf(lona,lata,-1*Ua_mean_highpassed,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([-70 10]);
ylim([-30 30]);
caxis([-50,50]);
xticks([-80 -50 -20 0]);
xticklabels({'80¡ãW','50¡ãW','20¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('e)','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 12,'TickDir','out'); 
worldmap3Atl(1);


ss(5)=subplot('Position',[px(1) py(2) axw axh]);
contourf(lona,lata,-1*Ua_mean,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([-70 10]);
ylim([-30 30]);
xticks([-80 -50 -20 0]);
xticklabels({'80¡ãW','50¡ãW','20¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('d)','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 12,'TickDir','out'); 
worldmap3Atl(1);

ss(6)=subplot('Position',[px(3) py(2) axw1 axh1]);

ax1 = gca;

h1 = plot(ax1,U_alowlow, lata, 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_abandpassed, lata, 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ylim([-30 30]);
xlim([-50 50]);
xticks([-50,-40,-20,0,20,40,50]);
xlabel('Zonal Transport m^2/s');
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
set(gca,'FontSize', 12,'TickDir','out'); 
ax1.XLabel.Position(2) = ax1.XLabel.Position(2) + 8.8; 
ax1.XLabel.Position(1) = ax1.XLabel.Position(1) + 0.06;
ax1.XAxis.TickLabelGapOffset = -6;  % For bottom x-axis (blue)
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match dWSC/dy line color
hold(ax2,'on');
plot(ax2, 1e12*wscay_smeridian, lata, 'LineWidth', 1, 'color', [0 0.5 0]);
xlabel(ax2, 'dWSC/dy (10^{-12}m^2/s)');
xlim(ax2, [-0.5 0.5]);
ax2.XLabel.Position(2) = ax2.XLabel.Position(2) - 3.85; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1) + 0.06; 
ax2.XAxis.TickLabelGapOffset = -6;  % For top x-axis (orange)
xticks([-0.5,-0.40,-0.20,0,0.20,0.40,0.50]);
xticklabels({'-0.5','-0.4','-0.2','0','0.2','0.4','0.5'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
%legend([h1, h2, h3], {'4^o lowpassed', '4^o highpassed', 'dWSC/dy'},'Location','northwest');
title('f) 40¡ã-10¡ãW','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'FontSize', 12,'TickDir','out'); 

%%
%plotting part of 3*2 subplots
marg_h =[0.08 0.15]; 
marg_w =[0.05 0.05];
Nh = 1;
Nw = 5;
gap=[0.08 0.04];

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-2*gap(2))/2.5;

axh1 = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw1 = (1-sum(marg_w)-2*gap(2))/5;

px = [marg_w(1) marg_w(1)+axw+gap(2) marg_w(1)+axw*2+2*gap(2)];
py = [1-marg_h(1)-axh]; 


ss(1)=subplot('Position',[px(2) py(1) axw axh]);
contourf(lonp,latp,-1*Us_mean_highpassed,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([120 290]);
ylim([-30 30]);
caxis([-80,80]);
xticks([120 150 180 210 240 270]);
xticklabels({'120¡ãE','150¡ãE','180¡ã','150¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({});
%yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('b)','Units','normalized','Position',[0,  1.045],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 22,'TickDir','out'); 
worldmap3(1);

ss(2)=subplot('Position',[px(1) py(1) axw axh]);
contourf(lonp,latp,-1*Us_mean,-200:1:200,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([120 290]);
ylim([-30 30]);
xticks([120 150 180 210 240 270]);
xticklabels({'120¡ãE','150¡ãE','180¡ã','150¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('a)','Units','normalized','Position',[0,  1.045],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 22,'TickDir','out'); 
worldmap3(1);

ss(3)=subplot('Position',[px(3) py(1) axw1 axh1]);
ax1 = gca;

h1 = plot(ax1,U_slowlow,latp, 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_sbandpassed,latp, 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ylim([-30 30]);
xlim([-60 60]);
xticks([-60,-40,-20,0,20,40,60]);
xticklabels({'-60','-40','-20','0','20','40','60'});
xlabel('Zonal Transport m^2/s');
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({});
%yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
set(gca,'FontSize', 22,'TickDir','out'); 
% ax1.XLabel.Position(2) = ax1.XLabel.Position(2) + 8.5; 
% ax1.XLabel.Position(1) = ax1.XLabel.Position(1) + 0.05; 
% ax1.XAxis.TickLabelGapOffset = -6;  % For bottom x-axis (blue)
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match dWSC/dy line color
hold(ax2,'on');
plot(ax2, 1e12*wscsy_smeridian, latp, 'LineWidth', 1, 'color', [0 0.5 0]);
xlabel(ax2, 'dWSC/dy (Nm^{-4})');
xlim(ax2, [-0.6 0.6]);
ax2.XLabel.Position(2) = ax2.XLabel.Position(2) - 0.75; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1) + 0.09; 
ax2.XAxis.TickLabelGapOffset = -5;  % For top x-axis (orange)
xticks([-0.6,-0.40,-0.20,0,0.20,0.40,0.60]);
xticklabels({'-0.6','-0.4','-0.2','0','0.2','0.4','0.6'});
title('c)','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'FontSize', 22,'TickDir','out'); 
legend([h1, h2, h3], {'4^o lowpassed', '4^o highpassed', 'dWSC/dy'},'Location','northwest','FontSize',18,'Box','off');



ss(4)=subplot('Position',[px(2) py(1) axw axh]);
contourf(lona,lata,-1*Ua_mean_highpassed,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([-70 10]);
ylim([-30 30]);
caxis([-50,50]);
xticks([-80 -50 -20 0]);
xticklabels({'80¡ãW','50¡ãW','20¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({});
%yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('b)','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 22,'TickDir','out'); 
worldmap3Atl(1);


ss(5)=subplot('Position',[px(1) py(1) axw axh]);
contourf(lona,lata,-1*Ua_mean,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([-70 10]);
ylim([-30 30]);
xticks([-80 -50 -20 0]);
xticklabels({'80¡ãW','50¡ãW','20¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('a)','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 22,'TickDir','out'); 
worldmap3Atl(1);

ss(6)=subplot('Position',[px(3) py(1) axw1 axh1]);

ax1 = gca;

h1 = plot(ax1,U_alowlow, lata, 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_abandpassed, lata, 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ylim([-30 30]);
xlim([-50 50]);
xticks([-40,-20,0,20,40]);
xlabel('Zonal Transport m^2/s');
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
set(gca,'FontSize', 22,'TickDir','out'); 
% ax1.XLabel.Position(2) = ax1.XLabel.Position(2) + 8.8; 
% ax1.XLabel.Position(1) = ax1.XLabel.Position(1) + 0.06;
% ax1.XAxis.TickLabelGapOffset = -6;  % For bottom x-axis (blue)
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match dWSC/dy line color
hold(ax2,'on');
plot(ax2, 1e12*wscay_smeridian, lata, 'LineWidth', 1, 'color', [0 0.5 0]);
xlabel(ax2, 'dWSC/dy (Nm^{-4})');
xlim(ax2, [-0.5 0.5]);
ax2.XLabel.Position(2) = ax2.XLabel.Position(2) - 0.75; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1) + 0.06; 
ax2.XAxis.TickLabelGapOffset = -5;  % For top x-axis (orange)
xticks([-0.40,-0.20,0,0.20,0.40]);
xticklabels({'-0.4','-0.2','0','0.2','0.4'});
% xticks([-0.5,-0.40,-0.20,0,0.20,0.40,0.50]);
% xticklabels({'-0.5','-0.4','-0.2','0','0.2','0.4','0.5'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({});
%yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('c)','Units','normalized','Position',[0, 1.045],'HorizontalAlignment','left');
set(gca,'FontSize', 22,'TickDir','out'); 
legend([h1, h2, h3], {'4^o lowpassed', '4^o highpassed', 'dWSC/dy'},'Location','northwest','FontSize',18,'Box','off');

%%
marg_h =[0.06 0.06]; 
marg_w =[0.035 0.02];
Nh = 2;
Nw = 5;
gap=[0.07 0.048];

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-2*gap(2))/2.5;

axh1 = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw1 = (1-sum(marg_w)-2*gap(2))/5;

px = [marg_w(1) marg_w(1)+axw+gap(2) marg_w(1)+axw*2+2*gap(2)];
py = [1-marg_h(1)-axh 1-marg_h(1)-axh*2-gap(1)]; 


set(gcf,'units','centimeters','position',[50,-5,40,25]);
ss(1)=subplot('Position',[px(2) py(1) axw axh]);
contourf(lonp,latp,Us_mean_highpassed,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([120 290]);
ylim([-30 30]);
caxis([-80,80]);
xticks([120 160 180 210 240 270]);
xticklabels({'120¡ãE','160¡ãE','180¡ã','160¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('b)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 16,'TickDir','out'); 
ax1=gca;
ax1.XAxis.TickLabelGapOffset = -6;
ax1.YAxis.TickLabelGapOffset = -4;
worldmap3(1);

ss(2)=subplot('Position',[px(1) py(1) axw axh]);
contourf(lonp,latp,Us_mean,-200:1:200,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
caxis([-200 200]);
xlim([120 290]);
ylim([-30 30]);
xticks([120 160 180 210 240 270]);
xticklabels({'120¡ãE','160¡ãE','180¡ã','160¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('a)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 16,'TickDir','out'); 
ax1=gca;
ax1.XAxis.TickLabelGapOffset = -6;
ax1.YAxis.TickLabelGapOffset = -4;
worldmap3(1);

ss(3)=subplot('Position',[px(3) py(1) axw1 axh1]);
ax1 = gca;
h1 = plot(ax1,U_slowlow(16:225),latp(16:225), 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_sbandpassed(16:225),latp(16:225), 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ylim([-30 30]);
xlim([-60 60]);
xticks([-60,-40,-20,0,20,40,60]);
xticklabels({'-60','-40','-20','0','20','40','60'});
%xlabel('Zonal Transport m^2/s');
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
set(gca,'FontSize', 16,'TickDir','out'); 
ax1.XLabel.Position(2) = ax1.XLabel.Position(2) + 4; 
ax1.XLabel.Position(1) = ax1.XLabel.Position(1) + 0.05; 
ax1.XAxis.TickLabelGapOffset = -6;  % For bottom x-axis (blue)
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match dWSC/dy line color
hold(ax2,'on');
plot(ax2, 1e12*wscsy_smeridian(16:225), latp(16:225), 'LineWidth', 1, 'color', [0 0.5 0]);
xlabel(ax2, '$\partial [\mathrm{\nabla}\times\vec{\tau}] / \partial y$ (10$^{-12}$m$^{2}$/s)','Interpreter','latex');
xlim(ax2, [-0.6 0.6]);
ax2.XLabel.Position(2) = ax2.XLabel.Position(2) - 0.9; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1) + 0.09; 
ax2.XAxis.TickLabelGapOffset = -6;  % For top x-axis (orange)
xticks([-0.40,-0.20,0,0.20,0.40]);
xticklabels({'-0.4','-0.2','0','0.2','0.4'});
title('c)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'FontSize', 16,'TickDir','out'); 
lgd=legend([h1, h2, h3], {'4$^{o}$ LP', '4$^{o}$ HP', '$\frac{\partial [\mathrm{\nabla}\times\vec{\tau}]}{ \partial y}$'},'Interpreter','latex','FontSize',16,'Box','off','Position',[0.795, 0.841, 0.1, 0.1],'Units','normalized');


ss(4)=subplot('Position',[px(2) py(2) axw axh]);
contourf(lona,lata,-1*Ua_mean_highpassed,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([-70 10]);
ylim([-30 30]);
caxis([-50,50]);
xticks([-80 -50 -20 0]);
xticklabels({'80¡ãW','50¡ãW','20¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('e)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 16,'TickDir','out'); 
ax1=gca;
ax1.XAxis.TickLabelGapOffset = -6;
ax1.YAxis.TickLabelGapOffset = -4;
worldmap3Atl(1);


ss(5)=subplot('Position',[px(1) py(2) axw axh]);
contourf(lona,lata,-1*Ua_mean,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
caxis([-100 100]);
xlim([-70 10]);
ylim([-30 30]);
xticks([-80 -50 -20 0]);
xticklabels({'80¡ãW','50¡ãW','20¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('d)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 16,'TickDir','out'); 
ax1=gca;
ax1.XAxis.TickLabelGapOffset = -6;
ax1.YAxis.TickLabelGapOffset = -4;
worldmap3Atl(1);

ss(6)=subplot('Position',[px(3) py(2) axw1 axh1]);

ax1 = gca;
h1 = plot(ax1,U_alowlow(16:225), lata(16:225), 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_abandpassed(16:225), lata(16:225), 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ylim([-30 30]);
xlim([-25 25]);
xticks([-20,-10,0,10,20]);
xlabel('Zonal Transport m^2/s');
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
set(gca,'FontSize', 16,'TickDir','out'); 
ax1.XLabel.Position(2) = ax1.XLabel.Position(2) + 4; 
ax1.XLabel.Position(1) = ax1.XLabel.Position(1) + 0.06;
ax1.XAxis.TickLabelGapOffset = -6;  % For bottom x-axis (blue)
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match dWSC/dy line color
hold(ax2,'on');
plot(ax2, 1e12*wscay_smeridian(16:225), lata(16:225), 'LineWidth', 1, 'color', [0 0.5 0]);
%xlabel(ax2, '$\partial [\mathrm{\nabla}\times\vec{\tau}] / \partial y$ (10$^{-12}$m$^{2}$/s)','Interpreter','latex');
xlim(ax2, [-0.27 0.27]);
ax2.XLabel.Position(2) = ax2.XLabel.Position(2) - 1.3; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1) + 0.06; 
ax2.XAxis.TickLabelGapOffset = -6;  % For top x-axis (orange)
xticks([-0.20,-0.1,0,0.10,0.20]);
xticklabels({'-0.2','-0.1','0','0.1','0.2'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('f)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'FontSize', 16,'TickDir','out'); 

%%
%cutoff at 25deg
marg_h =[0.06 0.06]; 
marg_w =[0.035 0.02];
Nh = 2;
Nw = 5;
gap=[0.07 0.048];

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw = (1-sum(marg_w)-2*gap(2))/2.5;

axh1 = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
axw1 = (1-sum(marg_w)-2*gap(2))/5;

px = [marg_w(1) marg_w(1)+axw+gap(2) marg_w(1)+axw*2+2*gap(2)];
py = [1-marg_h(1)-axh 1-marg_h(1)-axh*2-gap(1)]; 


set(gcf,'units','centimeters','position',[50,-5,40,25]);
ss(1)=subplot('Position',[px(2) py(1) axw axh]);
contourf(lonp,latp,Us_mean_highpassed,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([120 290]);
ylim([-25 25]);
caxis([-80,80]);
xticks([120 160 180 210 240 270]);
xticklabels({'120¡ãE','160¡ãE','180¡ã','160¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('b)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 16,'TickDir','out'); 
ax1=gca;
ax1.XAxis.TickLabelGapOffset = -6;
ax1.YAxis.TickLabelGapOffset = -4;
worldmap3(1);

ss(2)=subplot('Position',[px(1) py(1) axw axh]);
contourf(lonp,latp,Us_mean,-200:1:200,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
caxis([-200 200]);
xlim([120 290]);
ylim([-25 25]);
xticks([120 160 180 210 240 270]);
xticklabels({'120¡ãE','160¡ãE','180¡ã','160¡ãW','120¡ãW','90¡ãW'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('a)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 16,'TickDir','out'); 
ax1=gca;
ax1.XAxis.TickLabelGapOffset = -6;
ax1.YAxis.TickLabelGapOffset = -4;
worldmap3(1);

ss(3)=subplot('Position',[px(3) py(1) axw1 axh1]);
ax1 = gca;
h1 = plot(ax1,U_slowlow,latp, 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_sbandpassed,latp, 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ylim([-25 25]);
xlim([-60 60]);
xticks([-60,-40,-20,0,20,40,60]);
xticklabels({'-60','-40','-20','0','20','40','60'});
%xlabel('Zonal Transport m^2/s');
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
set(gca,'FontSize', 16,'TickDir','out'); 
ax1.XLabel.Position(2) = ax1.XLabel.Position(2) + 4; 
ax1.XLabel.Position(1) = ax1.XLabel.Position(1) + 0.05; 
ax1.XAxis.TickLabelGapOffset = -6;  % For bottom x-axis (blue)
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match dWSC/dy line color
hold(ax2,'on');
plot(ax2, 1e12*wscsy_smeridian, latp, 'LineWidth', 1, 'color', [0 0.5 0]);
xlabel(ax2, '$\partial [\mathrm{\nabla}\times\vec{\tau}] / \partial y$ (10$^{-12}$m$^{2}$/s)','Interpreter','latex');
xlim(ax2, [-0.6 0.6]);
ylim(ax2, [-25 25]);
ax2.XLabel.Position(2) = ax2.XLabel.Position(2) - 0.5; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1) + 0.09; 
ax2.XAxis.TickLabelGapOffset = -6;  % For top x-axis (orange)
xticks([-0.40,-0.20,0,0.20,0.40]);
xticklabels({'-0.4','-0.2','0','0.2','0.4'});
title('c)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'FontSize', 16,'TickDir','out'); 
lgd=legend([h1, h2, h3], {'4$^{o}$ LP', '4$^{o}$ HP', '$\frac{\partial [\mathrm{\nabla}\times\vec{\tau}]}{ \partial y}$'},'Interpreter','latex','FontSize',16,'Box','off','Position',[0.795, 0.841, 0.1, 0.1],'Units','normalized');


ss(4)=subplot('Position',[px(2) py(2) axw axh]);
contourf(lona,lata,-1*Ua_mean_highpassed,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
xlim([-70 10]);
ylim([-25 25]);
caxis([-50,50]);
xticks([-80 -50 -20 0]);
xticklabels({'80¡ãW','50¡ãW','20¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('e)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 16,'TickDir','out'); 
ax1=gca;
ax1.XAxis.TickLabelGapOffset = -6;
ax1.YAxis.TickLabelGapOffset = -4;
worldmap3Atl(1);


ss(5)=subplot('Position',[px(1) py(2) axw axh]);
contourf(lona,lata,-1*Ua_mean,-100:1:100,'LineStyle','None','edgecolor','None');hold on
cm = redblue(101);
colorbar();
colormap(cm);
caxis([-100 100]);
xlim([-70 10]);
ylim([-25 25]);
xticks([-80 -50 -20 0]);
xticklabels({'80¡ãW','50¡ãW','20¡ãW','0¡ã'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30¡ãS','20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN','30¡ãN'});
title('d)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'Color',0.6*[1 1 1],'FontSize', 16,'TickDir','out'); 
ax1=gca;
ax1.XAxis.TickLabelGapOffset = -6;
ax1.YAxis.TickLabelGapOffset = -4;
worldmap3Atl(1);

ss(6)=subplot('Position',[px(3) py(2) axw1 axh1]);

ax1 = gca;
h1 = plot(ax1,U_alowlow(16:225), lata(16:225), 'LineWidth', 1, 'color', 'blue');
hold on;
h2 = plot(ax1,U_abandpassed(16:225), lata(16:225), 'LineWidth', 1,'color', 'red');
% Dummy plot for legend entry (invisible on ax1)
h3 = plot(ax1,nan, nan, 'LineWidth', 1, 'color', [0 0.5 0]);
grid on;
ylim([-25 25]);
xlim([-25 25]);
xticks([-20,-10,0,10,20]);
xlabel('Zonal Transport m^2/s');
yticks([flip(-1*[0 10 20 ]) [10 20 ]]);
yticklabels({'20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN'});
set(gca,'FontSize', 16,'TickDir','out'); 
ax1.XLabel.Position(2) = ax1.XLabel.Position(2) + 3.5; 
ax1.XLabel.Position(1) = ax1.XLabel.Position(1) + 0.06;
ax1.XAxis.TickLabelGapOffset = -6;  % For bottom x-axis (blue)
% Secondary axes for wscy_smeridian (top x-axis)
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'Color', 'none', ...
           'YColor', 'none', ...
           'XColor', [0 0.5 0]); % Match dWSC/dy line color
hold(ax2,'on');
plot(ax2, 1e12*wscay_smeridian(16:225), lata(16:225), 'LineWidth', 1, 'color', [0 0.5 0]);
%xlabel(ax2, '$\partial [\mathrm{\nabla}\times\vec{\tau}] / \partial y$ (10$^{-12}$m$^{2}$/s)','Interpreter','latex');
xlim(ax2, [-0.27 0.27]);
ylim(ax2, [-25 25]);
ax2.XLabel.Position(2) = ax2.XLabel.Position(2) - 1.3; 
ax2.XLabel.Position(1) = ax2.XLabel.Position(1) + 0.06; 
ax2.XAxis.TickLabelGapOffset = -6;  % For top x-axis (orange)
xticks([-0.20,-0.1,0,0.10,0.20]);
xticklabels({'-0.2','-0.1','0','0.1','0.2'});
yticks([flip(-1*[0 10 20]) [10 20]]);
yticklabels({'20¡ãS','10¡ãS','0¡ã','10¡ãN','20¡ãN'});
title('f)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
set(gca,'FontSize', 16,'TickDir','out'); 

