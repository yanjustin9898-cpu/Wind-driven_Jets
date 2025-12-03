load('IfremerdailyTAUandSSTAllPacific.mat');

taux(isnan(sst))=NaN;
tauy(isnan(sst))=NaN;
sst=sst-273.15;

taux1=movmean(taux,4,1,'omitnan');
taux=movmean(taux1,4,2,'omitnan');
tauy1=movmean(tauy,4,1,'omitnan');
tauy=movmean(tauy1,4,2,'omitnan');
sst1=movmean(sst,4,1,'omitnan');
sst=movmean(sst1,4,2,'omitnan');
clear taux1 tauy1 sst1

[I,J,K]=size(taux);
A = find(lons>260);
B = find(lats>0);
for j=1:length(B)
    nanloc = find(taux(A,B(j),1)~=taux(A,B(j),1));
    taux(A(nanloc(1)):I,B(j),:)=NaN;
    tauy(A(nanloc(1)):I,B(j),:)=NaN;
    sst(A(nanloc(1)):I,B(j),:)=NaN;
end

load etopo_globe.mat
landmaskp =size(taux(:,:,1));
for i=1:I
    for j=1:J
        lonind = find(lon>=lons(i));  
        latind = find(lat<=lats(j));
        if topo(latind(1),lonind(1))<=0
            landmaskp(i,j)=topo(latind(1),lonind(1));
        else
            landmaskp(i,j)=NaN;
        end
    end
end    
clear topo

for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    taux(i,nanlocy,:)=NaN;
    tauy(i,nanlocy,:)=NaN;
    sst(i,nanlocy,:)=NaN;
end


taux = permute(taux,[2,1,3]);
tauy = permute(tauy,[2,1,3]);
sst = permute(sst,[2,1,3]);
tauxmean=nanmean(taux,3);
tauymean=nanmean(tauy,3);

% partial derivatives of taux and tauy
[txx,txy]=zh_grad2(tauxmean,lons,lats);
clear txx 
[tyx,tyy]=zh_grad2(tauymean,lons,lats);
clear tyy 
% curl of tau
wsc=tyx-txy;
clear txy tyx

sstmean=nanmean(sst,3);
wscmeans=nanmean(wsc,3);
% gradient of sst
[sstx,ssty]=zh_grad2(sst,lons,lats);

[sstxx,sstxy]=zh_grad2(sstx,lons,lats);
clear sstxy;
[sstyx,sstyy]=zh_grad2(ssty,lons,lats);
clear sstyx;
% SST laplacian
sstlap=sstxx+sstyy;clear sstxx sstyy
sstlapmean_p=nanmean(sstlap,3);clear sstlap

sstlapmean_lowpassed1=movmean(sstlapmean_p,16,1,'omitnan');
sstlapmean_lowpassed2=movmean(sstlapmean_lowpassed1,16,2,'omitnan');clear sstlapmean_lowpassed1
sstlapmean_lowpassed3=movmean(sstlapmean_p-sstlapmean_lowpassed2,16,1,'omitnan');
sstlapmean_lowpassed4=movmean(sstlapmean_lowpassed3,16,2,'omitnan');clear sstlapmean_lowpassed3
sstlapmeanp_highpassed=sstlapmean_p-(sstlapmean_lowpassed4+sstlapmean_lowpassed2);clear sstlapmean_lowpassed2 sstlapmean_lowpassed4

taudir=atan2(tauy,taux);
taumag=sqrt(taux.^2+tauy.^2);clear taux tauy
% sst gradient direction and magnitude
sstgraddir=atan2(ssty,sstx);
sstgradmag=sqrt(sstx.^2+ssty.^2);clear sstx ssty
%crosswind sst gradient
crosswindsstgrad = sstgradmag.*sin(taudir-sstgraddir);
cwsgs_mean=nanmean(crosswindsstgrad,3);clear crosswindsstgrad

wscmeans_lowpassed1=movmean(wscmeans,16,1,'omitnan');
wscmeans_lowpassed2=movmean(wscmeans_lowpassed1,16,2,'omitnan');clear wscmeans_lowpassed1
wscmeans_lowpassed3=movmean(wscmeans-wscmeans_lowpassed2,16,1,'omitnan');
wscmeans_lowpassed4=movmean(wscmeans_lowpassed3,16,2,'omitnan');clear wscmeans_lowpassed3
wscmeans_highpassed=wscmeans-(wscmeans_lowpassed4+wscmeans_lowpassed2);


sstmeans_lowpassed1=movmean(sstmean,16,1,'omitnan');
sstmeans_lowpassed2=movmean(sstmeans_lowpassed1,16,2,'omitnan');clear sstmeans_lowpassed1
sstmeans_lowpassed3=movmean(sstmean-sstmeans_lowpassed2,16,1,'omitnan');
sstmeans_lowpassed4=movmean(sstmeans_lowpassed3,16,2,'omitnan');clear sstmeans_lowpassed3
sstmeans_highpassed=sstmean-(sstmeans_lowpassed4+sstmeans_lowpassed2);



cwsgsmean_lowpassed1=movmean(cwsgs_mean,16,1,'omitnan');
cwsgsmean_lowpassed2=movmean(cwsgsmean_lowpassed1,16,2,'omitnan');
cwsgsmean_lowpassed3=movmean(cwsgs_mean-cwsgsmean_lowpassed2,16,1,'omitnan');
cwsgsmean_lowpassed4=movmean(cwsgsmean_lowpassed3,16,2,'omitnan');
cwsgsmean_highpassed=cwsgs_mean-(cwsgsmean_lowpassed4+cwsgsmean_lowpassed2);

wscmeans_1d=reshape(wscmeans_highpassed,[],1);
cwsgsmean_1d=reshape(cwsgsmean_highpassed,[],1);
sstlapmeanp_1d=reshape(sstlapmeanp_highpassed,[],1);

[J,I,K]=size(cwsgsmean_highpassed);
A = find(lons>260);
B = find(lats>0);
for j=1:length(B)
    nanloc = find(cwsgsmean_highpassed(B(j),A)~=cwsgsmean_highpassed(B(j),A));
    cwsgsmean_highpassed(B(j),A(nanloc(1)):I)=0;
    wscmeans_highpassed(B(j),A(nanloc(1)):I)=0;
    sstlapmeanp_highpassed(B(j),A(nanloc(1)):I)=0;
end

for i =1:I
    nanlocy =  find(landmaskp(i,:)~=landmaskp(i,:));
    cwsgsmean_highpassed(nanlocy,i)=NaN;
    wscmeans_highpassed(nanlocy,i)=NaN;
    sstlapmeanp_highpassed(nanlocy,i)=NaN;
end
lonp = lons;
latp = lats;



%%
load('IfremerdailyWindstressAndSSTAtlantic.mat');
[I,J,K]=size(taux);
taux(isnan(sst))=NaN;
tauy(isnan(sst))=NaN;

lona = lons;clear lons
lata = lats;clear lats

taux = permute(taux,[2,1,3]);
tauy = permute(tauy,[2,1,3]);
sst = permute(sst,[2,1,3]);

taux1=movmean(taux,4,1,'omitnan');
taux=movmean(taux1,4,2,'omitnan');
tauy1=movmean(tauy,4,1,'omitnan');
tauy=movmean(tauy1,4,2,'omitnan');
sst1=movmean(sst,4,1,'omitnan');
sst=movmean(sst1,4,2,'omitnan');
clear taux1 tauy1 sst1
tauxmean=nanmean(taux,3);
tauymean=nanmean(tauy,3);

% partial derivatives of taux and tauy
[txx,txy]=zh_grad2(tauxmean,lona,lata);
clear txx 
[tyx,tyy]=zh_grad2(tauymean,lona,lata);
clear tyy 
% curl of tau
wsc=tyx-txy;
clear txy tyx

sst=sst-273.15;
sstmean=nanmean(sst,3);
wscmeana=nanmean(wsc,3);
% gradient of sst
[sstx,ssty]=zh_grad2(sst,lona,lata);

[sstxx,sstxy]=zh_grad2(sstx,lona,lata);
clear sstxy;
[sstyx,sstyy]=zh_grad2(ssty,lona,lata);
clear sstyx;
% SST laplacian
sstlap=sstxx+sstyy;clear sstxx sstyy
sstlapmean_a=nanmean(sstlap,3);clear sstlap

sstlapmean_lowpassed1=movmean(sstlapmean_a,16,1,'omitnan');
sstlapmean_lowpassed2=movmean(sstlapmean_lowpassed1,16,2,'omitnan');clear sstlapmean_lowpassed1
sstlapmean_lowpassed3=movmean(sstlapmean_a-sstlapmean_lowpassed2,16,1,'omitnan');
sstlapmean_lowpassed4=movmean(sstlapmean_lowpassed3,16,2,'omitnan');clear sstlapmean_lowpassed3
sstlapmeana_highpassed=sstlapmean_a-(sstlapmean_lowpassed4+sstlapmean_lowpassed2);clear sstlapmean_lowpassed2 sstlapmean_lowpassed4

taudir=atan2(tauy,taux);
taumag=sqrt(taux.^2+tauy.^2);clear taux tauy
% sst gradient direction and magnitude
sstgraddir=atan2(ssty,sstx);
sstgradmag=sqrt(sstx.^2+ssty.^2);clear sstx ssty
%crosswind sst gradient
crosswindsstgrad = sstgradmag.*sin(taudir-sstgraddir);clear taudir sstgraddir sstgradmag
cwsga_mean=nanmean(crosswindsstgrad,3);clear crosswindsstgrad

wscmeana_lowpassed1=movmean(wscmeana,16,1,'omitnan');
wscmeana_lowpassed2=movmean(wscmeana_lowpassed1,16,2,'omitnan');clear wscmeann_lowpassed1
wscmeana_lowpassed3=movmean(wscmeana-wscmeana_lowpassed2,16,1,'omitnan');
wscmeana_lowpassed4=movmean(wscmeana_lowpassed3,16,2,'omitnan');clear wscmeann_lowpassed3
wscmeana_highpassed=wscmeana-(wscmeana_lowpassed4+wscmeana_lowpassed2);clear wscmeana_lowpassed2 wscmeana_lowpassed4

cwsgamean_lowpassed1=movmean(cwsga_mean,16,1,'omitnan');
cwsgamean_lowpassed2=movmean(cwsgamean_lowpassed1,16,2,'omitnan');clear cwsgamean_lowpassed1
cwsgamean_lowpassed3=movmean(cwsga_mean-cwsgamean_lowpassed2,16,1,'omitnan');
cwsgamean_lowpassed4=movmean(cwsgamean_lowpassed3,16,2,'omitnan');clear cwsgamean_lowpassed3
cwsgamean_highpassed=cwsga_mean-(cwsgamean_lowpassed4+cwsgamean_lowpassed2);clear cwsgamean_lowpassed2 cwsgamean_lowpassed4

wscmeana_1d=reshape(wscmeana_highpassed,[],1);
cwsgamean_1d=reshape(cwsgamean_highpassed,[],1);
sstlapamean_1d=reshape(sstlapmeana_highpassed,[],1);
%plot crosswind sst gradient
ss(1)=subplot(1,2,1);
contourf(lona,lata,wscmeana_highpassed,-1e-7:1e-9:1e-7,'edgecolor','None');
hold on;
xlabel('Longitude');
ylabel('Latiutude');
cm = redblue(101);
colormap(cm);
colorbar();
worldmap3Atl(2);
xlim([-60 10]);
ylim([-30 30]);
title('wscmeana highpassed');

ss(2)=subplot(1,2,2);
contourf(lona,lata,cwsgamean_highpassed,-1e-5:1e-7:1e-5,'LineStyle','None');
hold on;
ylabel('Latitude');
xlabel('Longitude');
cm = redblue(101);
colormap(cm);
colorbar();
worldmap3Atl(2);
xlim([-60 10]);
ylim([-30 30]);
title('Mean Crosswind SST gradient highpassed');

%%
marg_h =[0.06 0.04]; 
marg_w =[0.035 0.02];
Nh = 3;
Nw = 4;
gap=[0.055 0.035 0.045];

axh = (1-sum(marg_h)-1*gap(1))/1.5; 
axw = (1-sum(marg_w)-2*gap(2)-gap(3))/4;

axh1 = (1-sum(marg_h)-1*gap(1))/Nh; 
axw1 = (1-sum(marg_w)-2*gap(2)-gap(3))/4;

px = [marg_w(1) marg_w(1)+axw+gap(2)  marg_w(1)+axw+gap(2)+gap(3)+axw1 marg_w(1)+2*axw+2*gap(2)+gap(3)+axw1];
py = [1-marg_h(2)-axh 1-marg_h(2)-axh-gap(1)-axh1]; 

set(gcf,'units','centimeters','position',[50,-5,38.5,25.5]);
%plotting part of 4 subplots
ss(1)=subplot('Position',[px(1) py(1) axw axh]);
contourf(lonp,latp,wscmeans_highpassed,-5e-7:0.5e-9:5e-7,'LineStyle','None');
hold on;
%contour(lonp,latp,sstlapmeanp_highpassed,[-1e-11,0] ,'LineColor','black','LineStyle','--');
contour(lonp,latp,sstlapmeanp_highpassed,0:100:100,'LineColor','black','LineWidth',1.2);
%contour(lonp,latp,sstlapmeanp_highpassed,0:1e-11:1e-11,'LineColor','black');
cm=redblue(101);
colorbar();
colormap(gca,cm);
worldmap3(2);
xlim([210 290]);
ylim([-30 30]);
caxis([-3e-8,3e-8])
%xlabel('Longitude','FontSize',11);  
%ylabel('Latitude','FontSize',11);
xticks([210 240 270]);
xticklabels({'150°„W','120°„W','90°„W'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
ax = gca;
set(ax, 'FontSize', 15); 
ax.YAxis.TickLabelGapOffset = -2.5;
title('a)','Units','normalized','Position',[0, 1.005],'HorizontalAlignment','left','FontSize',15);
set(ax,'Color',0.6*[1 1 1],'TickDir','out');
ax.XAxis.TickLabelGapOffset = -4;

ss(2)=subplot('Position',[px(1) py(2) axw1 axh1]);
histogram2(wscmeans_1d.*1e8,sstlapmeanp_1d.*1e11,1500,'Normalization','probability','Displaystyle','tile');
hold on;
clb = colorbar();
clb.Ticks = [0.005,0.01,0.015,0.02];
clb.TickLabels={'0.5%', '1%','1.5%','2%'};
cm1=whitebluered(40);
colormap(gca,cm1);
ylim([-4,4]);
xlim([-3,3]);
grid off;
ax = gca;
ax.XAxis.TickLabelGapOffset = -2.5;
ax.XLabel.Position(2) = ax.XLabel.Position(2)-0.18; 
ax.YLabel.Position(1) = ax.YLabel.Position(1); 
ax.YLabel.Position(2) = ax.YLabel.Position(2)-0.04; 
set(ax, 'FontSize', 15,'TickDir','out');
ylabel('$\nabla^{2} T$ (10$^{-10}$ $^\circ\mathrm{C}$ m$^{-2}$)','Interpreter','latex');
xlabel('$\nabla \times \vec{\tau}$ (10$^{-8}$Nm$^{-3}$)','Interpreter','latex');
title('e)','Units','normalized','Position',[0, 1.005],'HorizontalAlignment','left','FontSize',15);

ss(3)=subplot('Position',[px(3) py(1) axw axh]);
contourf(lonp,latp,wscmeans_highpassed,-5e-7:1e-9:5e-7,'LineStyle','None');
hold on;
contour(lonp,latp,cwsgsmean_highpassed,[-2e-7,0] ,'LineColor','black','LineStyle','--');
contour(lonp,latp,cwsgsmean_highpassed,0:100:100,'LineColor','black','LineWidth',1.2);
contour(lonp,latp,cwsgsmean_highpassed,0:2e-7:2e-7,'LineColor','black');
cm=redblue(101);
colorbar();
colormap(gca,cm);
worldmap3(2);
xlim([210 290]);
ylim([-30 30]);
caxis([-3e-8,3e-8])
%xlabel('Longitude','FontSize',11);  
%ylabel('Latitude','FontSize',11);
xticks([210 240 270]);
xticklabels({'150°„W','120°„W','90°„W'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
title('c)','Units','normalized','Position',[0, 1.005],'HorizontalAlignment','left','FontSize',11);
set(gca,'Color',0.6*[1 1 1],'TickDir','out');
ax = gca;
ax.YAxis.TickLabelGapOffset = -2.5;
ax.XAxis.TickLabelGapOffset = -4;
set(ax, 'FontSize', 15); 

ss(4)=subplot('Position',[px(3) py(2) axw1 axh1]);
histogram2(wscmeans_1d.*1e8,cwsgsmean_1d.*1e6,800,'Normalization','probability','Displaystyle','tile');
hold on;
clb = colorbar();
clb.Ticks = [0.005,0.01,0.015,0.02];
clb.TickLabels={'0.5%', '1%','1.5%','2%'};
cm1=whitebluered(40);
colormap(gca,cm1);
ylabel(' crosswind $\nabla T$ ($10^{-6}$ $^\circ\mathrm{C}$ m$^{-1}$)','Interpreter','latex');
xlabel('$\nabla \times \vec{\tau}$ (10$^{-8}$Nm$^{-3}$)','Interpreter','latex');
ylim([-1,1]);
xlim([-3,3]);
grid off;
ax = gca;
ax.XAxis.TickLabelGapOffset = -2.5;
ax.XLabel.Position(2) = ax.XLabel.Position(2)-0.04; 
ax.YLabel.Position(1) = ax.YLabel.Position(1); 
ax.YLabel.Position(2) = ax.YLabel.Position(2)-0.03; 
set(ax, 'FontSize', 15,'TickDir','out'); 
title('g)','Units','normalized','Position',[0, 1.005],'HorizontalAlignment','left');

ss(5)=subplot('Position',[px(2) py(1) axw axh]);
contourf(lona,lata,wscmeana_highpassed,-5e-7:1e-9:5e-7,'LineStyle','None');
hold on;
%contour(lona,lata,cwsgamean_highpassed,[-2e-7,0] ,'LineColor','black','LineStyle','--');
contour(lona,lata,sstlapmeana_highpassed,0:100:100,'LineColor','black','LineWidth',1.2);
%contour(lona,lata,cwsgamean_highpassed,0:2e-7:2e-7,'LineColor','black');
cm=redblue(101);
colorbar();
colormap(gca,cm);
worldmap3Atl(2);
xlim([-60 10]);
ylim([-30 30]);
caxis([-3e-8,3e-8])
%xlabel('Longitude');
%ylabel('Latitude');
xticks([-60 -30 0]);
xticklabels({'60°„W','30°„W','0°„'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
title('b)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
ax = gca;
ax.YAxis.TickLabelGapOffset = -2.5;
ax.XAxis.TickLabelGapOffset = -4;
set(ax, 'FontSize', 15); 
set(gca,'Color',0.6*[1 1 1],'TickDir','out');

ss(6)=subplot('Position',[px(2) py(2) axw1 axh1]);
histogram2(wscmeana_1d.*1e8,sstlapamean_1d.*1e11,800,'Normalization','probability','Displaystyle','tile');
hold on;
clb = colorbar();
clb.Ticks = [0.005,0.01,0.015,0.02];
clb.TickLabels={'0.5%', '1%','1.5%','2%'};
cm1=whitebluered(40);
colormap(gca,cm1);
%ylabel('SST laplacian (10^{-10} ^oC m^{-2})','FontSize',11);
xlabel('$\nabla \times \vec{\tau}$ (10$^{-8}$Nm$^{-3}$)','Interpreter','latex');
ylim([-4,4]);
xlim([-3,3]);
grid off;
title('f)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
ax = gca;
ax.XAxis.TickLabelGapOffset = -2.5;
ax.XLabel.Position(2) = ax.XLabel.Position(2)-0.18; 
set(ax, 'FontSize', 15,'TickDir','out'); 

ss(7)=subplot('Position',[px(4) py(1) axw axh]);
contourf(lona,lata,wscmeana_highpassed,-5e-7:1e-9:5e-7,'LineStyle','None');
hold on;
contour(lona,lata,cwsgamean_highpassed,[-2e-7,0] ,'LineColor','black','LineStyle','--');
contour(lona,lata,cwsgamean_highpassed,0:100:100,'LineColor','black','LineWidth',1.2);
contour(lona,lata,cwsgamean_highpassed,0:2e-7:2e-7,'LineColor','black');
cm=redblue(101);
colorbar();
colormap(gca,cm);
worldmap3Atl(2);
xlim([-60 10]);
ylim([-30 30]);
caxis([-3e-8,3e-8])
%xlabel('Longitude');
%ylabel('Latitude');
xticks([-60 -30 0]);
xticklabels({'60°„W','30°„W','0°„'});
yticks([flip(-1*[0 10 20 30]) [10 20 30]]);
yticklabels({'30°„S','20°„S','10°„S','0°„','10°„N','20°„N','30°„N'});
title('d)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
ax = gca;
ax.YAxis.TickLabelGapOffset = -2.5;
ax.XAxis.TickLabelGapOffset = -4;
set(ax, 'FontSize', 15); 
set(gca,'Color',0.6*[1 1 1],'TickDir','out');

ss(8)=subplot('Position',[px(4) py(2) axw1 axh1]);
histogram2(wscmeana_1d.*1e8,cwsgamean_1d.*1e6,500,'Normalization','probability','Displaystyle','tile');
hold on;
clb = colorbar();
clb.Ticks = [0.005,0.01,0.015];
clb.TickLabels={'0.5%', '1%','1.5%'};
cm1=whitebluered(40);
colormap(gca,cm1);
%ylabel('crosswind SST gradient (10^{-6} ^oCm^{-1})','FontSize',11);
xlabel('$\nabla \times \vec{\tau}$ (10$^{-8}$Nm$^{-3}$)','Interpreter','latex');
ylim([-1,1]);
xlim([-3,3]);
grid off;
title('h)','Units','normalized','Position',[0, 1],'HorizontalAlignment','left');
ax = gca;
ax.XAxis.TickLabelGapOffset = -2.5;
ax.XLabel.Position(2) = ax.XLabel.Position(2)-0.04; 
set(ax, 'FontSize', 15,'TickDir','out'); 