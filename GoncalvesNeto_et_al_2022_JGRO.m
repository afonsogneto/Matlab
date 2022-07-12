%% GoncalvesNeto_et_al_2022_JGRO

% AGN 20211102

% This script provides the code to plot the figures in Goncalves Neto et al
% (2022) submitted to JGR-Oceans

%% Figure 1
% Schematic of the circulation in the Northwest Atlantic Ocean and the water 
% column temperature change in the period 2009-2018 relative to 2001-2007. The 
% filled red (blue) contours indicate warming (cooling) of the vertically-
% averaged ocean temperature from the EN4 objective analysis (Good et al., 2013)
% to 2,000 m or the seafloor if it is shallower than 2,000 m in 0.5oC increments
% (change in the unshaded region is <0.25oC). Background in grayscale shows the
% bathymetry of the region, with darker shades representing shallower areas. The
% main surface circulation patterns associated with the Gulf Stream (red) and 
% the Labrador Current (blue) systems are identified with arrows. Coastal and 
% shelf areas of interest are indicated. TGB stands for Tail of the Grand Banks 
% and SENR, for Southeastern Newfoundland Ridge.

% Set  directory to save data
data_dir = '/Users/afonso/Documents/Research/TGB/data/';
% Set directory to save figures and tables
output_dir = '/Users/afonso/Documents/Research/TGB/output/';

% Bathymetry (first layer)

load([output_dir,'etopo1_TGB.mat'],'lat_z','lon_z','z')

figure('Position', [100, 200, 800, 600])
subplot('Position',[0.06 0.09 0.88 0.85])
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
hold on
pcolorm(lat_z,lon_z,z)
shading flat
colormap(cmocean('gray','negative'))
caxis([-4500 2000])
geoshow('GSHHS_i_L1.shp','FaceColor','k')%[.7 .7 .7])
plabel('PLabelLocation',5,'fontsize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south','fontsize',14)
framem
tightmap
set(gcf, 'Color', 'w')
saveas(gcf,[output_dir,'GoncalvesNeto_et_al_2020_Nature/GoncalvesNeto_et_al_2020_Nature_Figure1_layer1'],'epsc')
%export_fig(gcf,[output_dir,'GoncalvesNeto_et_al_2020_Nature/GoncalvesNeto_et_al_2020_Nature_Figure1_layer1'],'-eps','-transparent')

% Temperature change (second layer)

load([data_dir,'en4_NWA.mat'],'t_en4','s_en4','lon_en4','lat_en4','time_en4')

[~,in1] = min(abs(time_en4 - datenum('20010115','yyyymmdd')));
[~,out1] = min(abs(time_en4 - datenum('20071215','yyyymmdd')));
[~,in2] = min(abs(time_en4 - datenum('20090115','yyyymmdd')));
[~,out2] = min(abs(time_en4 - datenum('20181215','yyyymmdd')));

t_en4_2001_2007 = nanmean(nanmean(t_en4(:,:,1:30,in1:out1),4),3);
t_en4_2009_2018 = nanmean(nanmean(t_en4(:,:,1:30,in2:out2),4),3);

figure('Position', [100, 200, 800, 600])
subplot('Position',[0.06 0.09 0.88 0.85])
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
hold on
t_diff = t_en4_2009_2018' - t_en4_2001_2007';
pcolorm(lat_en4,lon_en4,t_diff)
shading flat
fillm([53,53,55,55],[-67,-63,-63,-67],'k')
%rectangle('Position',[0.2,0.8,0.1,0.1],'FaceColor','k','EdgeColor','k')
%rectangle('Position',[52,-68,2,4],'FaceColor','k','EdgeColor','k')
colormap([5/256,113/256,175/256;...
    5/256,113/256,176/256;...
    247/256,247/256,247/256;...
    247/256,247/256,247/256;...
    244/256,165/256,130/256;...
    244/256,165/256,130/256;...
    214/256,96/256,77/256;...
    214/256,96/256,77/256;...
    178/256,24/256,43/256;...
    178/256,24/256,43/256])

set(gca,'FontSize',14)
geoshow('GSHHS_i_L1.shp','FaceColor','k')%[.7 .7 .7])
plabel('PLabelLocation',5,'Fontsize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',14)
framem
tightmap
%cbarrow

c = colorbar('Orientation','Horizontal','Position',[0.2 0.7 0.27 0.04],'Color',[1-eps,1,1],'Fontsize',10);%colorbar over land
c.Label.String = 'Temperature Difference (^\circC)';
c.Ticks = [-0.75 -0.25 0.25 0.75 1.25 1.75];
set(c,'Fontsize',14)
caxis([-0.75 1.75])

set(gcf, 'Color', 'w')
saveas(gcf,[output_dir,'GoncalvesNeto_et_al_2020_Nature/GoncalvesNeto_et_al_2020_Nature_Figure1_layer2'],'epsc')

clearvars -except *_dir
close all

%% Figure 2

% Maps of SSH from satellite and model

% Load HYCOM ATL
cd /Users/afonso/Documents/Research/TGB2/Data/

dir_output = '/Users/afonso/Documents/Research/TGB2/output/';
%load([dir_output,'etopo1_TGB.mat'],'lat_z','lon_z','z')
load('map_GSS_ATLg08.mat')


% Load HYCOM ATL
dir_input = 'Hycom_ATLg08/';
cd(dir_input)
dir_list = dir('*.mat');

ssh_hycom_atl = NaN(size(x,1),size(x,2),length(dir_list));
for ii = 1:length(dir_list)
    a = load(dir_list(ii).name);
    ssh_hycom_atl(:,:,ii) = a.ssh;
    clear a
end
clear ii

% Crop data to <47N to avoid irregular grid
ssh_hycom_atl = ssh_hycom_atl(1:215,:,:);
lat_hycom_atl = y(1:215,:);
lon_hycom_atl = x(1:215,:);
pdep = pdep(1:215,:);
clear x y

% Load Altimetry
cd(dir_output)
%cd /Users/afonso/Documents/Research/TGB/output

load('SSH_TGB','time_month','lat_ssh','lon_ssh','adt')
time_month = time_month(1:end-12);

adt = adt(lon_ssh >= 360-77 & lon_ssh <= 360-40,lat_ssh >= 34 & lat_ssh <= 47,1:300);
lat_ssh = lat_ssh(lat_ssh >= 34 & lat_ssh <= 47);
lon_ssh = lon_ssh(lon_ssh >= 360-77 & lon_ssh <= 360-40);
lon_ssh = lon_ssh - 360;
[lon_ssh_mesh, lat_ssh_mesh] = meshgrid(lon_ssh,lat_ssh);

% Interp2 hycom to altimetry grid
ssh_hycom_atl1 = NaN(length(lat_ssh),length(lon_ssh),length(dir_list));
for ii = 1:length(dir_list)
    ssh_hycom_atl1(:,:,ii) = interp2(lon_hycom_atl,lat_hycom_atl,ssh_hycom_atl(:,:,ii),lon_ssh_mesh,lat_ssh_mesh);
end
ssh_hycom_atl1 = permute(ssh_hycom_atl1,[2,1,3]);

clear ii
clear ssh_hycom_atl lat_hycom_atl lon_hycom_atl pdep
clear dir_list

% Load bathymetry from ETOPO1
x = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','x');
lon_z = x(x<=-40 & x>=-77);
lon_z = lon_z(1:5:end);
y = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','y');
lat_z = y(y<=47 & y>=34);
lat_z = lat_z(1:5:end);
z_aux = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','z');
z = z_aux(x<=-40 & x>=-77,y<=47 & y>=34); z = z';
z = z(1:5:end,1:5:end);
clear x y z_aux

ssh_tgb_alt = squeeze(nanmean(nanmean(adt(lon_ssh >= -53 & lon_ssh <= -50,...
    lat_ssh >= 41 & lat_ssh <=43,:),2),1));
ssh_tgb_hycom = squeeze(nanmean(nanmean(ssh_hycom_atl1(lon_ssh >= -53 & lon_ssh <= -50,...
    lat_ssh >= 41 & lat_ssh <=43,:),2),1));
ssh_tgb_alt2 = squeeze(nanmean(nanmean(adt(lon_ssh >= -51 & lon_ssh <= -46,...
    lat_ssh >= 39 & lat_ssh <=43,:),2),1));
ssh_tgb_hycom2 = squeeze(nanmean(nanmean(ssh_hycom_atl1(lon_ssh >= -51 & lon_ssh <= -46,...
    lat_ssh >= 39 & lat_ssh <=43,:),2),1));

corrcoef(ssh_tgb_hycom,ssh_tgb_alt)
ipt_hycom = findchangepts(ssh_tgb_hycom)
datestr(time_month(ipt_hycom))
ipt_alt = findchangepts(ssh_tgb_alt)
datestr(time_month(ipt_alt))

corrcoef(ssh_tgb_hycom,ssh_tgb_hycom2)
corrcoef(ssh_tgb_alt,ssh_tgb_alt2)
corrcoef(ssh_tgb_hycom2,ssh_tgb_alt2)
ipt_hycom = findchangepts(ssh_tgb_hycom2)
ipt_hycom = findchangepts(ssh_tgb_alt2)

ipt_hycom = NaN(length(lon_ssh),length(lat_ssh));
ipt_alt = NaN(length(lon_ssh),length(lat_ssh));
for ii = 1:length(lon_ssh)
    for jj = 1:length(lat_ssh)
        if sum(~isnan(squeeze(ssh_hycom_atl1(ii,jj,:))))>0
        ipt_hycom(ii,jj) = findchangepts(squeeze(ssh_hycom_atl1(ii,jj,:)));
        ipt_alt(ii,jj) = findchangepts(squeeze(adt(ii,jj,:)));
        end
    end
end

% Spatial correlation

% reshape
ssh_alt_aux = reshape(squeeze(nanmean(adt(:,:,1:end)*100,3)),[],1);
ssh_hycom_aux = reshape(squeeze(nanmean(ssh_hycom_atl1(:,:,1:end),3)),[],1);
% remove NaNs
ssh_alt = ssh_alt_aux;
ssh_hycom = ssh_hycom_aux;
ssh_alt(isnan(ssh_alt_aux) | isnan(ssh_hycom_aux)) = [];
ssh_hycom(isnan(ssh_alt_aux) | isnan(ssh_hycom_aux)) = [];
% calculate coefficients
[ssh_corr,ssh_p] = corrcoef(ssh_alt,ssh_hycom);
ssh_poly = polyfit(ssh_alt,ssh_hycom,1);
ssh_hycom2 = polyval(ssh_poly,ssh_alt);
%plot
figure
hold on
scatter(ssh_alt,ssh_hycom,10,rgb('gray'),'filled')
plot(ssh_alt,ssh_hycom2,'linewidth',2,'Color',rgb('Black'))
ylabel('HYCOM','FontWeight','Bold','FontSize',14)
xlabel('Altimetry','FontWeight','Bold','FontSize',14)
title(['SSH Spatial Correlation (rho = ',num2str(round(ssh_corr(1,2),2)),', p < 0.001)'],...
    'FontSize',18)
saveas(gcf,[dir_output,'manuscript_v3/TGB2_Figure2_v3_corr'],'png')


colormap_aux = readtable([dir_output,'5w_BRgpb.csv'], 'HeaderLines',1);
colormap_wave = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'4w_ROTB.csv'], 'HeaderLines',1);
colormap_wave2 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'3-wave-yellow-grey-blue.csv'], 'HeaderLines',1);
colormap_wave3 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'5-step-melow-wave.csv'], 'HeaderLines',1);
colormap_wave4 = table2array(colormap_aux(:,3:5));
clear colormap_aux

figure('Position', [100, 200, 900, 700])
subplot('Position',[0.05 0.53 0.93 0.46])
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-77 -40])
hold on
pcolorm(lat_ssh,lon_ssh,squeeze(nanmean(adt(:,:,1:end)*100,3))')
shading interp
%colormap(crameri('-lapaz'))
colormap(flipud(colormap_wave3(3:43,:)))
c = colorbar('Orientation','Horizontal','Position',[0.08 0.88 0.16 0.04],...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'SSH (cm)';
set(c,'Fontsize',15)%,'Fontweight','bold')
caxis([-80 80])
% contourm(lat_ssh,lon_ssh,squeeze(nanmean(adt(:,:,ipt:end)*100,3))' - ...
%     squeeze(nanmean(adt(:,:,1:ipt)*100,3))' - ...
%     nanmean(nanmean(squeeze(nanmean(adt(:,:,ipt:end)*100,3))' - ...
%     squeeze(nanmean(adt(:,:,1:ipt)*100,3))')),[5 5],'LineStyle','-','Color',rgb('Black'),'LineWidth',1.5)
% contourm(lat_ssh,lon_ssh,squeeze(nanmean(adt(:,:,ipt:end)*100,3))' - ...
%     squeeze(nanmean(adt(:,:,1:ipt)*100,3))' - ...
%     nanmean(nanmean(squeeze(nanmean(adt(:,:,ipt:end)*100,3))' - ...
%     squeeze(nanmean(adt(:,:,1:ipt)*100,3))')),[-5 -5],'LineStyle','--','Color',rgb('Black'),'LineWidth',1.5)
contourm(lat_z,lon_z,z,[-1000 -1000],'Color',rgb('DarkSlateGray'),'LineWidth',1)
contourm(lat_z,lon_z,z,[-4000 -4000],'Color',rgb('DarkSlateGray'),'LineWidth',1)
textm(42,-76.5,'(a)','Color',rgb('Snow'),'Fontsize',20,'Fontweight','bold')
textm(46,-57,'Altimetry','Color',rgb('Black'),'Fontsize',20)%,'Fontweight','bold')
f1 = fillm([41,43,43,41],[-53,-53,-50,-50],'w');
f1.FaceColor = 'w';
f1.EdgeColor = 'k';
f1.FaceAlpha = 0;
f1.EdgeAlpha = 1;
f1.LineWidth = 2.5;
%textm(43.5,-74.8,['Layer ',nlayer],'Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%title('HYCOM Atl Layer 20 - Mean Speed (2003-2012)')
plabel('PLabelLocation',5,'Fontsize',15)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

%figure('Position', [100, 200, 900, 367])
subplot('Position',[0.05 0.06 0.93 0.46])
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-77 -40])
hold on
pcolorm(lat_ssh,lon_ssh,squeeze(nanmean(ssh_hycom_atl1(:,:,1:end),3))')
shading interp
%colormap(crameri('-lapaz'))
colormap(flipud(colormap_wave3(3:43,:)))
c = colorbar('Orientation','Horizontal','Position',[0.08 0.41 0.16 0.04],...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'SSH (cm)';
set(c,'Fontsize',15)%,'Fontweight','bold')
caxis([-80 80])
% contourm(lat_ssh,lon_ssh,squeeze(nanmean(ssh_hycom_atl1(:,:,ipt:end)*100,3)) - ...
%     squeeze(nanmean(ssh_hycom_atl1(:,:,1:ipt)*100,3)) - ...
%     nanmean(nanmean(squeeze(nanmean(ssh_hycom_atl1(:,:,ipt:end)*100,3)) - ...
%     squeeze(nanmean(ssh_hycom_atl1(:,:,1:ipt)*100,3)))),[5 5],'LineStyle','-','Color',rgb('Black'),'LineWidth',1.5)
% contourm(lat_ssh,lon_ssh,squeeze(nanmean(ssh_hycom_atl1(:,:,ipt:end)*100,3)) - ...
%     squeeze(nanmean(ssh_hycom_atl1(:,:,1:ipt)*100,3)) - ...
%     nanmean(nanmean(squeeze(nanmean(ssh_hycom_atl1(:,:,ipt:end)*100,3)) - ...
%     squeeze(nanmean(ssh_hycom_atl1(:,:,1:ipt)*100,3)))),[-5 -5],'LineStyle','--','Color',rgb('Black'),'LineWidth',1.5)
contourm(lat_z,lon_z,z,[-1000 -1000],'Color',rgb('DarkSlateGray'),'LineWidth',1)
contourm(lat_z,lon_z,z,[-4000 -4000],'Color',rgb('DarkSlateGray'),'LineWidth',1)
textm(42,-76.5,'(b)','Color',rgb('Snow'),'Fontsize',20,'Fontweight','bold')
textm(46,-57,'HYCOM','Color',rgb('Black'),'Fontsize',20)%,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%title('HYCOM Atl Layer 20 - Mean Speed (2003-2012)')
plabel('PLabelLocation',5,'Fontsize',15)
mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

saveas(gcf,[dir_output,'manuscript_v3/TGB2_Figure2_v3'],'png')

figure
hold on
plot(time_month,ssh_tgb_hycom)
plot(time_month,ssh_tgb_alt)


figure('Position', [100, 200, 900, 700])
subplot('Position',[0.05 0.53 0.93 0.46])
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-77 -40])
hold on
pcolorm(lat_ssh,lon_ssh,ipt_alt')
shading interp
%colormap(crameri('-lapaz'))
%colormap(flipud(colormap_wave3(3:43,:)))
c = colorbar('Orientation','Horizontal','Position',[0.08 0.88 0.16 0.04],...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'SSH (cm)';
set(c,'Fontsize',15)%,'Fontweight','bold')
caxis([50 250])
contourm(lat_ssh,lon_ssh,ipt_alt',[181 181],'Color',rgb('DarkSlateGray'),'LineWidth',2)
contourm(lat_ssh,lon_ssh,ipt_alt',[193 193],'Color',rgb('DarkSlateGray'),'LineWidth',2)
% contourm(lat_ssh,lon_ssh,squeeze(nanmean(adt(:,:,ipt:end)*100,3))' - ...
%     squeeze(nanmean(adt(:,:,1:ipt)*100,3))' - ...
%     nanmean(nanmean(squeeze(nanmean(adt(:,:,ipt:end)*100,3))' - ...
%     squeeze(nanmean(adt(:,:,1:ipt)*100,3))')),[5 5],'LineStyle','-','Color',rgb('Black'),'LineWidth',1.5)
% contourm(lat_ssh,lon_ssh,squeeze(nanmean(adt(:,:,ipt:end)*100,3))' - ...
%     squeeze(nanmean(adt(:,:,1:ipt)*100,3))' - ...
%     nanmean(nanmean(squeeze(nanmean(adt(:,:,ipt:end)*100,3))' - ...
%     squeeze(nanmean(adt(:,:,1:ipt)*100,3))')),[-5 -5],'LineStyle','--','Color',rgb('Black'),'LineWidth',1.5)
contourm(lat_z,lon_z,z,[-1000 -1000],'Color',rgb('DarkSlateGray'),'LineWidth',1)
contourm(lat_z,lon_z,z,[-4000 -4000],'Color',rgb('DarkSlateGray'),'LineWidth',1)
textm(42,-76.5,'(a)','Color',rgb('Snow'),'Fontsize',20,'Fontweight','bold')
textm(46,-57,'Altimetry','Color',rgb('Black'),'Fontsize',20)%,'Fontweight','bold')
f1 = fillm([41,43,43,41],[-53,-53,-50,-50],'w');
f1.FaceColor = 'w';
f1.EdgeColor = 'k';
f1.FaceAlpha = 0;
f1.EdgeAlpha = 1;
f1.LineWidth = 2.5;
%textm(43.5,-74.8,['Layer ',nlayer],'Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%title('HYCOM Atl Layer 20 - Mean Speed (2003-2012)')
plabel('PLabelLocation',5,'Fontsize',15)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

%figure('Position', [100, 200, 900, 367])
subplot('Position',[0.05 0.06 0.93 0.46])
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-77 -40])
hold on
pcolorm(lat_ssh,lon_ssh,ipt_hycom')
shading interp
%colormap(crameri('-lapaz'))
%colormap(flipud(colormap_wave3(3:43,:)))
c = colorbar('Orientation','Horizontal','Position',[0.08 0.88 0.16 0.04],...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'SSH (cm)';
set(c,'Fontsize',15)%,'Fontweight','bold')
caxis([50 250])
contourm(lat_ssh,lon_ssh,ipt_hycom',[181 181],'Color',rgb('DarkSlateGray'),'LineWidth',2)
contourm(lat_ssh,lon_ssh,ipt_hycom',[193 193],'Color',rgb('DarkSlateGray'),'LineWidth',2)
% contourm(lat_ssh,lon_ssh,squeeze(nanmean(ssh_hycom_atl1(:,:,ipt:end)*100,3)) - ...
%     squeeze(nanmean(ssh_hycom_atl1(:,:,1:ipt)*100,3)) - ...
%     nanmean(nanmean(squeeze(nanmean(ssh_hycom_atl1(:,:,ipt:end)*100,3)) - ...
%     squeeze(nanmean(ssh_hycom_atl1(:,:,1:ipt)*100,3)))),[5 5],'LineStyle','-','Color',rgb('Black'),'LineWidth',1.5)
% contourm(lat_ssh,lon_ssh,squeeze(nanmean(ssh_hycom_atl1(:,:,ipt:end)*100,3)) - ...
%     squeeze(nanmean(ssh_hycom_atl1(:,:,1:ipt)*100,3)) - ...
%     nanmean(nanmean(squeeze(nanmean(ssh_hycom_atl1(:,:,ipt:end)*100,3)) - ...
%     squeeze(nanmean(ssh_hycom_atl1(:,:,1:ipt)*100,3)))),[-5 -5],'LineStyle','--','Color',rgb('Black'),'LineWidth',1.5)
contourm(lat_z,lon_z,z,[-1000 -1000],'Color',rgb('DarkSlateGray'),'LineWidth',1)
contourm(lat_z,lon_z,z,[-4000 -4000],'Color',rgb('DarkSlateGray'),'LineWidth',1)
textm(42,-76.5,'(b)','Color',rgb('Snow'),'Fontsize',20,'Fontweight','bold')
textm(46,-57,'HYCOM','Color',rgb('Black'),'Fontsize',20)%,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%title('HYCOM Atl Layer 20 - Mean Speed (2003-2012)')
plabel('PLabelLocation',5,'Fontsize',15)
mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

% contourm(lat,lon,squeeze(nanmean(thick,1))',[60 60],...
%     '-.k','LineWidth',0.5)%,'LevelStep',10,'ShowText','on')

%% Figure 3

% Cross-section of deployment location: mean sigma2, T, S and u-across.

dir_output = '/Users/afonso/Documents/Research/TGB2/output/';

%ncdisp([dir_output,'ATLg008_eulerian_LC.nc'])

% lat and lon
lat = ncread([dir_output,'ATLg008_eulerian_LC.nc'],'latitude');
lon = ncread([dir_output,'ATLg008_eulerian_LC.nc'],'longitude');
lat = lat';
lon = lon' - 360;

% Define lat and lon at LC
lat_lc = linspace(45.5,45,7);
lon_lc = linspace(-48.5,-48,7);

% length of section
gsw_distance([-48.5,-48],[45.5,45],0) / 1000

% Create variable for time
year = 2003:2012;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(120,1);
for yy = 1:10
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end


% thickness
thick = ncread([dir_output,'ATLg008_eulerian_LC.nc'],'thickness');
thick = thick/9806; % convert from pressure to metric

% density
density = ncread([dir_output,'ATLg008_eulerian_LC.nc'],'density');
density = density + 34; %offset

% temperature
temp = ncread([dir_output,'ATLg008_eulerian_LC.nc'],'temp');

% salinity
salt = ncread([dir_output,'ATLg008_eulerian_LC.nc'],'salt');

% u
u = ncread([dir_output,'ATLg008_eulerian_LC.nc'],'u');

% v
v = ncread([dir_output,'ATLg008_eulerian_LC.nc'],'v');

% ke
ke = ncread([dir_output,'ATLg008_eulerian_LC.nc'],'ke');


% calculate depth
depth = NaN(size(thick));
for dd = 1:size(thick,2)
    depth(:,dd,:,:) = nansum(thick(:,1:dd,:,:),2);
end
    
% interpolate variables to LC line
thick_lc = NaN(size(thick,1),size(thick,2),length(lat_lc));
depth_lc = NaN(size(thick,1),size(thick,2),length(lat_lc));
density_lc = NaN(size(thick,1),size(density,2),length(lat_lc));
temp_lc = NaN(size(thick,1),size(density,2),length(lat_lc));
salt_lc = NaN(size(thick,1),size(density,2),length(lat_lc));
u_lc = NaN(size(thick,1),size(density,2),length(lat_lc));
v_lc = NaN(size(thick,1),size(density,2),length(lat_lc));
ke_lc = NaN(size(thick,1),size(density,2),length(lat_lc));
for tt = 1:size(thick,1)
    for dd = 1:size(thick,2)
        thick_lc(tt,dd,:) = interp2(lon,lat,squeeze(thick(tt,dd,:,:))',...
            lon_lc,lat_lc);
        depth_lc(tt,dd,:) = interp2(lon,lat,squeeze(depth(tt,dd,:,:))',...
            lon_lc,lat_lc);
        density_lc(tt,dd,:) = interp2(lon,lat,squeeze(density(tt,dd,:,:))',...
            lon_lc,lat_lc);
        temp_lc(tt,dd,:) = interp2(lon,lat,squeeze(temp(tt,dd,:,:))',...
            lon_lc,lat_lc);
        salt_lc(tt,dd,:) = interp2(lon,lat,squeeze(salt(tt,dd,:,:))',...
            lon_lc,lat_lc);
        u_lc(tt,dd,:) = interp2(lon,lat,squeeze(u(tt,dd,:,:))',...
            lon_lc,lat_lc);
        v_lc(tt,dd,:) = interp2(lon,lat,squeeze(v(tt,dd,:,:))',...
            lon_lc,lat_lc);
        ke_lc(tt,dd,:) = interp2(lon,lat,squeeze(ke(tt,dd,:,:))',...
            lon_lc,lat_lc);
    end
end
clear tt dd

% Rotate velocity

uv_lc = sqrt(u_lc.^2 + v_lc.^2);

ucross_lc = NaN(size(u_lc));
ualong_lc = NaN(size(u_lc));
for tt = 1:size(u_lc,1)
    for dd = 1:size(u_lc,2)
        for ll = 1:size(u_lc,3)
            
            % direction of input vector (angle relative to north, increases clockwise)
            if u_lc(tt,dd,ll) >= 0 & v_lc(tt,dd,ll) >= 0
                theta_in = atand(u_lc(tt,dd,ll)/v_lc(tt,dd,ll));
                %theta_in_aux = theta_in
            elseif u_lc(tt,dd,ll) >= 0 & v_lc(tt,dd,ll) < 0
                theta_in = 180 + atand(u_lc(tt,dd,ll)/v_lc(tt,dd,ll));
                %theta_in_aux = theta_in - 90
            elseif u_lc(tt,dd,ll) < 0 & v_lc(tt,dd,ll) < 0
                theta_in = 180 + atand(u_lc(tt,dd,ll)/v_lc(tt,dd,ll));
                %theta_in_aux = theta_in - 180
            elseif u_lc(tt,dd,ll) < 0 & v_lc(tt,dd,ll) >= 0
                theta_in = 360 + atand(u_lc(tt,dd,ll)/v_lc(tt,dd,ll));
                %theta_in_aux = theta_in - 270
            end
            
            % direction of speed (angle relative to north, increases clockwise)
            % lat decreases and lon increases
            theta_rot = 180 + atand((lon_lc(end)-lon_lc(1))/(lat_lc(end)-lat_lc(1)));
            
            % direction of speed relative to the new coordinate system
            if theta_in > theta_rot
                theta_out = theta_in - theta_rot;
            else
                theta_out = 360 - theta_rot + theta_in;
            end
            
            ucross_lc(tt,dd,ll) = uv_lc(tt,dd,ll) * sind(theta_out);
            ualong_lc(tt,dd,ll) = uv_lc(tt,dd,ll) * cosd(theta_out);
            
            clear theta*
        end
    end
end
clear tt dd ll

% plot mean fields


% plot_location = [0.02 0.02 0.17 0.96;...
%     0.22 0.02 0.17 0.96;...
%     0.42 0.02 0.17 0.96;...
%     0.62 0.02 0.17 0.96;...
%     0.82 0.02 0.17 0.96];
% 
% cbar_location = [0.079 0.05 0.02 0.2;...
%     0.279 0.05 0.02 0.2;...
%     0.479 0.05 0.02 0.2;...
%     0.679 0.05 0.02 0.2;...
%     0.879 0.05 0.02 0.2];

plot_location = [0.03 0.1 0.21 0.88;...
    0.275 0.1 0.21 0.88;...
    0.52 0.1 0.21 0.88;...
    0.765 0.1 0.21 0.88];

cbar_location = [0.11 0.13 0.02 0.28;...
    0.355 0.13 0.02 0.28;...
    0.60 0.13 0.02 0.28;...
    0.845 0.13 0.02 0.28];

colormap_aux = readtable([dir_output,'5w_BRgpb.csv'], 'HeaderLines',1);
colormap_wave = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'4w_ROTB.csv'], 'HeaderLines',1);
colormap_wave2 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'3-wave-yellow-grey-blue.csv'], 'HeaderLines',1);
colormap_wave3 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'5-step-melow-wave.csv'], 'HeaderLines',1);
colormap_wave4 = table2array(colormap_aux(:,3:5));
clear colormap_aux

figure('Position', [0, 0, 1280, 740])
s1 = subplot('Position',plot_location(1,:));
hold on
rectangle('Position',[-48.5 0 0.5 3000],'FaceColor','k','EdgeColor','k')
pcolor(repmat(lon_lc,32,1),squeeze(nanmean(depth_lc,1))/1000,...
    squeeze(nanmean(density_lc,1)))
%pcolor(repmat(lon_lc,32,1),squeeze(depth_lc(1,:,:)),squeeze(density_lc(1,:,:)))
shading interp
set(gca,'ydir','reverse')
colormap(s1,flipud(colormap_wave(3:82,:)))
c = colorbar('Orientation','Vertical','Position',cbar_location(1,:),...
    'Color',[1-eps,1,1],'Fontsize',10); %colorbar over land
%c = colorbar('Location','manual','Orientation','Horizontal','Position',[0.17 0.7 0.25 0.03],'Color','w','Fontsize',10)%colorbar over land
%c = colorbar;
c.Label.String = 'Sigma 2 (kg m^{-3})';
%ylabel(h,'SST difference (^\circC)','Fontsize',14)
set(c,'Fontsize',30)
caxis([35 37])
ylim([0 3])
xlim([-48.5 -48])
set(gca,'ytick',[],'xtick',[-48.5,-48.25,-48],...
    'xticklabels',{'48.5','48.25','48'})
for dd = 1:size(depth_lc,2)
    line(lon_lc,squeeze(nanmean(depth_lc(:,dd,:),1))/1000,...
        'color',rgb('Gray'),'LineWidth',0.5)
end
line(lon_lc,squeeze(nanmean(depth_lc(:,18,:),1))/1000,'color','k','LineWidth',2)
line(lon_lc,squeeze(nanmean(depth_lc(:,19,:),1))/1000,'color','k','LineWidth',2)
%line(lon_lc,squeeze(nanmean(depth_lc(:,20,:),1))/1000,'color','k','LineWidth',2)
ylabel('Depth (km)')
xlabel('Longitude (°W)')
text(-48.47,1.5,'(a)','color','w','FontSize',30,'FontWeight','Bold')
set(gca,'FontSize',25)%,'FontWeight','Bold')

s2 = subplot('Position',plot_location(2,:));
hold on
rectangle('Position',[-48.5 0 0.5 3000],'FaceColor','k','EdgeColor','k')
pcolor(repmat(lon_lc,32,1),squeeze(nanmean(depth_lc,1))/1000,...
    squeeze(nanmean(temp_lc,1)))
%pcolor(repmat(lon_lc,32,1),squeeze(depth_lc(1,:,:)),squeeze(density_lc(1,:,:)))
shading interp
set(gca,'ydir','reverse')
colormap(s2,flipud(colormap_wave(3:82,:)))
c = colorbar('Orientation','Vertical','Position',cbar_location(2,:),...
    'Color',[1-eps,1,1],'Fontsize',10); %colorbar over land
%c = colorbar('Location','manual','Orientation','Horizontal','Position',[0.17 0.7 0.25 0.03],'Color','w','Fontsize',10)%colorbar over land
%c = colorbar;
c.Label.String = 'Temperature (^\circC)';
%ylabel(h,'SST difference (^\circC)','Fontsize',14)
set(c,'Fontsize',30)
caxis([2 10])
ylim([0 3])
xlim([-48.5 -48])
%set(gca,'xtick',[],'ytick',[])
set(gca,'xtick',[-48.5,-48.25,-48],...
    'xticklabels',{'48.5','48.25','48'})
for dd = 1:size(depth_lc,2)
    line(lon_lc,squeeze(nanmean(depth_lc(:,dd,:),1))/1000,...
        'color',rgb('Gray'),'LineWidth',0.5)
end
line(lon_lc,squeeze(nanmean(depth_lc(:,18,:),1))/1000,'color','k','LineWidth',2)
line(lon_lc,squeeze(nanmean(depth_lc(:,19,:),1))/1000,'color','k','LineWidth',2)
%line(lon_lc,squeeze(nanmean(depth_lc(:,20,:),1))/1000,'color','k','LineWidth',2)
%ylabel('Depth (m)')
xlabel('Longitude (°W)')
text(-48.47,1.5,'(b)','color','w','FontSize',30,'FontWeight','Bold')
set(gca,'FontSize',25)%,'FontWeight','Bold')

s3 = subplot('Position',plot_location(3,:));
hold on
rectangle('Position',[-48.5 0 0.5 3000],'FaceColor','k','EdgeColor','k')
pcolor(repmat(lon_lc,32,1),squeeze(nanmean(depth_lc,1))/1000,...
    squeeze(nanmean(salt_lc,1)))
%pcolor(repmat(lon_lc,32,1),squeeze(depth_lc(1,:,:)),squeeze(density_lc(1,:,:)))
shading interp
set(gca,'ydir','reverse')
colormap(s3,flipud(colormap_wave(3:82,:)))
c = colorbar('Orientation','Vertical','Position',cbar_location(3,:),...
    'Color',[1-eps,1,1],'Fontsize',10); %colorbar over land
%c = colorbar('Location','manual','Orientation','Horizontal','Position',[0.17 0.7 0.25 0.03],'Color','w','Fontsize',10)%colorbar over land
%c = colorbar;
c.Label.String = 'Salinity';
%ylabel(h,'SST difference (^\circC)','Fontsize',14)
set(c,'Fontsize',30)
caxis([33 35])
ylim([0 3])
xlim([-48.5 -48])
%set(gca,'xtick',[],'ytick',[])
set(gca,'xtick',[-48.5,-48.25,-48],...
    'xticklabels',{'48.5','48.25','48'})
for dd = 1:size(depth_lc,2)
    line(lon_lc,squeeze(nanmean(depth_lc(:,dd,:),1))/1000,...
        'color',rgb('Gray'),'LineWidth',0.5)
end
line(lon_lc,squeeze(nanmean(depth_lc(:,18,:),1))/1000,'color','k','LineWidth',2)
line(lon_lc,squeeze(nanmean(depth_lc(:,19,:),1))/1000,'color','k','LineWidth',2)
%line(lon_lc,squeeze(nanmean(depth_lc(:,20,:),1))/1000,'color','k','LineWidth',2)
%ylabel('Depth (m)')
xlabel('Longitude (°W)')
text(-48.47,1.5,'(c)','color','w','FontSize',30,'FontWeight','Bold')
set(gca,'FontSize',25)%,'FontWeight','Bold')

s4 = subplot('Position',plot_location(4,:));
hold on
rectangle('Position',[-48.5 0 0.5 3000],'FaceColor','k','EdgeColor','k')
pcolor(repmat(lon_lc,32,1),squeeze(nanmean(depth_lc,1))/1000,...
    squeeze(nanmean(ucross_lc,1)))
%pcolor(repmat(lon_lc,32,1),squeeze(depth_lc(1,:,:)),squeeze(density_lc(1,:,:)))
shading interp
set(gca,'ydir','reverse')
colormap(s4,flipud(colormap_wave2(4:70,:)))
c = colorbar('Orientation','Vertical','Position',cbar_location(4,:),...
    'Color',[1-eps,1,1],'Fontsize',10); %colorbar over land
%c = colorbar('Location','manual','Orientation','Horizontal','Position',[0.17 0.7 0.25 0.03],'Color','w','Fontsize',10)%colorbar over land
%c = colorbar;
c.Label.String = 'U across (m s^{-1})';
%ylabel(h,'SST difference (^\circC)','Fontsize',14)
set(c,'Fontsize',30)
caxis([-0.37 0.3])
ylim([0 3])
xlim([-48.5 -48])
%set(gca,'xtick',[],'ytick',[])
set(gca,'xtick',[-48.5,-48.25,-48],...
    'xticklabels',{'48.5','48.25','48'})
for dd = 1:size(depth_lc,2)
    line(lon_lc,squeeze(nanmean(depth_lc(:,dd,:),1))/1000,...
        'color',rgb('Gray'),'LineWidth',0.5)
end
line(lon_lc,squeeze(nanmean(depth_lc(:,18,:),1))/1000,'color','k','LineWidth',2)
line(lon_lc,squeeze(nanmean(depth_lc(:,19,:),1))/1000,'color','k','LineWidth',2)
%line(lon_lc,squeeze(nanmean(depth_lc(:,20,:),1))/1000,'color','k','LineWidth',2)
%ylabel('Depth (m)')
xlabel('Longitude (°W)')
text(-48.47,1.5,'(d)','color','w','FontSize',30,'FontWeight','Bold')
set(gca,'FontSize',25)%,'FontWeight','Bold')
set(gcf, 'Color', 'w')
saveas(gcf,[dir_output,'manuscript_v3/TGB2_Figure3_v3'],'png')

% Statistics

% Mean T (all-time)
nanmean(nanmean(temp_lc(:,19,:),3))% 4.8158
nanstd(nanmean(temp_lc(:,19,:),3))% 0.3095
% Mean S (all-time)
nanmean(nanmean(salt_lc(:,19,:),3))% 34.8802
nanstd(nanmean(salt_lc(:,19,:),3))% 0.0361
% Mean thickness (all-time)
nanmean(nanmean(thick_lc(:,19,:),3))% 63.1416
nanstd(nanmean(thick_lc(:,19,:),3))% 6.1839
% Mean speed (all-time)
nanmean(nanmean(ucross_lc(:,19,:),3))% 0.0767
nanstd(nanmean(ucross_lc(:,19,:),3))% 0.0430

% Yearly
temp_lc_yr = NaN(10,2);
salt_lc_yr = NaN(10,2);
thick_lc_yr = NaN(10,2);
ucross_lc_yr = NaN(10,2);
for yy = 1:10
    temp_lc_yr(yy,1) = nanmean(nanmean(temp_lc(yy*12-11:yy*12,19,:),3));
    temp_lc_yr(yy,2) = nanstd(nanmean(temp_lc(yy*12-11:yy*12,19,:),3));
    salt_lc_yr(yy,1) = nanmean(nanmean(salt_lc(yy*12-11:yy*12,19,:),3));
    salt_lc_yr(yy,2) = nanstd(nanmean(salt_lc(yy*12-11:yy*12,19,:),3));
    thick_lc_yr(yy,1) = nanmean(nanmean(thick_lc(yy*12-11:yy*12,19,:),3));
    thick_lc_yr(yy,2) = nanstd(nanmean(thick_lc(yy*12-11:yy*12,19,:),3));
    ucross_lc_yr(yy,1) = nanmean(nanmean(ucross_lc(yy*12-11:yy*12,19,:),3));
    ucross_lc_yr(yy,2) = nanstd(nanmean(ucross_lc(yy*12-11:yy*12,19,:),3));
end

% Monthly
temp_lc_mo = NaN(12,2);
salt_lc_mo = NaN(12,2);
thick_lc_mo = NaN(12,2);
ucross_lc_mo = NaN(12,2);
for mm = 1:12
    temp_lc_mo(mm,1) = nanmean(nanmean(temp_lc(mm:12:9*12+mm,19,:),3));
    temp_lc_mo(mm,2) = nanstd(nanmean(temp_lc(mm:12:9*12+mm,19,:),3));
    salt_lc_mo(mm,1) = nanmean(nanmean(salt_lc(mm:12:9*12+mm,19,:),3));
    salt_lc_mo(mm,2) = nanstd(nanmean(salt_lc(mm:12:9*12+mm,19,:),3));
    thick_lc_mo(mm,1) = nanmean(nanmean(thick_lc(mm:12:9*12+mm,19,:),3));
    thick_lc_mo(mm,2) = nanstd(nanmean(thick_lc(mm:12:9*12+mm,19,:),3));
    ucross_lc_mo(mm,1) = nanmean(nanmean(ucross_lc(mm:12:9*12+mm,19,:),3));
    ucross_lc_mo(mm,2) = nanstd(nanmean(ucross_lc(mm:12:9*12+mm,19,:),3));
end

% Fig 3 Statistics Layer 19 (1993-2017)

dir_input = '/Users/afonso/Documents/Research/TGB2/Data/Eulerian/';
dir_output = '/Users/afonso/Documents/Research/TGB2/output/';

lat = ncread([dir_output,'ATLg008_eulerian_layer19.nc'],'latitude');
lon = ncread([dir_output,'ATLg008_eulerian_layer19.nc'],'longitude');
lat = lat';
lon = lon' - 360;

% Define lat and lon at LC
lat_lc = linspace(45.5,45,7);
lon_lc = linspace(-48.5,-48,7);

% length of section
gsw_distance([-48.5,-48],[45.5,45],0) / 1000

% Create variable for time
year = 1993:2017;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(300,1);
for yy = 1:25
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

d1 = dir('/Users/afonso/Documents/Research/TGB2/Data/Eulerian');
d1 = d1(4:end);
thick = [];
depth = [];
temp = [];
salt = [];
for ii=1:length(d1)
    load([dir_input,d1(ii).name]);
    thick = cat(3,thick,hh(817:1033,263:727));
    depth = cat(3,depth,zu(817:1033,263:727));
    temp = cat(3,temp,tt(817:1033,263:727));
    salt = cat(3,salt,ss(817:1033,263:727));
    clear hh ss tt zu
end
clear lat_aux lon_aux ii d1
temp(temp >= 10^10) = NaN;% remove values on land
temp(thick <= 1) = NaN;% remove values on shelf
salt(salt >= 10^10) = NaN;% remove values on land
salt(thick <= 1) = NaN;% remove values on shelf

% Density
SA = gsw_SA_from_SP(salt,1000,lat_lc(1),lon_lc(1));
TC = gsw_CT_from_t(salt,temp,1000);
density = gsw_sigma2(SA,TC);
clear SA TC

% U
u = ncread([dir_output,'ATLg_008_layer19_monthly.nc'],'uvel');
u = u(263:727,817:1033,:);
u = permute(u,[2,1,3]);
u(u >= 10^10) = NaN;% remove values on land
u(thick <= 1) = NaN;% remove values on shelf
u = u .* 100;% cm/s

% V
v = ncread([dir_output,'ATLg_008_layer19_monthly.nc'],'vvel');
v = v(263:727,817:1033,:);
v = permute(v,[2,1,3]);
v(v >= 10^10) = NaN;% remove values on land
v(thick <= 1) = NaN;% remove values on shelf
v = v .* 100;% cm/s

% interpolate variables to LC line
thick_lc = NaN(length(lat_lc),size(thick,3));
depth_lc = NaN(length(lat_lc),size(thick,3));
density_lc = NaN(length(lat_lc),size(thick,3));
temp_lc = NaN(length(lat_lc),size(thick,3));
salt_lc = NaN(length(lat_lc),size(thick,3));
u_lc = NaN(length(lat_lc),size(thick,3));
v_lc = NaN(length(lat_lc),size(thick,3));
for tt = 1:size(thick,3)
    thick_lc(:,tt) = interp2(lon(1:end-1,:),lat(1:end-1,:),thick(1:end-1,:,tt),...
        lon_lc,lat_lc);
    depth_lc(:,tt) = interp2(lon(1:end-1,:),lat(1:end-1,:),squeeze(depth(1:end-1,:,tt)),...
        lon_lc,lat_lc);
    density_lc(:,tt) = interp2(lon(1:end-1,:),lat(1:end-1,:),squeeze(density(1:end-1,:,tt)),...
        lon_lc,lat_lc);
    temp_lc(:,tt) = interp2(lon(1:end-1,:),lat(1:end-1,:),squeeze(temp(1:end-1,:,tt)),...
        lon_lc,lat_lc);
    salt_lc(:,tt) = interp2(lon(1:end-1,:),lat(1:end-1,:),squeeze(salt(1:end-1,:,tt)),...
        lon_lc,lat_lc);
    u_lc(:,tt) = interp2(lon(1:end-1,:),lat(1:end-1,:),squeeze(u(1:end-1,:,tt)),...
        lon_lc,lat_lc);
    v_lc(:,tt) = interp2(lon(1:end-1,:),lat(1:end-1,:),squeeze(v(1:end-1,:,tt)),...
        lon_lc,lat_lc);
end
clear tt

% Rotate velocity

uv_lc = sqrt(u_lc.^2 + v_lc.^2);

ucross_lc = NaN(size(u_lc));
ualong_lc = NaN(size(u_lc));
for tt = 1:size(u_lc,2)
        for ll = 2:size(u_lc,1)
            
            % direction of input vector (angle relative to north, increases clockwise)
            if u_lc(ll,tt) >= 0 & v_lc(ll,tt) >= 0
                theta_in = atand(u_lc(ll,tt)/v_lc(ll,tt));
                %theta_in_aux = theta_in
            elseif u_lc(ll,tt) >= 0 & v_lc(ll,tt) < 0
                theta_in = 180 + atand(u_lc(ll,tt)/v_lc(ll,tt));
                %theta_in_aux = theta_in - 90
            elseif u_lc(ll,tt) < 0 & v_lc(ll,tt) < 0
                theta_in = 180 + atand(u_lc(ll,tt)/v_lc(ll,tt));
                %theta_in_aux = theta_in - 180
            elseif u_lc(ll,tt) < 0 & v_lc(ll,tt) >= 0
                theta_in = 360 + atand(u_lc(ll,tt)/v_lc(ll,tt));
                %theta_in_aux = theta_in - 270
            end
            
            % direction of speed (angle relative to north, increases clockwise)
            % lat decreases and lon increases
            theta_rot = 180 + atand((lon_lc(end)-lon_lc(1))/(lat_lc(end)-lat_lc(1)));
            
            % direction of speed relative to the new coordinate system
            if theta_in > theta_rot
                theta_out = theta_in - theta_rot;
            else
                theta_out = 360 - theta_rot + theta_in;
            end
            
            ucross_lc(ll,tt) = uv_lc(ll,tt) * sind(theta_out);
            ualong_lc(ll,tt) = uv_lc(ll,tt) * cosd(theta_out);
            
            clear theta*
        end
end
clear tt ll

% Statistics

% Mean T (all-time)
nanmean(nanmean(temp_lc))% 4.9572
nanstd(nanmean(temp_lc))% 0.4149
% Mean S (all-time)
nanmean(nanmean(salt_lc))% 34.9075
nanstd(nanmean(salt_lc))% 0.0361
% Mean thickness (all-time)
nanmean(nanmean(thick_lc))% 63.9210
nanstd(nanmean(thick_lc))% 7.8016
% Mean speed (all-time)
nanmean(nanmean(ucross_lc))% 7.8356
nanstd(nanmean(ucross_lc))% 4.8943

% Yearly
temp_lc_yr = NaN(10,2);
salt_lc_yr = NaN(10,2);
thick_lc_yr = NaN(10,2);
ucross_lc_yr = NaN(10,2);
for yy = 1:25
    temp_lc_yr(yy,1) = nanmean(nanmean(temp_lc(:,yy*12-11:yy*12)));
    temp_lc_yr(yy,2) = nanstd(nanmean(temp_lc(:,yy*12-11:yy*12)));
    salt_lc_yr(yy,1) = nanmean(nanmean(salt_lc(:,yy*12-11:yy*12)));
    salt_lc_yr(yy,2) = nanstd(nanmean(salt_lc(:,yy*12-11:yy*12)));
    thick_lc_yr(yy,1) = nanmean(nanmean(thick_lc(:,yy*12-11:yy*12)));
    thick_lc_yr(yy,2) = nanstd(nanmean(thick_lc(:,yy*12-11:yy*12)));
    ucross_lc_yr(yy,1) = nanmean(nanmean(ucross_lc(:,yy*12-11:yy*12)));
    ucross_lc_yr(yy,2) = nanstd(nanmean(ucross_lc(:,yy*12-11:yy*12)));
end

% Monthly
temp_lc_mo = NaN(12,2);
salt_lc_mo = NaN(12,2);
thick_lc_mo = NaN(12,2);
ucross_lc_mo = NaN(12,2);
for mm = 1:12
    temp_lc_mo(mm,1) = nanmean(nanmean(temp_lc(:,mm:12:24*12+mm)));
    temp_lc_mo(mm,2) = nanstd(nanmean(temp_lc(:,mm:12:24*12+mm)));
    salt_lc_mo(mm,1) = nanmean(nanmean(salt_lc(:,mm:12:24*12+mm)));
    salt_lc_mo(mm,2) = nanstd(nanmean(salt_lc(:,mm:12:24*12+mm)));
    thick_lc_mo(mm,1) = nanmean(nanmean(thick_lc(:,mm:12:24*12+mm)));
    thick_lc_mo(mm,2) = nanstd(nanmean(thick_lc(:,mm:12:24*12+mm)));
    ucross_lc_mo(mm,1) = nanmean(nanmean(ucross_lc(:,mm:12:24*12+mm)));
    ucross_lc_mo(mm,2) = nanstd(nanmean(ucross_lc(:,mm:12:24*12+mm)));
end

%% Figure 4

% Mean probability map with 50 random trajectories

dir_output = '/Users/afonso/Documents/Research/TGB2/output/';

%file_field = [dir_output,'hycom_glo_surface1_particles1-1.nc'];
file_traj = [dir_output,'ATLg_008_layer19_LC80'];

ncdisp([file_traj,'.nc'])

pset.trajectory = ncread([file_traj,'.nc'],'trajectory');
pset.lon = ncread([file_traj,'.nc'],'lon');
pset.lat = ncread([file_traj,'.nc'],'lat');
pset.time = ncread([file_traj,'.nc'],'time');
pset.time2 = pset.time/(60*60*24)+datenum(1993,1,1,12,0,0);
clear pset.time

% Load bathymetry from ETOPO1
x = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','x');
lon_z = x(x<=-10 & x>=-80);
lon_z = lon_z(1:end);
y = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','y');
lat_z = y(y<=65 & y>=30);
lat_z = lat_z(1:end);
z_aux = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','z');
z = z_aux(x<=-10 & x>=-80,y<=65 & y>=30); z = z';
z = z(1:end,1:end);
clear x y z_aux

% Define run time
time_run = 1;%years
time_total = 25;%years

%Subset
time_length = floor(size(pset.time,1)*(time_run/time_total));
traj_length = floor(size(pset.time,2)*((time_total-time_run)/time_total));

lat_prob = 30:1:65;
lon_prob = -80:1:-10;

prob1 = NaN(length(lat_prob)-1,length(lon_prob)-1);
for ii = 1:length(lon_prob)-1
    disp([num2str(100*ii/(length(lon_prob)-1)),'%'])
    for jj = 1:length(lat_prob)-1
        prob1_aux = zeros(traj_length,1);
        for kk = 1:traj_length
            prob1_aux(kk,1) = sum(find(...
                pset.lat(1:time_length,kk) >= lat_prob(jj) & ...
                pset.lat(1:time_length,kk) < lat_prob(jj+1) & ...
                pset.lon(1:time_length,kk) >= lon_prob(ii) & ...
                pset.lon(1:time_length,kk) < lon_prob(ii+1)));
        end
        prob1(jj,ii) = 100 * sum(prob1_aux>0)/traj_length;
%             prob1(jj,ii) = 100*length(find(...
%                 pset.lat(1:time_length,1:traj_length) >= lat_prob(jj) & ...
%                 pset.lat(1:time_length,1:traj_length) < lat_prob(jj+1) & ...
%                 pset.lon(1:time_length,1:traj_length) >= lon_prob(ii) & ...
%                 pset.lon(1:time_length,1:traj_length) < lon_prob(ii+1)))/traj_length;
        
    end
end

prob2 = prob1;
prob2(prob2 == 0) = NaN;

% How many particles?
traj_number = 50;

traj_plot = randi([1 traj_length],traj_number,1);

cbar_bathy = [5/256,112/256,176/256;...
    54/256,144/256,192/256;...
    116/256,169/256,207/256;...
    166/256,189/256,219/256;...
    208/256,209/256,230/256;...
    236/256,231/256,242/256;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0];

colormap_aux = readtable([dir_output,'5w_BRgpb.csv'], 'HeaderLines',1);
colormap_wave = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'4w_ROTB.csv'], 'HeaderLines',1);
colormap_wave2 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'3-wave-yellow-grey-blue.csv'], 'HeaderLines',1);
colormap_wave3 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'5-step-melow-wave.csv'], 'HeaderLines',1);
colormap_wave4 = table2array(colormap_aux(:,3:5));
clear colormap_aux

figure('Position',[100,200,800,500])
%figure('Units','Inches','Position', [0, 0, 13.33, 7.5])
p1 = subplot('Position',[0.06 0.06 0.9 0.92]);
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30 65],'MapLonLimit',[-80 -10])
hold on
%pcolorm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end));
pcolorm(lat_z,lon_z,z);
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',14)
shading flat
colormap(p1,cbar_bathy)
%contourm(lat_z(1:5:end),lon_z(1:5:end),z(1:5:end,1:5:end),[0 0],'color','k','LineWidth',1.5)
%colormap(cmocean('ice',50))
%cb = colorbar;
caxis([-6000 6000])
alpha 0.4
plabel('PLabelLocation',10,'FontSize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%framem
tightmap

p2 = axes('Position',[0.06 0.06 0.9 0.92]);
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30 65],'MapLonLimit',[-80 -10])
hold on
pcolorm(lat_prob(1:end-1),lon_prob(1:end-1),prob2)
shading flat
caxis([0 41])
%colormap(crameri('bilbao',25))
colormap(p2,flipud(colormap_wave2(1:end-9,:)))
c = colorbar('Orientation','Horizontal','Position',[0.11 0.6 0.16 0.05],...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'Probability (%)';
set(c,'Fontsize',15,'Fontweight','bold')
alpha 0.95
contourm(lat_prob(1:end-1),lon_prob(1:end-1),prob2,[10 10],'color',rgb('snow'),'LineWidth',4)
for kk = 1:traj_number
    linem([pset.lat(1,traj_plot(kk)) pset.lat(time_length,traj_plot(kk))],...
        [pset.lon(1,traj_plot(kk)) pset.lon(time_length,traj_plot(kk))],...
        'color',rgb('Black'),'LineWidth',1.5)
end
% contourm(lat_z,lon_z,z,...
%     [-100 -100],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-1000 -1000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-3000 -3000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-4000 -4000],'color',[77/256 77/256 77/256],'linewidth',0.5)
geoshow('GSHHS_i_L1.shp','FaceColor',rgb('Black'))
%geoshow('landareas.shp','FaceColor','k')
plabel('PLabelLocation',10,'FontSize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%title('Probability map (Labrador Current, Forward 1 year, 2003-2012)')
set(p2,'color','none')
framem
tightmap
saveas(gcf,[dir_output,'manuscript_v3/TGB2_Figure4_v4'],'png')

%% Figure 5

% Probability map with 50 random trajectories (2008 minus 2004)

dir_output = '/Users/afonso/Documents/Research/TGB2/output/';

%file_field = [dir_output,'hycom_glo_surface1_particles1-1.nc'];
file_traj = [dir_output,'ATLg_008_layer19_LC80'];

ncdisp([file_traj,'.nc'])

pset.trajectory = ncread([file_traj,'.nc'],'trajectory');
pset.lon = ncread([file_traj,'.nc'],'lon');
pset.lat = ncread([file_traj,'.nc'],'lat');
pset.time = ncread([file_traj,'.nc'],'time');
pset.time2 = pset.time/(60*60*24)+datenum(1993,1,1,12,0,0);

% Load bathymetry from ETOPO1
x = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','x');
lon_z = x(x<=-10 & x>=-80);
lon_z = lon_z(1:end);
y = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','y');
lat_z = y(y<=65 & y>=30);
lat_z = lat_z(1:end);
z_aux = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','z');
z = z_aux(x<=-10 & x>=-80,y<=65 & y>=30); z = z';
z = z(1:end,1:end);
clear x y z_aux

month_deploy = str2num(datestr(pset.time2(1,:),'mm'));
year_deploy = str2num(datestr(pset.time2(1,:),'yyyy'));

% Define run time
time_run = 1;%years
time_total = 25;%years
%Subset
time_length = floor(size(pset.time,1)*(time_run/time_total));
traj_length = floor(size(pset.time,2)*((time_total-time_run)/time_total));

lat_prob = 30:1:65;
lon_prob = -80:1:-10;

prob2 = NaN(length(lat_prob)-1,length(lon_prob)-1);
for ii = 1:length(lon_prob)-1
    disp([num2str(100*ii/(length(lon_prob)-1)),'%'])
    for jj = 1:length(lat_prob)-1
            traj_pos = find(year_deploy == 2004);
            traj_length = length(pset.time2(year_deploy == 2004));
            prob2_aux = zeros(traj_length,1);
            for kk = 1:traj_length
                prob2_aux(kk,1) = sum(find(...
                    pset.lat(1:time_length,traj_pos(kk,1)) >= lat_prob(jj) & ...
                    pset.lat(1:time_length,traj_pos(kk,1)) < lat_prob(jj+1) & ...
                    pset.lon(1:time_length,traj_pos(kk,1)) >= lon_prob(ii) & ...
                    pset.lon(1:time_length,traj_pos(kk,1)) < lon_prob(ii+1)));
            end
            prob2(jj,ii) = 100 * sum(prob2_aux>0,1)/traj_length;
            clear traj_pos traj_length
%             prob1(jj,ii) = 100*length(find(...
%                 pset.lat(1:time_length,1:traj_length) >= lat_prob(jj) & ...
%                 pset.lat(1:time_length,1:traj_length) < lat_prob(jj+1) & ...
%                 pset.lon(1:time_length,1:traj_length) >= lon_prob(ii) & ...
%                 pset.lon(1:time_length,1:traj_length) < lon_prob(ii+1)))/traj_length;
        
    end
end
clear ii jj prob2_aux

prob3 = NaN(length(lat_prob)-1,length(lon_prob)-1);
for ii = 1:length(lon_prob)-1
    disp([num2str(100*ii/(length(lon_prob)-1)),'%'])
    for jj = 1:length(lat_prob)-1
            traj_pos = find(year_deploy == 2008);
            traj_length = length(pset.time2(year_deploy == 2008));
            prob3_aux = zeros(traj_length,1);
            for kk = 1:traj_length
                prob3_aux(kk,1) = sum(find(...
                    pset.lat(1:time_length,traj_pos(kk,1)) >= lat_prob(jj) & ...
                    pset.lat(1:time_length,traj_pos(kk,1)) < lat_prob(jj+1) & ...
                    pset.lon(1:time_length,traj_pos(kk,1)) >= lon_prob(ii) & ...
                    pset.lon(1:time_length,traj_pos(kk,1)) < lon_prob(ii+1)));
            end
            prob3(jj,ii) = 100 * sum(prob3_aux>0,1)/traj_length;
            clear traj_pos traj_length
%             prob1(jj,ii) = 100*length(find(...
%                 pset.lat(1:time_length,1:traj_length) >= lat_prob(jj) & ...
%                 pset.lat(1:time_length,1:traj_length) < lat_prob(jj+1) & ...
%                 pset.lon(1:time_length,1:traj_length) >= lon_prob(ii) & ...
%                 pset.lon(1:time_length,1:traj_length) < lon_prob(ii+1)))/traj_length;
        
    end
end
clear ii jj prob3_aux

prob4 = prob3 - prob2;
prob4a = prob4;
prob4a(prob4 == 0) = NaN;
prob4a(prob4a < 1 & prob4a > -1) = NaN;

% How many particles?
traj_number = 100;

traj_plot_2004 = randi([find(year_deploy==2004,1,'first') ...
    find(year_deploy==2004,1,'last')],traj_number,1);
traj_plot_2008 = randi([find(year_deploy==2008,1,'first') ...
    find(year_deploy==2008,1,'last')],traj_number,1);

% probability map (March-August)
prob5 = NaN(length(lat_prob)-1,length(lon_prob)-1);
for ii = 1:length(lon_prob)-1
    disp([num2str(100*ii/(length(lon_prob)-1)),'%'])
    for jj = 1:length(lat_prob)-1
            traj_pos = find(month_deploy >= 3 & month_deploy <= 8);
            traj_length = length(pset.time2(month_deploy >= 3 & month_deploy <= 8));
            prob5_aux = zeros(traj_length,1);
            for kk = 1:traj_length
                prob5_aux(kk,1) = sum(find(...
                    pset.lat(1:time_length,traj_pos(kk,1)) >= lat_prob(jj) & ...
                    pset.lat(1:time_length,traj_pos(kk,1)) < lat_prob(jj+1) & ...
                    pset.lon(1:time_length,traj_pos(kk,1)) >= lon_prob(ii) & ...
                    pset.lon(1:time_length,traj_pos(kk,1)) < lon_prob(ii+1)));
            end
            prob5(jj,ii) = 100 * sum(prob5_aux>0,1)/traj_length;
            clear traj_pos traj_length
%             prob1(jj,ii) = 100*length(find(...
%                 pset.lat(1:time_length,1:traj_length) >= lat_prob(jj) & ...
%                 pset.lat(1:time_length,1:traj_length) < lat_prob(jj+1) & ...
%                 pset.lon(1:time_length,1:traj_length) >= lon_prob(ii) & ...
%                 pset.lon(1:time_length,1:traj_length) < lon_prob(ii+1)))/traj_length;
        
    end
end
clear ii jj prob5_aux

prob6 = NaN(length(lat_prob)-1,length(lon_prob)-1);
for ii = 1:length(lon_prob)-1
    disp([num2str(100*ii/(length(lon_prob)-1)),'%'])
    for jj = 1:length(lat_prob)-1
            traj_pos = find(month_deploy <= 2 | month_deploy >= 9);
            traj_length = length(pset.time2(month_deploy <= 2 | month_deploy >= 9));
            prob6_aux = zeros(traj_length,1);
            for kk = 1:traj_length
                prob6_aux(kk,1) = sum(find(...
                    pset.lat(1:time_length,traj_pos(kk,1)) >= lat_prob(jj) & ...
                    pset.lat(1:time_length,traj_pos(kk,1)) < lat_prob(jj+1) & ...
                    pset.lon(1:time_length,traj_pos(kk,1)) >= lon_prob(ii) & ...
                    pset.lon(1:time_length,traj_pos(kk,1)) < lon_prob(ii+1)));
            end
            prob6(jj,ii) = 100 * sum(prob6_aux>0,1)/traj_length;
            clear traj_pos traj_length
%             prob1(jj,ii) = 100*length(find(...
%                 pset.lat(1:time_length,1:traj_length) >= lat_prob(jj) & ...
%                 pset.lat(1:time_length,1:traj_length) < lat_prob(jj+1) & ...
%                 pset.lon(1:time_length,1:traj_length) >= lon_prob(ii) & ...
%                 pset.lon(1:time_length,1:traj_length) < lon_prob(ii+1)))/traj_length;
        
    end
end
clear ii jj prob6_aux

prob7 = prob5 - prob6;
prob7a = prob7;
prob7a(prob7 == 0) = NaN;
prob7a(prob7a < 1 & prob7a > -1) = NaN;

% How many particles?
traj_number = 100;

traj_plot_summer = randi([find(month_deploy >= 3 & month_deploy <= 8,1,'first') ...
    find(month_deploy >= 3 & month_deploy <= 8,1,'last')],traj_number,1);
traj_plot_winter = randi([find(month_deploy <= 2 | month_deploy >= 9,1,'first') ...
    find(month_deploy <= 2 | month_deploy >= 9,1,'last')],traj_number,1);

load([dir_output,'etopo1_TGB.mat'],'lat_z','lon_z','z')

figure('Position', [100, 200, 800, 600])
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
hold on
c1_aux = contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
linem([36 36],[360-75 360-53],'k')
linem([46 46],[360-75 360-53],'k')
linem([36 46],[360-75 360-75],'k')
linem([36 46],[360-53 360-53],'k')
geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])

close all

c1 = NaN(2,1);
j = 1;
for i = 2:length(c1_aux)
    if c1_aux(1,i)>=360-55 & c1_aux(1,i)<=360-52 & c1_aux(2,i)>=36 & c1_aux(2,i)<=46
        c1(:,j) = c1_aux(:,i);
        if j>1
            if abs(c1_aux(1,i)-c1(1,j-1))<0.2
                j = j+1;
            else
                c1(:,j) = NaN;
            end
        else
            j = j+1;
        end
    end
end
clear i j

if isnan(c1(1,end))
    c1 = c1(:,1:end-1);% lon,lat of 1000m isobath
end

c1(3,:) = zeros;
c1_dist = gsw_distance(c1(1,:),c1(2,:));
for i = 2:length(c1)
    c1(3,i) = sum(c1_dist(1:i-1));
end
clear c1_aux c1_dist i

% Line 1 degree to the south of the 1,000-m isobath

c2 = c1;
c2(2,:) = c1(2,:) - 1;
c3 = c1;
c3(2,:) = c1(2,:) - 2;

% figure('Position', [100, 200, 800, 600])
% axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
% hold on
% contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem(c1(2,:),c1(1,:),'k','linewidth',2)
% linem(c2(2,:),c2(1,:),'--k','linewidth',1)
% linem(c3(2,:),c3(1,:),'k','linewidth',2)
% %contourm(c1(2,:),c1(2,:),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem([36 36],[360-75 360-52],'k')
% linem([46 46],[360-75 360-52],'k')
% linem([36 46],[360-75 360-75],'k')
% linem([36 46],[360-52 360-52],'k')
% geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
% 
% close all

close all

% Define run time
time_run = 1;%years
time_total = 25;%years
%Subset
time_length = floor(size(pset.time,1)*(time_run/time_total));
traj_length = floor(size(pset.time,2)*((time_total-time_run)/time_total));

traj_slope = NaN(traj_length,1);
for ii = 1:traj_length
    traj_slope(ii,1) = sum(inpolygon(pset.lon(1:time_length,ii),pset.lat(1:time_length,ii),...
        [c1(1,:)-360 fliplr(c3(1,:)-360) c1(1,1)-360],...
        [c1(2,:) fliplr(c3(2,:)) c1(2,1)]));
end

% Number of particles per release time (time series)
slope_releasetime = pset.time2(1,traj_slope>=9);
slope_releasetime_yr = str2num(datestr(slope_releasetime,'yyyy'));
slope_releasetime_mo = str2num(datestr(slope_releasetime,'mm'));
perc1_slope = 100*length(slope_releasetime)/traj_length; % 15% of floats make it to slope

yr = 1993:2016;
slope_releasetime_yr1 = NaN(length(yr),1);
for ii = 1:length(yr)
    slope_releasetime_yr1(ii,1) = sum(slope_releasetime_yr==yr(ii));
end
perc_yr_slope = 100*slope_releasetime_yr1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(yr)
    perc_yr1_slope(ii,1) = 100*slope_releasetime_yr1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii));
end

mo = 1:12;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

slope_releasetime_mo1 = NaN(length(mo),1);
for ii = 1:length(mo)
    slope_releasetime_mo1(ii,1) = sum(slope_releasetime_mo==mo(ii));
end
perc_mo_slope = 100*slope_releasetime_mo1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(mo)
    perc_mo1_slope(ii,1) = 100*slope_releasetime_mo1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(ii));
end

% Area along the NWA Slope
lat_nwa = pset.lat(1:time_length,traj_slope>=9);
lon_nwa = pset.lon(1:time_length,traj_slope>=9);
time2_nwa = pset.time2(1:time_length,traj_slope>=9);

% Releases per month
slope_release_number = NaN(288,4);
for ii = 1:length(yr)
    for jj = 1:length(mo)
        slope_release_number(12*(ii-1)+jj,1) = yr(ii);
        slope_release_number(12*(ii-1)+jj,2) = mo(jj);
        slope_release_number(12*(ii-1)+jj,3) = ...
            sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii)...
            & str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(jj));
        slope_release_number(12*(ii-1)+jj,4) = sum(str2num(datestr(time2_nwa(1,:),'mm'))...
            ==mo(jj) & str2num(datestr(time2_nwa(1,:),'yyyy'))==yr(ii));
    end
end




% Define run time
time_run = 1;%years
time_total = 25;%years
%Subset
time_length = floor(size(pset.time,1)*(time_run/time_total));
traj_length = floor(size(pset.time,2)*((time_total-time_run)/time_total));

traj_nac = NaN(traj_length,1);
for ii = 1:traj_length
%    traj_nac(ii,1) = sum(inpolygon(pset.lon(1:time_length,ii),pset.lat(1:time_length,ii),...
%        [-43 -43 -40 -40 -43],[49 51 51 49 49]));
    traj_nac(ii,1) = sum(inpolygon(pset.lon(1:time_length,ii),pset.lat(1:time_length,ii),...
        [-44 -44 -41 -41 -44],[44 46 46 44 44]));
end

% Number of particles per release time (time series)
nac_releasetime = pset.time2(1,traj_nac>=3);
nac_releasetime_yr = str2num(datestr(nac_releasetime,'yyyy'));
nac_releasetime_mo = str2num(datestr(nac_releasetime,'mm'));
perc1_nac = 100*length(nac_releasetime)/traj_length; % 15% of floats make it to slope

yr = 1993:2016;
nac_releasetime_yr1 = NaN(length(yr),1);
for ii = 1:length(yr)
    nac_releasetime_yr1(ii,1) = sum(nac_releasetime_yr==yr(ii));
end
perc_yr_nac = 100*nac_releasetime_yr1./length(nac_releasetime); % 15% of floats make it to slope

for ii = 1:length(yr)
    perc_yr1_nac(ii,1) = 100*nac_releasetime_yr1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii));
end

mo = 1:12;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

nac_releasetime_mo1 = NaN(length(mo),1);
for ii = 1:length(mo)
    nac_releasetime_mo1(ii,1) = sum(nac_releasetime_mo==mo(ii));
end
perc_mo_nac = 100*nac_releasetime_mo1./length(nac_releasetime); % 15% of floats make it to slope

for ii = 1:length(mo)
    perc_mo1_nac(ii,1) = 100*nac_releasetime_mo1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(ii));
end

% Area along the NWA Slope
lat_nac = pset.lat(1:time_length,traj_nac>=3);
lon_nac = pset.lon(1:time_length,traj_nac>=3);
time2_nac = pset.time2(1:time_length,traj_nac>=3);

% Releases per month
nac_release_number = NaN(288,4);
for ii = 1:length(yr)
    for jj = 1:length(mo)
        nac_release_number(12*(ii-1)+jj,1) = yr(ii);
        nac_release_number(12*(ii-1)+jj,2) = mo(jj);
        nac_release_number(12*(ii-1)+jj,3) = ...
            sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii)...
            & str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(jj));
        nac_release_number(12*(ii-1)+jj,4) = sum(str2num(datestr(time2_nac(1,:),'mm'))...
            ==mo(jj) & str2num(datestr(time2_nac(1,:),'yyyy'))==yr(ii));
    end
end


% % Load bathymetry from ETOPO1
x = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','x');
lon_z = x(x<=-10 & x>=-80);
lon_z = lon_z(1:10:end);
y = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','y');
lat_z = y(y<=65 & y>=30);
lat_z = lat_z(1:10:end);
z_aux = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','z');
z = z_aux(x<=-10 & x>=-80,y<=65 & y>=30); z = z';
z = z(1:10:end,1:10:end);
clear x y z_aux

cbar_bathy = [5/256,112/256,176/256;...
    54/256,144/256,192/256;...
    116/256,169/256,207/256;...
    166/256,189/256,219/256;...
    208/256,209/256,230/256;...
    236/256,231/256,242/256;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0];

aa = find(str2num(datestr(time2_nwa(1,:),'yyyy'))==2004);
bb = find(str2num(datestr(time2_nac(1,:),'yyyy'))==2004);
aa1 = find(str2num(datestr(time2_nwa(1,:),'yyyy'))==2008);
bb1 = find(str2num(datestr(time2_nac(1,:),'yyyy'))==2008);

colormap_aux = readtable([dir_output,'blue-orange-div.csv'], 'HeaderLines',1);
colormap_wave = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'4w_ROTB.csv'], 'HeaderLines',1);
colormap_wave2 = table2array(colormap_aux(:,3:5));
clear colormap_aux


cbar_bathy = [5/256,112/256,176/256;...
    54/256,144/256,192/256;...
    116/256,169/256,207/256;...
    166/256,189/256,219/256;...
    208/256,209/256,230/256;...
    236/256,231/256,242/256;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0];

% Statistics

% How many particles moved southward upon deployment? 
lat1 = pset.lat(1,1:traj_length);% lat at t=1
lat2 = pset.lat(2,1:traj_length);% lat at t=2
lat3 = pset.lat(3,1:traj_length);% lat at t=2
a = find(lat2>lat1); % instances where particles moved northward
a1 = find(lat3>lat1); % instances where particles moved northward
100 * length(a) / length(lat1) % 19.7108% (1993-2016), 10.9479% (2003-2011) northward
100 - (100 * length(a) / length(lat1)) % 80.2892% (1993-2016), 89.0521% (2003-2011 southward


% position of northward flow (mostly offshore)
a3 = NaN(80,1);
for ii = 1:80
    a3(ii,1) = length(find(rem(a,80) == (ii-1)));
end
clear lat1 lat2 lat3 a a1 a3

% Total distance traveled by each particle in 1 year
lat_diff = pset.lat(time_length,1:traj_length)' - pset.lat(1,1:traj_length)';
lon_diff = pset.lon(time_length,1:traj_length)' - pset.lon(1,1:traj_length)';
particle_dist = NaN(traj_length,time_length);
for ii = 1:traj_length
    for jj = 1:time_length-1
        particle_dist(ii,jj) = gsw_distance([pset.lat(jj,ii) pset.lat(jj+1,ii)],...
            [pset.lon(jj,ii) pset.lon(jj+1,ii)],0);
    end
end
% (2003-2011) Northeast: 3,923 +- 798 km (45.6%, 12,028 particles)
% (1993-2016) Northeast: 5,042 +- 1,204 km (47.0%, 33,011 particles)
length(particle_dist(lat_diff > 0 & lon_diff > 0,1))
100*length(particle_dist(lat_diff > 0 & lon_diff > 0,1))/traj_length
nanmean(nansum(particle_dist(lat_diff > 0 & lon_diff > 0,:),2))
nanstd(nansum(particle_dist(lat_diff > 0 & lon_diff > 0,:),2))

% (2003-2011) Southeast: 3,365 +- 659 km (38.4%, 10,123 particles)
% (1993-2016) Southeast: 4,085 +- 947 km (38.1%, 26,710 particles)
length(particle_dist(lat_diff < 0 & lon_diff > 0))
100*length(particle_dist(lat_diff < 0 & lon_diff > 0))/traj_length
nanmean(nansum(particle_dist(lat_diff < 0 & lon_diff > 0,:),2))
nanstd(nansum(particle_dist(lat_diff < 0 & lon_diff > 0,:),2))

% (2003-2011) Southwest: 2,618 +- 585 km (14.3%, 3,759 particles)
% (1993-2016) Southwest: 3,234 +- 783 km (13.1%, 9,194 particles)
length(particle_dist(lat_diff < 0 & lon_diff < 0))
100*length(particle_dist(lat_diff < 0 & lon_diff < 0))/traj_length
nanmean(nansum(particle_dist(lat_diff < 0 & lon_diff < 0,:),2))
nanstd(nansum(particle_dist(lat_diff < 0 & lon_diff < 0,:),2))

% (2003-2011) Northwest: 2,908 +- 995 km (1.7%, 442 particles)
% (1993-2016) Northwest: 3,655 +- 1,475 km (1.8%, 1,245 particles)
length(particle_dist(lat_diff > 0 & lon_diff < 0))
100*length(particle_dist(lat_diff > 0 & lon_diff < 0))/traj_length
nanmean(nansum(particle_dist(lat_diff > 0 & lon_diff < 0,:),2))
nanstd(nansum(particle_dist(lat_diff > 0 & lon_diff < 0,:),2))

% Distance between initial and final positions
% Northeast: 1,183 +- 565 km (12,028 particles)
length(particle_dist(lat_diff > 0 & lon_diff > 0,1))
nanmean(particle_dist(lat_diff > 0 & lon_diff > 0,1))
nanstd(particle_dist(lat_diff > 0 & lon_diff > 0,1))

% Southeast: 863 +- 441 km (10,123 particles)
length(particle_dist(lat_diff < 0 & lon_diff > 0))
nanmean(particle_dist(lat_diff < 0 & lon_diff > 0))
nanstd(particle_dist(lat_diff < 0 & lon_diff > 0))

% Southwest: 769 +- 415 km (3,759 particles)
length(particle_dist(lat_diff < 0 & lon_diff < 0))
nanmean(particle_dist(lat_diff < 0 & lon_diff < 0))
nanstd(particle_dist(lat_diff < 0 & lon_diff < 0))

% Northwest: 884 +- 336 km (442 particles)
length(particle_dist(lat_diff > 0 & lon_diff < 0))
nanmean(particle_dist(lat_diff > 0 & lon_diff < 0))
nanstd(particle_dist(lat_diff > 0 & lon_diff < 0))

plot_location = [0.06 0.53 0.45 0.44;...
    0.53 0.53 0.45 0.44;...
    0.06 0.07 0.45 0.44;...
    0.53 0.07 0.45 0.44];


figure('Position',[100,100,820,700])
%figure('Units','Inches','Position', [0, 0, 13.33, 7.5])
p1 = subplot('Position',plot_location(1,:));
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30.1 60],'MapLonLimit',[-65 -25])
hold on
%pcolorm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end));
pcolorm(lat_z,lon_z,z);
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',14)
shading flat
colormap(p1,cbar_bathy)
%contourm(lat_z(1:5:end),lon_z(1:5:end),z(1:5:end,1:5:end),[0 0],'color','k','LineWidth',1.5)
%colormap(cmocean('ice',50))
%cb = colorbar;
caxis([-6000 6000])
alpha 0.4
plabel('PLabelLocation',10,'FontSize',18)
%mlabel('MLabelLocation',10,'MLabelParallel','south')
%framem
tightmap

p2 = axes('Position',plot_location(1,:));
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30.1 60],'MapLonLimit',[-65 -25])
hold on
pcolorm(lat_prob(1:end-1),lon_prob(1:end-1),prob7a)
shading flat
caxis([-32 32])
%colormap(crameri('bilbao',25))
colormap(p2,colormap_wave)
c = colorbar('Orientation','Horizontal','Position',[0.35 0.61 0.1 0.02],...
    'Color','k','Ticks',[-30 0 30]);%colorbar over land
c.Label.String = 'Probability (%)';
set(c,'Fontsize',18,'Fontweight','bold')
alpha 0.95
textm(36,-56.8,'Fall/Winter','Color',rgb('midnightblue'),'FontSize',21,'Fontweight','bold')
textm(54,-44,'Spring/Summer','Color',rgb('DarkRed'),'FontSize',21,'Fontweight','bold')
textm(58.7,-62,'(a)','Color',rgb('Black'),'FontSize',26,'Fontweight','bold')
% for kk = 1:traj_number
%     linem(pset.lat(1:time_length,traj_plot_2004(kk)),...
%         pset.lon(1:time_length,traj_plot_2004(kk)),...
%         'color',rgb('SteelBlue'),'LineWidth',1.5)
%     linem(pset.lat(1:time_length,traj_plot_2008(kk)),...
%         pset.lon(1:time_length,traj_plot_2008(kk)),...
%         'color',rgb('DarkRed'),'LineWidth',1.5)
% end
% contourm(lat_z,lon_z,z,...
%     [-100 -100],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-1000 -1000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-3000 -3000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-4000 -4000],'color',[77/256 77/256 77/256],'linewidth',0.5)
geoshow('GSHHS_i_L1.shp','FaceColor',rgb('Black'))
%geoshow('landareas.shp','FaceColor','k')
plabel('PLabelLocation',10,'FontSize',18)
%mlabel('MLabelLocation',10,'MLabelParallel','south')
%title('Probability map (Labrador Current, Forward 1 year, 2003-2012)')
set(p2,'color','none')
framem
tightmap

%figure('Units','Inches','Position', [0, 0, 13.33, 7.5])
p3 = subplot('Position',plot_location(2,:));
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30.1 60],'MapLonLimit',[-65 -25])
hold on
%pcolorm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end));
pcolorm(lat_z,lon_z,z);
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',14)
shading flat
colormap(p3,cbar_bathy)
%contourm(lat_z(1:5:end),lon_z(1:5:end),z(1:5:end,1:5:end),[0 0],'color','k','LineWidth',1.5)
%colormap(cmocean('ice',50))
%cb = colorbar;
caxis([-6000 6000])
alpha 0.4
plabel('off')
%mlabel('MLabelLocation',10,'MLabelParallel','south','FontSize',18)
%framem
tightmap

p4 = axes('Position',plot_location(2,:));
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30.1 60],'MapLonLimit',[-65 -25])
hold on
pcolorm(lat_prob(1:end-1),lon_prob(1:end-1),prob4a)
shading flat
caxis([-32 32])
%colormap(crameri('bilbao',25))
colormap(p4,colormap_wave)
c = colorbar('Orientation','Horizontal','Position',[0.82 0.61 0.1 0.02],...
    'Color','k','Ticks',[-30 0 30]);%colorbar over land
c.Label.String = 'Probability (%)';
set(c,'Fontsize',18,'Fontweight','bold')
alpha 0.95
textm(36,-56.8,'2004','Color',rgb('midnightblue'),'FontSize',24,'Fontweight','bold')
textm(54,-38,'2008','Color',rgb('DarkRed'),'FontSize',24,'Fontweight','bold')
textm(58.7,-62,'(b)','Color',rgb('Black'),'FontSize',26,'Fontweight','bold')
plotm([44 46 46 44 44],[-44 -44 -41 -41 -44],'k','linewidth',6) 
plotm([c1(2,:) fliplr(c3(2,:)) c1(2,1)],[c1(1,:) fliplr(c3(1,:)) c1(1,1)],'k','linewidth',6)

% for kk = 1:traj_number
%     linem(pset.lat(1:time_length,traj_plot_2004(kk)),...
%         pset.lon(1:time_length,traj_plot_2004(kk)),...
%         'color',rgb('SteelBlue'),'LineWidth',1.5)
%     linem(pset.lat(1:time_length,traj_plot_2008(kk)),...
%         pset.lon(1:time_length,traj_plot_2008(kk)),...
%         'color',rgb('DarkRed'),'LineWidth',1.5)
% end
% contourm(lat_z,lon_z,z,...
%     [-100 -100],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-1000 -1000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-3000 -3000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-4000 -4000],'color',[77/256 77/256 77/256],'linewidth',0.5)
geoshow('GSHHS_i_L1.shp','FaceColor',rgb('Black'))
%geoshow('landareas.shp','FaceColor','k')
%plabel('PLabelLocation',10,'FontSize',18)
%mlabel('MLabelLocation',10,'MLabelParallel','south')
%title('Probability map (Labrador Current, Forward 1 year, 2003-2012)')
set(p4,'color','none')
framem
tightmap

% 2004 vs 2008

p5 = subplot('Position',plot_location(3,:));
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30.1 60],'MapLonLimit',[-65 -25])
hold on
%pcolorm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end));
pcolorm(lat_z,lon_z,z);
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',14)
shading flat
colormap(p5,cbar_bathy)
%contourm(lat_z(1:5:end),lon_z(1:5:end),z(1:5:end,1:5:end),[0 0],'color','k','LineWidth',1.5)
%colormap(cmocean('ice',50))
%cb = colorbar;
caxis([-6000 6000])
alpha 0.4
plabel('PLabelLocation',10,'FontSize',18)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%framem
tightmap

p6 = axes('Position',plot_location(3,:));
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30.1 60],'MapLonLimit',[-65 -25])
hold on
for ii = 1:length(aa)
    linem(lat_nwa(:,aa(ii)),lon_nwa(:,aa(ii)),...
        'color',rgb('midnightblue'),'linewidth',0.5)
    scatterm(lat_nwa(1,aa(ii)),lon_nwa(1,aa(ii)),...
        50,rgb('midnightblue'))
    scatterm(lat_nwa(end,aa(ii)),lon_nwa(end,aa(ii)),...
        30,rgb('midnightblue'),'filled')
end
for ii = 1:length(bb)
    linem(lat_nac(:,bb(ii)),lon_nac(:,bb(ii)),...
        'color',rgb('FireBrick'),'linewidth',0.5)
    scatterm(lat_nac(1,bb(ii)),lon_nac(1,bb(ii)),...
        50,rgb('FireBrick'))
    scatterm(lat_nac(end,bb(ii)),lon_nac(end,bb(ii)),...
        30,rgb('FireBrick'),'filled')
end
scatterm(45.25,-48.25,500,'k','filled','v')
scatterm(45.25,-48.25,200,'w','filled','v')
scatterm(45.25,-48.25,50,'k','filled','v')
textm(59,-46.8,'2004','Color',rgb('midnightblue'),'FontSize',20,'Fontweight','bold')
textm(59,-62,'(c)','Color',rgb('Black'),'FontSize',26,'Fontweight','bold')
% for kk = 1:traj_number
%     linem(pset.lat(1:time_length,traj_plot_2004(kk)),...
%         pset.lon(1:time_length,traj_plot_2004(kk)),...
%         'color',rgb('SteelBlue'),'LineWidth',1.5)
%     linem(pset.lat(1:time_length,traj_plot_2008(kk)),...
%         pset.lon(1:time_length,traj_plot_2008(kk)),...
%         'color',rgb('DarkRed'),'LineWidth',1.5)
% end
% contourm(lat_z,lon_z,z,...
%     [-100 -100],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-1000 -1000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-3000 -3000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-4000 -4000],'color',[77/256 77/256 77/256],'linewidth',0.5)
geoshow('GSHHS_i_L1.shp','FaceColor',rgb('Black'))
%geoshow('landareas.shp','FaceColor','k')
plabel('PLabelLocation',10,'FontSize',18)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%title('Probability map (Labrador Current, Forward 1 year, 2003-2012)')
set(p6,'color','none')
framem
tightmap

%figure('Units','Inches','Position', [0, 0, 13.33, 7.5])
p7 = subplot('Position',plot_location(4,:));
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30.1 60],'MapLonLimit',[-65 -25])
hold on
%pcolorm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end));
pcolorm(lat_z,lon_z,z);
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',14)
shading flat
colormap(p7,cbar_bathy)
%contourm(lat_z(1:5:end),lon_z(1:5:end),z(1:5:end,1:5:end),[0 0],'color','k','LineWidth',1.5)
%colormap(cmocean('ice',50))
%cb = colorbar;
caxis([-6000 6000])
alpha 0.4
plabel('off')
mlabel('MLabelLocation',10,'MLabelParallel','south','FontSize',18)
%framem
tightmap

p8 = axes('Position',plot_location(4,:));
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30.1 60],'MapLonLimit',[-65 -25])
hold on
for ii = 1:length(aa1)
    linem(lat_nwa(:,aa1(ii)),lon_nwa(:,aa1(ii)),...
        'color',rgb('midnightblue'),'linewidth',0.5)
    scatterm(lat_nwa(1,aa1(ii)),lon_nwa(1,aa1(ii)),...
        50,rgb('midnightblue'))
    scatterm(lat_nwa(end,aa1(ii)),lon_nwa(end,aa1(ii)),...
        30,rgb('midnightblue'),'filled')
end
for ii = 1:length(bb1)
    linem(lat_nac(:,bb1(ii)),lon_nac(:,bb1(ii)),...
        'color',rgb('FireBrick'),'linewidth',0.5)
    scatterm(lat_nac(1,bb1(ii)),lon_nac(1,bb1(ii)),...
        50,rgb('FireBrick'))
    scatterm(lat_nac(end,bb1(ii)),lon_nac(end,bb1(ii)),...
        30,rgb('FireBrick'),'filled')
end
textm(59,-46.8,'2008','Color',rgb('DarkRed'),'FontSize',20,'Fontweight','bold')
textm(59,-62,'(d)','Color',rgb('Black'),'FontSize',26,'Fontweight','bold')
scatterm(45.25,-48.25,500,'k','filled','v')
scatterm(45.25,-48.25,200,'w','filled','v')
scatterm(45.25,-48.25,50,'k','filled','v')
% contourm(lat_z,lon_z,z,...
%     [-100 -100],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-1000 -1000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-3000 -3000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-4000 -4000],'color',[77/256 77/256 77/256],'linewidth',0.5)
geoshow('GSHHS_i_L1.shp','FaceColor',rgb('Black'))
%geoshow('landareas.shp','FaceColor','k')
%plabel('PLabelLocation',10,'FontSize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south','FontSize',18)
%title('Probability map (Labrador Current, Forward 1 year, 2003-2012)')
set(p8,'color','none')
framem
tightmap

saveas(gcf,[dir_output,'manuscript_v3/TGB2_Figure5_v3'],'png')

%% Figure 6

% Map of mean potential thickness with NWA Slope and NW Corner areas

%nlayer = '19';

dir_input = '/Users/afonso/Documents/Research/TGB2/Data/Eulerian/';
dir_output = '/Users/afonso/Documents/Research/TGB2/output/';

%ncdisp([dir_output,'hycom_atl_19_monthly.nc'])

% lat and lon (high resolution)
lat_aux = ncread([dir_output,'hycom_atl_19_monthly.nc'],'lat');
lon_aux = ncread([dir_output,'hycom_atl_19_monthly.nc'],'lon');
%lat_aux = lat_aux';
%lon_aux = lon_aux' - 360;
%[lon1,lat1] = meshgrid(lon_aux,lat_aux);

% lat and lon (low resolution)
lat = ncread([dir_output,'ATLg008_eulerian_layer19.nc'],'latitude');
lon = ncread([dir_output,'ATLg008_eulerian_layer19.nc'],'longitude');
lat = lat';
lon = lon' - 360;


% Bathymetry
%ncdisp([dir_output,'depth_ATLg008_10.nc'])
bathy = ncread([dir_output,'depth_ATLg008_10.nc'],'bathymetry');
bathy(bathy >= 10^10) = NaN;% remove values on land
bathy = bathy';

% Create variable for time
year = 1993:2017;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(300,1);
for yy = 1:25
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

% Eulerian fields

d1 = dir('/Users/afonso/Documents/Research/TGB2/Data/Eulerian');
d1 = d1(4:end);
thick = [];
depth = [];
for ii=1:length(d1)
    load([dir_input,d1(ii).name]);
    thick = cat(3,thick,hh(817:1033,263:727));
    depth = cat(3,depth,zu(817:1033,263:727));
    clear hh ss tt zu
end
clear lat_aux lon_aux ii d1

% Potential thickness (top of layer to bottom of ocean)
f = gsw_f(lat);
f_40 = gsw_f(40);
thick_pot = NaN(size(depth));
for tt = 1:length(time)
    thick_pot(:,:,tt) = f_40 * (bathy - squeeze(depth(:,:,tt))) ./ f;
end

thick_pot_anom = NaN(size(thick_pot));
for tt = 1:length(time)
    thick_pot_anom(:,:,tt) = squeeze(thick_pot(:,:,tt)) - squeeze(nanmean(thick_pot,3));
end


% % Potential thickness (top of layer to bottom of layer)
% f = gsw_f(lat);
% f_40 = gsw_f(40);
% thick_pot_layer = NaN(size(depth));
% for tt = 1:length(time)
%     thick_pot_layer(tt,:,:) = f_40 * squeeze(thick(tt,:,:)) ./ f';
% end

%file_field = [dir_output,'hycom_glo_surface1_particles1-1.nc'];
%file_traj = [dir_output,'ATLg008_layer19_LC_80'];
file_traj = [dir_output,'ATLg_008_layer19_LC80'];


ncdisp([file_traj,'.nc']);

pset.trajectory = ncread([file_traj,'.nc'],'trajectory');
pset.lon = ncread([file_traj,'.nc'],'lon');
pset.lat = ncread([file_traj,'.nc'],'lat');
pset.time = ncread([file_traj,'.nc'],'time');
pset.time2 = pset.time/(60*60*24)+datenum(1993,1,1,12,0,0);

load([dir_output,'etopo1_TGB.mat'],'lat_z','lon_z','z')

% % Create variable for time
% year = 1993:2017;
% month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
%     'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
% kk = 0;
% time = NaN(300,1);
% for yy = 1:25
%     for mm = 1:12
%         kk = kk + 1;
%         if mm < 10
%             time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
%                 'yyyymmdd');
%         else
%             time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
%                 'yyyymmdd');
%         end
%     end
% end
% clear kk yy mm

figure('Position', [100, 200, 800, 600])
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
hold on
c1_aux = contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
linem([36 36],[360-55 360-52],'k')
linem([46 46],[360-55 360-52],'k')
linem([36 46],[360-55 360-55],'k')
linem([36 46],[360-52 360-52],'k')
geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])

close all

c1 = NaN(2,1);
j = 1;
for i = 2:length(c1_aux)
    if c1_aux(1,i)>=360-55 & c1_aux(1,i)<=360-52 & c1_aux(2,i)>=36 & c1_aux(2,i)<=46
        c1(:,j) = c1_aux(:,i);
        if j>1
            if abs(c1_aux(1,i)-c1(1,j-1))<0.2
                j = j+1;
            else
                c1(:,j) = NaN;
            end
        else
            j = j+1;
        end
    end
end
clear i j

if isnan(c1(1,end))
    c1 = c1(:,1:end-1);% lon,lat of 1000m isobath
end

c1(3,:) = zeros;
c1_dist = gsw_distance(c1(1,:),c1(2,:));
for i = 2:length(c1)
    c1(3,i) = sum(c1_dist(1:i-1));
end
clear c1_aux c1_dist i

% Line 1 degree to the south of the 1,000-m isobath

c2 = c1;
c2(2,:) = c1(2,:) - 1;
c3 = c1;
c3(2,:) = c1(2,:) - 2;

% figure('Position', [100, 200, 800, 600])
% axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
% hold on
% contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem(c1(2,:),c1(1,:),'k','linewidth',2)
% linem(c2(2,:),c2(1,:),'--k','linewidth',1)
% linem(c3(2,:),c3(1,:),'k','linewidth',2)
% %contourm(c1(2,:),c1(2,:),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem([36 36],[360-75 360-52],'k')
% linem([46 46],[360-75 360-52],'k')
% linem([36 46],[360-75 360-75],'k')
% linem([36 46],[360-52 360-52],'k')
% geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
% 
% close all

close all

% Define run time
time_run = 1;%years
time_total = 25;%years
%Subset
time_length = floor(size(pset.time,1)*(time_run/time_total));
traj_length = floor(size(pset.time,2)*((time_total-time_run)/time_total));

% How many particles get to the Slope Region in 90 days
traj_slope = NaN(traj_length,1);
for ii = 1:traj_length
    traj_slope(ii,1) = sum(inpolygon(pset.lon(1:90,ii),pset.lat(1:90,ii),...
        [c1(1,:)-360 fliplr(c3(1,:)-360) c1(1,1)-360],...
        [c1(2,:) fliplr(c3(2,:)) c1(2,1)]));
end

% Number of particles per release time (time series)
slope_releasetime = pset.time2(1,traj_slope>=1);
slope_releasetime_yr = str2num(datestr(slope_releasetime,'yyyy'));
slope_releasetime_mo = str2num(datestr(slope_releasetime,'mm'));
perc1_slope = 100*length(slope_releasetime)/traj_length; % 15% of floats make it to slope

yr = 1993:2016;
slope_releasetime_yr1 = NaN(length(yr),1);
for ii = 1:length(yr)
    slope_releasetime_yr1(ii,1) = sum(slope_releasetime_yr==yr(ii));
end
perc_yr_slope = 100*slope_releasetime_yr1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(yr)
    perc_yr1_slope(ii,1) = 100*slope_releasetime_yr1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii));
end

mo = 1:12;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

slope_releasetime_mo1 = NaN(length(mo),1);
for ii = 1:length(mo)
    slope_releasetime_mo1(ii,1) = sum(slope_releasetime_mo==mo(ii));
end
perc_mo_slope = 100*slope_releasetime_mo1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(mo)
    perc_mo1_slope(ii,1) = 100*slope_releasetime_mo1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(ii));
end

% Area along the NWA Slope
lat_nwa = pset.lat(1:time_length,traj_slope>=1);
lon_nwa = pset.lon(1:time_length,traj_slope>=1);
time2_nwa = pset.time2(1:time_length,traj_slope>=1);

% Releases per month
slope_release_number = NaN(288,5);
for ii = 1:length(yr)
    for jj = 1:length(mo)
        slope_release_number(12*(ii-1)+jj,1) = yr(ii);
        slope_release_number(12*(ii-1)+jj,2) = mo(jj);
        slope_release_number(12*(ii-1)+jj,3) = ...
            sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii)...
            & str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(jj));
        slope_release_number(12*(ii-1)+jj,4) = sum(str2num(datestr(time2_nwa(1,:),'mm'))...
            ==mo(jj) & str2num(datestr(time2_nwa(1,:),'yyyy'))==yr(ii));
        slope_release_number(12*(ii-1)+jj,5) = 100 * slope_release_number(12*(ii-1)+jj,4)...
            /slope_release_number(12*(ii-1)+jj,3);
    end
end

% thick_pot_mean = squeeze(nanmean(thick_pot,1));
% interp2(lon(1,:),lat(:,1)',thick_pot_mean',-48.5,45.5)
% interp2(lon(1,:),lat(:,1)',thick_pot_mean',-48,45)
% 
% thick_pot_offshore = NaN(length(time),1);
% for tt = 1:length(time)
%     thick_pot_offshore(tt,1) = interp2(lon(1,:),lat(:,1)',...
%         squeeze(thick_pot(tt,:,:))',-48,45);
% end
% nanmean(thick_pot_offshore)
% nanstd(thick_pot_offshore)


% Load trajectories on mean velocity field
file_traj = [dir_output,'hycom_atl_19_for_LC_mean'];

ncdisp([file_traj,'.nc'])

pset.trajectory = ncread([file_traj,'.nc'],'trajectory');
pset.lon = ncread([file_traj,'.nc'],'lon');
pset.lat = ncread([file_traj,'.nc'],'lat');
pset.time = ncread([file_traj,'.nc'],'time');
pset.time2 = pset.time/(60*60*24)+datenum(1993,1,1,12,0,0);

% colormap

colormap_aux = readtable([dir_output,'4w_ROTB.csv'], 'HeaderLines',1);
colormap_wave2 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'3-wave-yellow-grey-blue.csv'], 'HeaderLines',1);
colormap_wave3 = table2array(colormap_aux(:,3:5));
clear colormap_aux

% slope of the potential thickness contours
thick_pot_mean = squeeze(nanmean(thick_pot,3));

% plot
figure('Position',[100,200,720,600])
%figure('Units','Inches','Position', [0, 0, 13.33, 7.5])
p1 = subplot('Position',[0.08 0.55 0.82 0.45]);
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[34 47],'MapLonLimit',[-77 -40])
hold on
pcolorm(lat,lon,squeeze(nanmean(thick_pot,3)))
shading interp
%colormap(p1,flipud(colormap_wave3(1:43,:)))
colormap(p1,flipud(colormap_wave2(1:end-9,:)))
c = colorbar('Orientation','Horizontal','Position',[0.11 0.84 0.2 0.05],...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'Pot Thickness (m)';
set(c,'Fontsize',18,'Fontweight','bold')
alpha 0.8
c1_aux = contourm(lat,lon,squeeze(nanmean(thick_pot,3)),[2300 2300],...
    'color',rgb('Black'),'linewidth',1.5);
c2_aux = contourm(lat,lon,squeeze(nanmean(thick_pot,3)),[3000 3000],...
    'color',rgb('Black'),'linewidth',1.5);
% Upstream
plotm([c1_aux(2,658),c2_aux(2,860)],[c1_aux(1,658),c2_aux(1,860)],'Color',rgb('DeepSkyBlue'),'LineWidth',4)
%scatterm(c1_aux(2,658),c1_aux(1,658),50,'r','filled')
%scatterm(c2_aux(2,860),c2_aux(1,860),50,'r','filled')
dist_upstream = gsw_distance([c1_aux(1,658),c2_aux(1,860)],[c1_aux(2,658),c2_aux(2,860)],0);
% TGB
plotm([c1_aux(2,620),c2_aux(2,730)],[c1_aux(1,620),c2_aux(1,730)],'Color',rgb('FireBrick'),'LineWidth',4)
%scatterm(c1_aux(2,620),c1_aux(1,620),50,'r','filled')
%scatterm(c2_aux(2,730),c2_aux(1,730),50,'r','filled')
dist_tgb = gsw_distance([c1_aux(1,620),c2_aux(1,730)],[c1_aux(2,620),c2_aux(2,730)],0);
%plotm([49 51 51 49 49],[-43 -43 -40 -40 -43],'k','linewidth',3) 
%plotm([44 46 46 44 44],[-44 -44 -41 -41 -44],'k','linewidth',3) 
plotm([c1(2,:) fliplr(c3(2,:)) c1(2,1)],[c1(1,:) fliplr(c3(1,:)) c1(1,1)],'k','linewidth',6)
%linem([45.5 45],[-48.5 -48],'k','LineWidth',12)
% for ii = 1:80
%     linem(pset.lat(:,ii),pset.lon(:,ii),...
%         'color',rgb('SteelBlue'),'linewidth',0.5)
% end
scatterm(45.25,-48.25,600,'k','filled','v')
scatterm(45.25,-48.25,240,'w','filled','v')
scatterm(45.25,-48.25,60,'k','filled','v')
textm(40.5,-76.8,'(a)','Color',rgb('Snow'),'FontSize',26,'Fontweight','bold')
f1 = fillm([41,43,43,41],[-53,-53,-50,-50],'w');
f1.FaceColor = 'w';
f1.EdgeColor = rgb('White');
f1.FaceAlpha = 0;
f1.EdgeAlpha = 1;
f1.LineWidth = 3;
% contourm(lat_z,lon_z,z,...
%     [-100 -100],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-1000 -1000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-3000 -3000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-4000 -4000],'color',[77/256 77/256 77/256],'linewidth',0.5)
geoshow('GSHHS_i_L1.shp','FaceColor',rgb('Black'))
%geoshow('landareas.shp','FaceColor','k')
plabel('PLabelLocation',5,'FontSize',18)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%title('Probability map (Labrador Current, Forward 1 year, 2003-2012)')
%set(p2,'color','none')
framem
tightmap

%figure('Position',[100,200,640,760])
p2 = subplot('Position',[0.08 0.09 0.82 0.37]);
hold on
colororder([rgb('SteelBlue');rgb('Black')])
yyaxis left
b = bar(time(1:end-12),slope_release_number(:,5));
b.EdgeColor = rgb('SteelBlue');
%plot(time(1:end-12),slope_release_number(:,5),'Color',rgb('MidnightBlue'),'LineWidth',1.5)
text(time(2),56,'(b)','Color',rgb('Black'),'FontSize',26,'Fontweight','bold')
ylabel('% of Slope Particles')
ylim([0 60])
set(gca,'FontSize',18)
datetick('x','yy')
xlabel('Time (yr)')
xlim([time(1) time(end)-365])
set(gca,'FontSize',18,'xtick',...
    [datenum('19940101','yyyymmdd'):730:datenum('20160101','yyyymmdd')],...
    'xticklabel',num2str(year(2:2:end)'))
box on
grid on

yyaxis right
thick_pot_anom_tgb = squeeze(nanmean(nanmean(thick_pot_anom...
    (lat(:,1)>=41 & lat(:,1) <=43,lon(1,:)>=-53 & lon(1,:)<=-50,:),1),2));
%thick_pot_anom_tgb2 = nanmean(nanmean(thick_pot_anom...
%    (lat(:,1)>=41 & lat(:,1) <=43,lon(1,:)>=-49 & lon(1,:)<=-47,:),1),2);
%plot(time(121:228),thick_pot_anom_tgb(1:108),'color',rgb('Black'),'LineWidth',2.5)
plot(time,thick_pot_anom_tgb,'color',rgb('Black'),'LineWidth',1.8)
ylabel('Pot Thick anom (m)')
xlabel('Time (yr)')
set(gca,'FontSize',18)
ylim([-150 150])
xlim([time(1) time(end)-365])
set(gca,'FontSize',18,'xtick',...
    [datenum('19940101','yyyymmdd'):730:datenum('20160101','yyyymmdd')],...
    'xticklabel',num2str(year(2:2:end)'))
grid on
saveas(gcf,[dir_output,'manuscript_v3/TGB2_Figure6_v4'],'png')

% Statistics
[r0 p0] = corrcoef(slope_release_number(:,5),thick_pot_anom_tgb(1:end-12))
[r1 p1] = corrcoef(slope_release_number(:,5),thick_pot_anom_tgb(2:end-11))
[r2 p2] = corrcoef(slope_release_number(:,5),thick_pot_anom_tgb(3:end-10))

% Change point
[ipt_slope, residuals_slope] = findchangepts(slope_release_number(:,5)) %158 (Feb/2006)
[ipt_thick, residuals_thick] = findchangepts(thick_pot_anom_tgb) %169 (Jan/2007)

figure('Position',[100,200,720,320])
%figure('Units','Inches','Position', [0, 0, 13.33, 7.5])
p1 = subplot('Position',[0.1 0.1 0.86 0.86]);
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[40 47],'MapLonLimit',[-60 -40])
hold on
pcolorm(lat,lon,squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(p1,flipud(colormap_wave3(1:43,:)))
colormap(p1,flipud(colormap_wave2(1:end-9,:)))
%c = colorbar('Orientation','Horizontal','Position',[0.11 0.84 0.2 0.05],...
%    'Color',[1-eps,1,1]);%colorbar over land
%c.Label.String = 'Pot Thickness (m)';
%set(c,'Fontsize',18,'Fontweight','bold')
alpha 0.3
%contourm(lat,lon,squeeze(nanmean(thick_pot,1))',[2300 2300],...
%    'color',rgb('Black'),'linewidth',1.5)
%contourm(lat,lon,squeeze(nanmean(thick_pot,1))',[3000 3000],...
%    'color',rgb('Black'),'linewidth',1.5)
%plotm([49 51 51 49 49],[-43 -43 -40 -40 -43],'k','linewidth',3) 
%plotm([44 46 46 44 44],[-44 -44 -41 -41 -44],'k','linewidth',3) 
%plotm([c1(2,:) fliplr(c3(2,:)) c1(2,1)],[c1(1,:) fliplr(c3(1,:)) c1(1,1)],'k','linewidth',6)
%linem([45.5 45],[-48.5 -48],'k','LineWidth',12)
for ii = 1:80
    linem(pset.lat(:,ii),pset.lon(:,ii),...
        'color',rgb('Black'),'linewidth',0.5)
end
scatterm(45.25,-48.25,600,'k','filled','v')
scatterm(45.25,-48.25,240,'w','filled','v')
scatterm(45.25,-48.25,60,'k','filled','v')
textm(40.5,-76.8,'(a)','Color',rgb('Snow'),'FontSize',26,'Fontweight','bold')
% f1 = fillm([41,43,43,41],[-53,-53,-50,-50],'w');
% f1.FaceColor = 'w';
% f1.EdgeColor = 'w';
% f1.FaceAlpha = 0;
% f1.EdgeAlpha = 1;
% f1.LineWidth = 5;
% contourm(lat_z,lon_z,z,...
%     [-100 -100],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-1000 -1000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-3000 -3000],'color',[77/256 77/256 77/256],'linewidth',0.5)
% contourm(lat_z,lon_z,z,...
%     [-4000 -4000],'color',[77/256 77/256 77/256],'linewidth',0.5)
%geoshow('GSHHS_i_L1.shp','FaceColor',rgb('Black'))
geoshow('landareas.shp','FaceColor','k')
plabel('PLabelLocation',5,'FontSize',18)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%title('Probability map (Labrador Current, Forward 1 year, 2003-2012)')
%set(p2,'color','none')
framem
tightmap
saveas(gcf,[dir_output,'manuscript_v2/TGB2_Figure6_aux'],'png')

% Figure mean traj

% Set  directory to save data
dir_input = '/Users/afonso/Documents/Research/TGB2/output/';
% Set directory to save figures and tables
dir_output = '/Users/afonso/Documents/Research/TGB2/output/';

% Create variable for time
year = 1993:2017;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(300,1);
for yy = 1:25
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

file_traj = [dir_output,'ATLg_008_layer19_LC80'];

ncdisp([file_traj,'.nc'])

pset.trajectory = ncread([file_traj,'.nc'],'trajectory');
pset.lon = ncread([file_traj,'.nc'],'lon');
pset.lat = ncread([file_traj,'.nc'],'lat');
pset.time = ncread([file_traj,'.nc'],'time');
pset.time2 = pset.time/(60*60*24)+datenum(1993,1,1,12,0,0);

load([dir_output,'etopo1_TGB.mat'],'lat_z','lon_z','z')

% Create variable for time
year = 1993:2017;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(300,1);
for yy = 1:25
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

figure('Position', [100, 200, 800, 600])
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
hold on
c1_aux = contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
linem([36 36],[360-55 360-52],'k')
linem([46 46],[360-55 360-52],'k')
linem([36 46],[360-55 360-55],'k')
linem([36 46],[360-52 360-52],'k')
geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])

close all

c1 = NaN(2,1);
j = 1;
for i = 2:length(c1_aux)
    if c1_aux(1,i)>=360-55 & c1_aux(1,i)<=360-52 & c1_aux(2,i)>=36 & c1_aux(2,i)<=46
        c1(:,j) = c1_aux(:,i);
        if j>1
            if abs(c1_aux(1,i)-c1(1,j-1))<0.2
                j = j+1;
            else
                c1(:,j) = NaN;
            end
        else
            j = j+1;
        end
    end
end
clear i j

if isnan(c1(1,end))
    c1 = c1(:,1:end-1);% lon,lat of 1000m isobath
end

c1(3,:) = zeros;
c1_dist = gsw_distance(c1(1,:),c1(2,:));
for i = 2:length(c1)
    c1(3,i) = sum(c1_dist(1:i-1));
end
clear c1_aux c1_dist i

% Line 1 degree to the south of the 1,000-m isobath

c2 = c1;
c2(2,:) = c1(2,:) - 1;
c3 = c1;
c3(2,:) = c1(2,:) - 2;

% figure('Position', [100, 200, 800, 600])
% axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
% hold on
% contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem(c1(2,:),c1(1,:),'k','linewidth',2)
% linem(c2(2,:),c2(1,:),'--k','linewidth',1)
% linem(c3(2,:),c3(1,:),'k','linewidth',2)
% %contourm(c1(2,:),c1(2,:),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem([36 36],[360-75 360-52],'k')
% linem([46 46],[360-75 360-52],'k')
% linem([36 46],[360-75 360-75],'k')
% linem([36 46],[360-52 360-52],'k')
% geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
% 
% close all

close all

% Define run time
time_run = 1;%years
time_total = 25;%years
%Subset
time_length = floor(size(pset.time,1)*(time_run/time_total));
traj_length = floor(size(pset.time,2)*((time_total-time_run)/time_total));

traj_slope = NaN(traj_length,1);
for ii = 1:traj_length
    traj_slope(ii,1) = sum(inpolygon(pset.lon(1:90,ii),pset.lat(1:90,ii),...
        [c1(1,:)-360 fliplr(c3(1,:)-360) c1(1,1)-360],...
        [c1(2,:) fliplr(c3(2,:)) c1(2,1)]));
end
clear ii

traj_slope1 = [traj_slope;zeros(size(pset.time2,2) - length(traj_slope),1)];

% Number of particles per release time (time series)
slope_releasetime = pset.time2(1,traj_slope>=1);
slope_releasetime_yr = str2num(datestr(slope_releasetime,'yyyy'));
slope_releasetime_mo = str2num(datestr(slope_releasetime,'mm'));
perc1_slope = 100*length(slope_releasetime)/traj_length; % 15% of floats make it to slope


yr = 1993:2016;
slope_releasetime_yr1 = NaN(length(yr),1);
for ii = 1:length(yr)
    slope_releasetime_yr1(ii,1) = sum(slope_releasetime_yr==yr(ii));
end
perc_yr_slope = 100*slope_releasetime_yr1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(yr)
    perc_yr1_slope(ii,1) = 100*slope_releasetime_yr1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii));
end

mo = 1:12;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

slope_releasetime_mo1 = NaN(length(mo),1);
for ii = 1:length(mo)
    slope_releasetime_mo1(ii,1) = sum(slope_releasetime_mo==mo(ii));
end
perc_mo_slope = 100*slope_releasetime_mo1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(mo)
    perc_mo1_slope(ii,1) = 100*slope_releasetime_mo1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(ii));
end

% Area along the NWA Slope
lat_nwa = pset.lat(1:time_length,traj_slope>=1);
lon_nwa = pset.lon(1:time_length,traj_slope>=1);
time2_nwa = pset.time2(1:time_length,traj_slope>=1);

% Releases per month
slope_release_number = NaN(288,5);
for ii = 1:length(yr)
    for jj = 1:length(mo)
        slope_release_number(12*(ii-1)+jj,1) = yr(ii);
        slope_release_number(12*(ii-1)+jj,2) = mo(jj);
        slope_release_number(12*(ii-1)+jj,3) = ...
            sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii)...
            & str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(jj));
        slope_release_number(12*(ii-1)+jj,4) = sum(str2num(datestr(time2_nwa(1,:),'mm'))...
            ==mo(jj) & str2num(datestr(time2_nwa(1,:),'yyyy'))==yr(ii));
        slope_release_number(12*(ii-1)+jj,5) = 100 * slope_release_number(12*(ii-1)+jj,4)...
            /slope_release_number(12*(ii-1)+jj,3);
    end
end


% Load bathymetry from ETOPO1
x = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','x');
lon_z1 = x(x<=-25 & x>=-90.6);
%lon_z = lon_z(1:60:end);
y = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','y');
lat_z1 = y(y<=60 & y>=28);
%lat_z = lat_z(1:60:end);
z_aux = ncread('/Users/afonso/Google Drive/Data/ETOPO1_Bed_g_gmt4.grd','z');
z1 = z_aux(x<=-25 & x>=-90.6,y<=60 & y>=28); z1 = z1';
%z = z(1:60:end,1:60:end);
clear x y z_aux

cbar_bathy = [5/256,112/256,176/256;...
    54/256,144/256,192/256;...
    116/256,169/256,207/256;...
    166/256,189/256,219/256;...
    208/256,209/256,230/256;...
    236/256,231/256,242/256;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0;...
    0,0,0];

% Load mean trajectories
file_traj = [dir_input,'ATLg_008_layer19_LC80_mean'];
file_traj = [dir_input,'ATLg_008_layer19_LC80_open'];
file_traj = [dir_input,'ATLg_008_layer19_LC80_closed'];

ncdisp([file_traj,'.nc']);

pset.trajectory = ncread([file_traj,'.nc'],'trajectory');
pset.lon = ncread([file_traj,'.nc'],'lon');
pset.lat = ncread([file_traj,'.nc'],'lat');
pset.time = ncread([file_traj,'.nc'],'time');
pset.time2 = pset.time/(60*60*24)+datenum(1993,1,1,12,0,0);

kk_aux = [];
for kk = 1:80
    if sum(pset.lon(:,kk)<-55) == 0
        kk_aux = cat(1,kk_aux,kk)
    end
end
100*24/80% 30% of the particles in the Open Valve experiment followed the shelf break

%figure('Units','Inches','Position', [0, 0, 13.33, 7.5])
%p1 = subplot('Position',[0.0 0.0 1.0 1.0]);
%axesm('MapProjection','eqdcylin','MapLatLimit',[28 60],'MapLonLimit',[-90.6 -25])
figure('Position',[100,200,800,500])
p1 = subplot('Position',[0.06 0.06 0.9 0.92]);
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30 60],'MapLonLimit',[-80 -25])
hold on
%pcolorm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end));
pcolorm(lat_z1,lon_z1,z1);
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',14)
shading flat
colormap(p1,cbar_bathy)
%contourm(lat_z(1:5:end),lon_z(1:5:end),z(1:5:end,1:5:end),[0 0],'color','k','LineWidth',1.5)
%colormap(cmocean('ice',50))
%cb = colorbar;
caxis([-6000 6000])
alpha 0.4
plabel('PLabelLocation',10,'FontSize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%framem
tightmap

p2 = axes('Position',[0.06 0.06 0.9 0.92]);
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -25])
hold on
for kk = 1:80
    linem(pset.lat(:,kk),...
        pset.lon(:,kk),...
        'color',rgb('Black'),'LineWidth',2)
    scatterm(pset.lat(end,kk),...
        pset.lon(end,kk),...
        50,rgb('Black'),'filled')
end
%plotm([44 46 46 44 44],[-44 -44 -41 -41 -44],'color',rgb('FireBrick'),'linewidth',6) 
%plotm([c1(2,:) fliplr(c3(2,:)) c1(2,1)],[c1(1,:) fliplr(c3(1,:)) c1(1,1)],'color',rgb('MidnightBlue'),'linewidth',6)
geoshow('GSHHS_i_L1.shp','FaceColor',rgb('Black'))
%geoshow('landareas.shp','FaceColor','k')
plabel('PLabelLocation',10,'FontSize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%title('Probability map (Labrador Current, Forward 1 year, 2003-2012)')
set(p2,'color','none')
%framem
tightmap
set(gcf, 'Color', 'w')    
saveas(gcf,[dir_output,'manuscript_v3/TGB2_FigureNew'],'png')

%figure('Units','Inches','Position', [0, 0, 13.33, 7.5])
%p1 = subplot('Position',[0.0 0.0 1.0 1.0]);
%axesm('MapProjection','eqdcylin','MapLatLimit',[28 60],'MapLonLimit',[-90.6 -25])
figure('Position',[100,200,800,500])
p1 = subplot('Position',[0.06 0.06 0.9 0.92]);
axesm('MapProjection','eqdcylin',...
    'MapLatLimit',[30 60],'MapLonLimit',[-80 -25])
hold on
%pcolorm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end));
pcolorm(lat_z1,lon_z1,z1);
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',14)
shading flat
colormap(p1,cbar_bathy)
%contourm(lat_z(1:5:end),lon_z(1:5:end),z(1:5:end,1:5:end),[0 0],'color','k','LineWidth',1.5)
%colormap(cmocean('ice',50))
%cb = colorbar;
caxis([-6000 6000])
alpha 0.4
plabel('PLabelLocation',10,'FontSize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%framem
tightmap

p2 = axes('Position',[0.06 0.06 0.9 0.92]);
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -25])
hold on
kk_aux = 0;
for kk = 1:80
    if sum(pset.lon(:,kk)<-55) == 1
        kk_aux = kk_aux + 1;
    linem(pset.lat(:,kk),...
        pset.lon(:,kk),...
        'color',rgb('Black'),'LineWidth',2)
    scatterm(pset.lat(end,kk),...
        pset.lon(end,kk),...
        50,rgb('Black'),'filled')
    end
end
%plotm([44 46 46 44 44],[-44 -44 -41 -41 -44],'color',rgb('FireBrick'),'linewidth',6) 
%plotm([c1(2,:) fliplr(c3(2,:)) c1(2,1)],[c1(1,:) fliplr(c3(1,:)) c1(1,1)],'color',rgb('MidnightBlue'),'linewidth',6)
geoshow('GSHHS_i_L1.shp','FaceColor',rgb('Black'))
%geoshow('landareas.shp','FaceColor','k')
plabel('PLabelLocation',10,'FontSize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south')
%title('Probability map (Labrador Current, Forward 1 year, 2003-2012)')
set(p2,'color','none')
%framem
tightmap
set(gcf, 'Color', 'w')    
saveas(gcf,[dir_output,'manuscript_v3/TGB2_FigureNew2'],'png')

%% Figure 7

nlayer = '19';

dir_output = '/Users/afonso/Documents/Research/TGB2/output/';

%ncdisp([dir_output,'ATLg008_eulerian_layer',nlayer,'.nc'])

% lat and lon
lat = ncread([dir_output,'ATLg008_eulerian_layer',nlayer,'.nc'],'latitude');
lon = ncread([dir_output,'ATLg008_eulerian_layer',nlayer,'.nc'],'longitude');
lat = lat';
lon = lon' - 360;

% Bathymetry
%ncdisp([dir_output,'depth_ATLg008_10.nc'])
bathy = ncread([dir_output,'depth_ATLg008_10.nc'],'bathymetry');
bathy(bathy >= 10^10) = NaN;% remove values on land

% Create variable for time
year = 2003:2012;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(120,1);
for yy = 1:10
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

% Thickness
thick = ncread([dir_output,'ATLg008_eulerian_layer',nlayer,'.nc'],'thickness');
thick = thick/9806; % convert from pressure to metric
thick(thick >= 10^10) = NaN;% remove values on land
%thick(thick <= 1) = NaN;% remove values on shelf

% Depth
depth = ncread([dir_output,'ATLg008_eulerian_layer',nlayer,'.nc'],'depth');


% Potential thickness (top of layer to bottom of ocean)
f = gsw_f(lat);
f_40 = gsw_f(40);
thick_pot = NaN(size(depth));
for tt = 1:length(time)
    thick_pot(tt,:,:) = f_40 * (bathy - squeeze(depth(tt,:,:))) ./ f';
end

%file_field = [dir_output,'hycom_glo_surface1_particles1-1.nc'];
file_traj = [dir_output,'ATLg008_layer19_LC_80'];

ncdisp([file_traj,'.nc'])

pset.trajectory = ncread([file_traj,'.nc'],'trajectory');
pset.lon = ncread([file_traj,'.nc'],'lon');
pset.lat = ncread([file_traj,'.nc'],'lat');
pset.time = ncread([file_traj,'.nc'],'time');
pset.time2 = pset.time/(60*60*24)+datenum(2003,1,1,12,0,0);

load([dir_output,'etopo1_TGB.mat'],'lat_z','lon_z','z')

% Create variable for time
year = 2003:2012;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(120,1);
for yy = 1:10
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

figure('Position', [100, 200, 800, 600])
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
hold on
c1_aux = contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
linem([36 36],[360-55 360-52],'k')
linem([46 46],[360-55 360-52],'k')
linem([36 46],[360-55 360-55],'k')
linem([36 46],[360-52 360-52],'k')
geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])

close all

c1 = NaN(2,1);
j = 1;
for i = 2:length(c1_aux)
    if c1_aux(1,i)>=360-55 & c1_aux(1,i)<=360-52 & c1_aux(2,i)>=36 & c1_aux(2,i)<=46
        c1(:,j) = c1_aux(:,i);
        if j>1
            if abs(c1_aux(1,i)-c1(1,j-1))<0.2
                j = j+1;
            else
                c1(:,j) = NaN;
            end
        else
            j = j+1;
        end
    end
end
clear i j

if isnan(c1(1,end))
    c1 = c1(:,1:end-1);% lon,lat of 1000m isobath
end

c1(3,:) = zeros;
c1_dist = gsw_distance(c1(1,:),c1(2,:));
for i = 2:length(c1)
    c1(3,i) = sum(c1_dist(1:i-1));
end
clear c1_aux c1_dist i

% Line 1 degree to the south of the 1,000-m isobath

c2 = c1;
c2(2,:) = c1(2,:) - 1;
c3 = c1;
c3(2,:) = c1(2,:) - 2;

% figure('Position', [100, 200, 800, 600])
% axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
% hold on
% contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem(c1(2,:),c1(1,:),'k','linewidth',2)
% linem(c2(2,:),c2(1,:),'--k','linewidth',1)
% linem(c3(2,:),c3(1,:),'k','linewidth',2)
% %contourm(c1(2,:),c1(2,:),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem([36 36],[360-75 360-52],'k')
% linem([46 46],[360-75 360-52],'k')
% linem([36 46],[360-75 360-75],'k')
% linem([36 46],[360-52 360-52],'k')
% geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
% 
% close all

close all

% Define run time
time_run = 1;

%Subset
time_length = floor(size(pset.time,1)*(time_run/10));
traj_length = floor(size(pset.time,2)*((10-time_run)/10));

traj_slope = NaN(traj_length,1);
for ii = 1:traj_length
    traj_slope(ii,1) = sum(inpolygon(pset.lon(1:10,ii),pset.lat(1:10,ii),...
        [c1(1,:)-360 fliplr(c3(1,:)-360) c1(1,1)-360],...
        [c1(2,:) fliplr(c3(2,:)) c1(2,1)]));
end
clear ii

traj_slope1 = [traj_slope;zeros(size(pset.time2,2) - length(traj_slope),1)];

% Number of particles per release time (time series)
slope_releasetime = pset.time2(1,traj_slope>=1);
slope_releasetime_yr = str2num(datestr(slope_releasetime,'yyyy'));
slope_releasetime_mo = str2num(datestr(slope_releasetime,'mm'));
perc1_slope = 100*length(slope_releasetime)/traj_length; % 15% of floats make it to slope


yr = 2003:2011;
slope_releasetime_yr1 = NaN(length(yr),1);
for ii = 1:length(yr)
    slope_releasetime_yr1(ii,1) = sum(slope_releasetime_yr==yr(ii));
end
perc_yr_slope = 100*slope_releasetime_yr1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(yr)
    perc_yr1_slope(ii,1) = 100*slope_releasetime_yr1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii));
end

mo = 1:12;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

slope_releasetime_mo1 = NaN(length(mo),1);
for ii = 1:length(mo)
    slope_releasetime_mo1(ii,1) = sum(slope_releasetime_mo==mo(ii));
end
perc_mo_slope = 100*slope_releasetime_mo1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(mo)
    perc_mo1_slope(ii,1) = 100*slope_releasetime_mo1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(ii));
end

% Area along the NWA Slope
lat_nwa = pset.lat(1:time_length,traj_slope>=1);
lon_nwa = pset.lon(1:time_length,traj_slope>=1);
time2_nwa = pset.time2(1:time_length,traj_slope>=1);

% Releases per month
slope_release_number = NaN(108,5);
for ii = 1:length(yr)
    for jj = 1:length(mo)
        slope_release_number(12*(ii-1)+jj,1) = yr(ii);
        slope_release_number(12*(ii-1)+jj,2) = mo(jj);
        slope_release_number(12*(ii-1)+jj,3) = ...
            sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii)...
            & str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(jj));
        slope_release_number(12*(ii-1)+jj,4) = sum(str2num(datestr(time2_nwa(1,:),'mm'))...
            ==mo(jj) & str2num(datestr(time2_nwa(1,:),'yyyy'))==yr(ii));
        slope_release_number(12*(ii-1)+jj,5) = 100 * slope_release_number(12*(ii-1)+jj,4)...
            /slope_release_number(12*(ii-1)+jj,3);
    end
end

slope_timestep = find(slope_release_number(:,5)==0);

% % split in two matrices to plot quiver in different colors
% u_low = NaN(size(u));
% u_high = NaN(size(u));
% v_low = NaN(size(v));
% v_high = NaN(size(v));
% 
% u_low(uv<=12) = u(uv<=12);
% u_high(uv>12) = u(uv>12);
% v_low(uv<=12) = v(uv<=12);
% v_high(uv>12) = v(uv>12);

plot_location = [0.01 0.71 0.32 0.29;...
    0.34 0.71 0.32 0.29;...
    0.67 0.71 0.32 0.29;...
    0.01 0.42 0.32 0.29;...
    0.34 0.42 0.32 0.29;...
    0.67 0.42 0.32 0.29;...
    0.01 0.13 0.32 0.29;...
    0.34 0.13 0.32 0.29;...
    0.67 0.13 0.32 0.29];

colormap_aux = readtable([dir_output,'5w_BRgpb.csv'], 'HeaderLines',1);
colormap_wave = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'4w_ROTB.csv'], 'HeaderLines',1);
colormap_wave2 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'3-wave-yellow-grey-blue.csv'], 'HeaderLines',1);
colormap_wave3 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'5-step-melow-wave.csv'], 'HeaderLines',1);
colormap_wave4 = table2array(colormap_aux(:,3:5));
clear colormap_aux 

colormap_aux = readtable([dir_output,'turqoise-olive.csv'], 'HeaderLines',1);
colormap_wave5 = table2array(colormap_aux(:,3:5));
clear colormap_aux 



%for ii = 11:length(slope_timestep)
ii = 23; % Dec 2007
figure('Position', [0, 0, 860, 600])
s1 = subplot('Position',plot_location(1,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[37 47],'MapLonLimit',[-59.9 -40])
hold on
pcolorm(lat,lon,squeeze(thick_pot(slope_timestep(ii),:,:))' - ...
    squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(cmocean('speed'))
%colormap(s1,flipud(colormap_wave3(1:43,:)))
colormap(s1,flipud(colormap_wave5))
cb = colorbar('Orientation','Horizontal','Position',[0.3,0.095,0.4 0.035],...
    'Color',rgb('Black'));%colorbar over land
cb.Label.String = 'Pontential Thickness Anomaly (m)';
set(cb,'Fontsize',22)%,'Fontweight','bold')
caxis([-250 250])
contourm(lat,lon,bathy',[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
contourm(lat,lon,bathy',[4000 4000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
alpha 0.9
% all particles started in ii month
a = find(str2num(datestr(pset.time2(1,:),'yyyy')) == ...
    slope_release_number(slope_timestep(ii),1) & ...
    str2num(datestr(pset.time2(1,:),'mm')) == ...
    slope_release_number(slope_timestep(ii),2));
% Slope trajectories started in ii month
b = find(str2num(datestr(pset.time2(1,:),'yyyy')) == ...
    slope_release_number(slope_timestep(ii),1) & ...
    str2num(datestr(pset.time2(1,:),'mm')) == ...
    slope_release_number(slope_timestep(ii),2) & ...
    traj_slope1 >= 1);
% All but slope trajectories started in ii month
c = setdiff(a,b);
% How many particles?
c_subsample = randi([c(1) c(end)],15,1);
%     for jj = 1:length(c)
%         linem(pset.lat(1:10,c(jj)),...
%             pset.lon(1:10,c(jj)),...
%             'color',rgb('Salmon'),'LineWidth',0.3)
%     end
%     for kk = 1:length(b)
%         linem(pset.lat(1:10,b(kk)),...
%             pset.lon(1:10,b(kk)),...
%             'color',rgb('SkyBlue'),'LineWidth',0.3)
%     end
for jj = 1:length(c_subsample)
    month_this = find(str2num(datestr(pset.time2(:,c_subsample(jj)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii),1) & ...
        str2num(datestr(pset.time2(:,c_subsample(jj)),'mm')) == ...
        slope_release_number(slope_timestep(ii),2));
    linem(pset.lat(month_this,c_subsample(jj)),...
        pset.lon(month_this,c_subsample(jj)),...
        'color',rgb('DarkRed'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),c(jj)),...
    %             pset.lon(month_this(1),c(jj)),...
    %             30,rgb('DarkRed'))
    scatterm(pset.lat(month_this(end),c_subsample(jj)),...
        pset.lon(month_this(end),c_subsample(jj)),...
        30,rgb('DarkRed'),'filled')
    clear month_this
end
%clear a b c jj kk
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[450 450],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[700 700],...
%    'Color','k','LineStyle','-','LineWidth',2)
textm(46.3,-59.8,'(a)','Color',rgb('Black'),'FontSize',26)%,'Fontweight','bold')
%textm(42.5,-74.8,'19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
%geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
textm(46.3,-57.6,[datestr(time(slope_timestep(ii)),'mmm'),' ',...
    datestr(time(slope_timestep(ii)),'yyyy')],...
    'Color',rgb('Black'),'FontSize',26)
%     title([datestr(time(slope_timestep(ii)),'mmm'),' ',...
%         datestr(time(slope_timestep(ii)),'yyyy'),', ',...
%         num2str(round(slope_release_number(slope_timestep(ii),5),1)),...
%         '% to Western GB'],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

s2 = subplot('Position',plot_location(2,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[37 47],'MapLonLimit',[-59.9 -40])
hold on
pcolorm(lat,lon,squeeze(thick_pot(slope_timestep(ii)+1,:,:))' - ...
    squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(cmocean('speed'))
colormap(s2,flipud(colormap_wave5))
%     cb = colorbar('Orientation','Horizontal','Position',cbar_location(2,:),...
%         'Color',[1-eps,1,1]);%colorbar over land
%     cb.Label.String = 'H Pot anom (m)';
%     set(cb,'Fontsize',14)%,'Fontweight','bold')
caxis([-250 250])
contourm(lat,lon,bathy',[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
contourm(lat,lon,bathy',[4000 4000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
alpha 0.9

%     for jj = 1:length(c)
%         linem(pset.lat(1:10,c(jj)),...
%             pset.lon(1:10,c(jj)),...
%             'color',rgb('Salmon'),'LineWidth',0.3)
%     end
%     for kk = 1:length(b)
%         linem(pset.lat(1:10,b(kk)),...
%             pset.lon(1:10,b(kk)),...
%             'color',rgb('SkyBlue'),'LineWidth',0.3)
%     end
for jj = 1:length(c_subsample)
    month_this = find(str2num(datestr(pset.time2(:,c_subsample(jj)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii)+1,1) & ...
        str2num(datestr(pset.time2(:,c_subsample(jj)),'mm')) == ...
        slope_release_number(slope_timestep(ii)+1,2));
    linem(pset.lat(1:month_this(1)-1,c_subsample(jj)),...
        pset.lon(1:month_this(1)-1,c_subsample(jj)),...
        'color',rgb('Salmon'),'LineWidth',0.8)
    linem(pset.lat(month_this(1)-1:month_this(end),c_subsample(jj)),...
        pset.lon(month_this(1)-1:month_this(end),c_subsample(jj)),...
        'color',rgb('DarkRed'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),c(jj)),...
    %             pset.lon(month_this(1),c(jj)),...
    %             30,rgb('DarkRed'))
    scatterm(pset.lat(month_this(end),c_subsample(jj)),...
        pset.lon(month_this(end),c_subsample(jj)),...
        30,rgb('DarkRed'),'filled')
    clear month_this
end

%clear a b c jj kk
%contourm(lat,lon,squeeze(nanmean(depth(6*12-11:6*12,:,:),1))',[450 450],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(depth(6*12-11:6*12,:,:),1))',[700 700],...
%    'Color','k','LineStyle','-','LineWidth',2)
textm(46.3,-59.8,'(b)','Color',rgb('Black'),'FontSize',26)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
%geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
textm(46.3,-57.6,[datestr(time(slope_timestep(ii)+1),'mmm'),' ',...
    datestr(time(slope_timestep(ii)+1),'yyyy')],...
    'Color',rgb('Black'),'FontSize',26)
%     title([datestr(time(slope_timestep(ii)+1),'mmm'),' ',...
%         datestr(time(slope_timestep(ii)+1),'yyyy')],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

s3 = subplot('Position',plot_location(3,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[37 47],'MapLonLimit',[-59.9 -40])
hold on
pcolorm(lat,lon,squeeze(thick_pot(slope_timestep(ii)+2,:,:))' - ...
    squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(cmocean('speed'))
colormap(s3,flipud(colormap_wave5))
%     cb = colorbar('Orientation','Horizontal','Position',cbar_location(3,:),...
%         'Color',[1-eps,1,1]);%colorbar over land
%     cb.Label.String = 'H Pot anom (m)';
%     set(cb,'Fontsize',14)%,'Fontweight','bold')
caxis([-250 250])
contourm(lat,lon,bathy',[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
contourm(lat,lon,bathy',[4000 4000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
alpha 0.9

%     for jj = 1:length(c)
%         linem(pset.lat(1:10,c(jj)),...
%             pset.lon(1:10,c(jj)),...
%             'color',rgb('Salmon'),'LineWidth',0.3)
%     end
%     for kk = 1:length(b)
%         linem(pset.lat(1:10,b(kk)),...
%             pset.lon(1:10,b(kk)),...
%             'color',rgb('SkyBlue'),'LineWidth',0.3)
%     end
for jj = 1:length(c_subsample)
    month_this = find(str2num(datestr(pset.time2(:,c_subsample(jj)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii)+2,1) & ...
        str2num(datestr(pset.time2(:,c_subsample(jj)),'mm')) == ...
        slope_release_number(slope_timestep(ii)+2,2));
    linem(pset.lat(1:month_this(1)-1,c_subsample(jj)),...
        pset.lon(1:month_this(1)-1,c_subsample(jj)),...
        'color',rgb('Salmon'),'LineWidth',0.8)
    linem(pset.lat(month_this(1)-1:month_this(end),c_subsample(jj)),...
        pset.lon(month_this(1)-1:month_this(end),c_subsample(jj)),...
        'color',rgb('DarkRed'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),c(jj)),...
    %             pset.lon(month_this(1),c(jj)),...
    %             30,rgb('DarkRed'))
    scatterm(pset.lat(month_this(end),c_subsample(jj)),...
        pset.lon(month_this(end),c_subsample(jj)),...
        30,rgb('DarkRed'),'filled')
    clear month_this
end
for kk = 1:length(b)
    month_this = find(str2num(datestr(pset.time2(:,b(kk)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii)+2,1) & ...
        str2num(datestr(pset.time2(:,b(kk)),'mm')) == ...
        slope_release_number(slope_timestep(ii)+2,2));
    linem(pset.lat(1:month_this(1)-1,b(kk)),...
        pset.lon(1:month_this(1)-1,b(kk)),...
        'color',rgb('SkyBlue'),'LineWidth',0.8)
    linem(pset.lat(month_this(1)-1:month_this(end),b(kk)),...
        pset.lon(month_this(1)-1:month_this(end),b(kk)),...
        'color',rgb('MidnightBlue'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),b(kk)),...
    %             pset.lon(month_this(1),b(kk)),...
    %             30,rgb('MidnightBlue'))
    %         scatterm(pset.lat(month_this(end),b(kk)),...
    %             pset.lon(month_this(end),b(kk)),...
    %             30,rgb('MidnightBlue'),'filled')
    clear month_this
end
clear a b c jj kk
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.05 35.05],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.1 35.1],...
%    'Color','k','LineStyle','-','LineWidth',2)
textm(46.3,-59.8,'(c)','Color',rgb('Black'),'FontSize',26)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
%geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
textm(46.3,-57.6,[datestr(time(slope_timestep(ii)+2),'mmm'),' ',...
    datestr(time(slope_timestep(ii)+2),'yyyy')],...
    'Color',rgb('Black'),'FontSize',26)
%     title([datestr(time(slope_timestep(ii)+2),'mmm'),' ',...
%         datestr(time(slope_timestep(ii)+2),'yyyy')],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

ii = 35;% Feb 2009
s4 = subplot('Position',plot_location(4,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[37 47],'MapLonLimit',[-59.9 -40])
hold on
pcolorm(lat,lon,squeeze(thick_pot(slope_timestep(ii),:,:))' - ...
    squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(cmocean('speed'))
colormap(s4,flipud(colormap_wave5))
% cb = colorbar('Orientation','Horizontal','Position',[0.3,0.08,0.4 0.04],...
%     'Color',rgb('Black'));%colorbar over land
% cb.Label.String = 'Pontential Thickness Anomaly (m)';
% set(cb,'Fontsize',22)%,'Fontweight','bold')
caxis([-250 250])
contourm(lat,lon,bathy',[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
contourm(lat,lon,bathy',[4000 4000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
alpha 0.9
% all particles started in ii month
a = find(str2num(datestr(pset.time2(1,:),'yyyy')) == ...
    slope_release_number(slope_timestep(ii),1) & ...
    str2num(datestr(pset.time2(1,:),'mm')) == ...
    slope_release_number(slope_timestep(ii),2));
% Slope trajectories started in ii month
b = find(str2num(datestr(pset.time2(1,:),'yyyy')) == ...
    slope_release_number(slope_timestep(ii),1) & ...
    str2num(datestr(pset.time2(1,:),'mm')) == ...
    slope_release_number(slope_timestep(ii),2) & ...
    traj_slope1 >= 1);
% All but slope trajectories started in ii month
c = setdiff(a,b);
% How many particles?
c_subsample = randi([c(1) c(end)],15,1);
%     for jj = 1:length(c)
%         linem(pset.lat(1:10,c(jj)),...
%             pset.lon(1:10,c(jj)),...
%             'color',rgb('Salmon'),'LineWidth',0.3)
%     end
%     for kk = 1:length(b)
%         linem(pset.lat(1:10,b(kk)),...
%             pset.lon(1:10,b(kk)),...
%             'color',rgb('SkyBlue'),'LineWidth',0.3)
%     end
for jj = 1:length(c_subsample)
    month_this = find(str2num(datestr(pset.time2(:,c_subsample(jj)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii),1) & ...
        str2num(datestr(pset.time2(:,c_subsample(jj)),'mm')) == ...
        slope_release_number(slope_timestep(ii),2));
    linem(pset.lat(month_this,c_subsample(jj)),...
        pset.lon(month_this,c_subsample(jj)),...
        'color',rgb('DarkRed'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),c(jj)),...
    %             pset.lon(month_this(1),c(jj)),...
    %             30,rgb('DarkRed'))
    scatterm(pset.lat(month_this(end),c_subsample(jj)),...
        pset.lon(month_this(end),c_subsample(jj)),...
        30,rgb('DarkRed'),'filled')
    clear month_this
end
%clear a b c jj kk
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[450 450],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[700 700],...
%    'Color','k','LineStyle','-','LineWidth',2)
textm(46.3,-59.8,'(d)','Color',rgb('Black'),'FontSize',26)%,'Fontweight','bold')
%textm(42.5,-74.8,'19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
%geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
textm(46.3,-57.6,[datestr(time(slope_timestep(ii)),'mmm'),' ',...
    datestr(time(slope_timestep(ii)),'yyyy')],...
    'Color',rgb('Black'),'FontSize',26)
%     title([datestr(time(slope_timestep(ii)),'mmm'),' ',...
%         datestr(time(slope_timestep(ii)),'yyyy'),', ',...
%         num2str(round(slope_release_number(slope_timestep(ii),5),1)),...
%         '% to Western GB'],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

s5 = subplot('Position',plot_location(5,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[37 47],'MapLonLimit',[-59.9 -40])
hold on
pcolorm(lat,lon,squeeze(thick_pot(slope_timestep(ii)+1,:,:))' - ...
    squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(cmocean('speed'))
colormap(s5,flipud(colormap_wave5))
%     cb = colorbar('Orientation','Horizontal','Position',cbar_location(2,:),...
%         'Color',[1-eps,1,1]);%colorbar over land
%     cb.Label.String = 'H Pot anom (m)';
%     set(cb,'Fontsize',14)%,'Fontweight','bold')
caxis([-250 250])
contourm(lat,lon,bathy',[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
contourm(lat,lon,bathy',[4000 4000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
alpha 0.9

%     for jj = 1:length(c)
%         linem(pset.lat(1:10,c(jj)),...
%             pset.lon(1:10,c(jj)),...
%             'color',rgb('Salmon'),'LineWidth',0.3)
%     end
%     for kk = 1:length(b)
%         linem(pset.lat(1:10,b(kk)),...
%             pset.lon(1:10,b(kk)),...
%             'color',rgb('SkyBlue'),'LineWidth',0.3)
%     end
for jj = 1:length(c_subsample)
    month_this = find(str2num(datestr(pset.time2(:,c_subsample(jj)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii)+1,1) & ...
        str2num(datestr(pset.time2(:,c_subsample(jj)),'mm')) == ...
        slope_release_number(slope_timestep(ii)+1,2));
    linem(pset.lat(1:month_this(1)-1,c_subsample(jj)),...
        pset.lon(1:month_this(1)-1,c_subsample(jj)),...
        'color',rgb('Salmon'),'LineWidth',0.8)
    linem(pset.lat(month_this(1)-1:month_this(end),c_subsample(jj)),...
        pset.lon(month_this(1)-1:month_this(end),c_subsample(jj)),...
        'color',rgb('DarkRed'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),c(jj)),...
    %             pset.lon(month_this(1),c(jj)),...
    %             30,rgb('DarkRed'))
    scatterm(pset.lat(month_this(end),c_subsample(jj)),...
        pset.lon(month_this(end),c_subsample(jj)),...
        30,rgb('DarkRed'),'filled')
    clear month_this
end

%clear a b c jj kk
%contourm(lat,lon,squeeze(nanmean(depth(6*12-11:6*12,:,:),1))',[450 450],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(depth(6*12-11:6*12,:,:),1))',[700 700],...
%    'Color','k','LineStyle','-','LineWidth',2)
textm(46.3,-59.8,'(e)','Color',rgb('Black'),'FontSize',26)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
%geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
textm(46.3,-57.6,[datestr(time(slope_timestep(ii)+1),'mmm'),' ',...
    datestr(time(slope_timestep(ii)+1),'yyyy')],...
    'Color',rgb('Black'),'FontSize',26)
%     title([datestr(time(slope_timestep(ii)+1),'mmm'),' ',...
%         datestr(time(slope_timestep(ii)+1),'yyyy')],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

s6 = subplot('Position',plot_location(6,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[37 47],'MapLonLimit',[-59.9 -40])
hold on
pcolorm(lat,lon,squeeze(thick_pot(slope_timestep(ii)+2,:,:))' - ...
    squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(cmocean('speed'))
colormap(s6,flipud(colormap_wave5))
%     cb = colorbar('Orientation','Horizontal','Position',cbar_location(3,:),...
%         'Color',[1-eps,1,1]);%colorbar over land
%     cb.Label.String = 'H Pot anom (m)';
%     set(cb,'Fontsize',14)%,'Fontweight','bold')
caxis([-250 250])
contourm(lat,lon,bathy',[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
contourm(lat,lon,bathy',[4000 4000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
alpha 0.9

%     for jj = 1:length(c)
%         linem(pset.lat(1:10,c(jj)),...
%             pset.lon(1:10,c(jj)),...
%             'color',rgb('Salmon'),'LineWidth',0.3)
%     end
%     for kk = 1:length(b)
%         linem(pset.lat(1:10,b(kk)),...
%             pset.lon(1:10,b(kk)),...
%             'color',rgb('SkyBlue'),'LineWidth',0.3)
%     end
for jj = 1:length(c_subsample)
    month_this = find(str2num(datestr(pset.time2(:,c_subsample(jj)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii)+2,1) & ...
        str2num(datestr(pset.time2(:,c_subsample(jj)),'mm')) == ...
        slope_release_number(slope_timestep(ii)+2,2));
    linem(pset.lat(1:month_this(1)-1,c_subsample(jj)),...
        pset.lon(1:month_this(1)-1,c_subsample(jj)),...
        'color',rgb('Salmon'),'LineWidth',0.8)
    linem(pset.lat(month_this(1)-1:month_this(end),c_subsample(jj)),...
        pset.lon(month_this(1)-1:month_this(end),c_subsample(jj)),...
        'color',rgb('DarkRed'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),c(jj)),...
    %             pset.lon(month_this(1),c(jj)),...
    %             30,rgb('DarkRed'))
    scatterm(pset.lat(month_this(end),c_subsample(jj)),...
        pset.lon(month_this(end),c_subsample(jj)),...
        30,rgb('DarkRed'),'filled')
    clear month_this
end
for kk = 1:length(b)
    month_this = find(str2num(datestr(pset.time2(:,b(kk)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii)+2,1) & ...
        str2num(datestr(pset.time2(:,b(kk)),'mm')) == ...
        slope_release_number(slope_timestep(ii)+2,2));
    linem(pset.lat(1:month_this(1)-1,b(kk)),...
        pset.lon(1:month_this(1)-1,b(kk)),...
        'color',rgb('SkyBlue'),'LineWidth',0.8)
    linem(pset.lat(month_this(1)-1:month_this(end),b(kk)),...
        pset.lon(month_this(1)-1:month_this(end),b(kk)),...
        'color',rgb('MidnightBlue'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),b(kk)),...
    %             pset.lon(month_this(1),b(kk)),...
    %             30,rgb('MidnightBlue'))
    %         scatterm(pset.lat(month_this(end),b(kk)),...
    %             pset.lon(month_this(end),b(kk)),...
    %             30,rgb('MidnightBlue'),'filled')
    clear month_this
end
clear a b c jj kk
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.05 35.05],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.1 35.1],...
%    'Color','k','LineStyle','-','LineWidth',2)
textm(46.3,-59.8,'(f)','Color',rgb('Black'),'FontSize',26)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
%geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
textm(46.3,-57.6,[datestr(time(slope_timestep(ii)+2),'mmm'),' ',...
    datestr(time(slope_timestep(ii)+2),'yyyy')],...
    'Color',rgb('Black'),'FontSize',26)
%     title([datestr(time(slope_timestep(ii)+2),'mmm'),' ',...
%         datestr(time(slope_timestep(ii)+2),'yyyy')],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

ii = 44;% Jan 2010
%figure('Position', [0, 0, 1340, 600])
s7 = subplot('Position',plot_location(7,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[37 47],'MapLonLimit',[-59.9 -40])
hold on
pcolorm(lat,lon,squeeze(thick_pot(slope_timestep(ii),:,:))' - ...
    squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(cmocean('speed'))
colormap(s7,flipud(colormap_wave5))
% cb = colorbar('Orientation','Horizontal','Position',[0.3,0.08,0.4 0.04],...
%     'Color',rgb('Black'));%colorbar over land
% cb.Label.String = 'Pontential Thickness Anomaly (m)';
% set(cb,'Fontsize',22)%,'Fontweight','bold')
caxis([-250 250])
contourm(lat,lon,bathy',[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
contourm(lat,lon,bathy',[4000 4000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
alpha 0.9
% all particles started in ii month
a = find(str2num(datestr(pset.time2(1,:),'yyyy')) == ...
    slope_release_number(slope_timestep(ii),1) & ...
    str2num(datestr(pset.time2(1,:),'mm')) == ...
    slope_release_number(slope_timestep(ii),2));
% Slope trajectories started in ii month
b = find(str2num(datestr(pset.time2(1,:),'yyyy')) == ...
    slope_release_number(slope_timestep(ii),1) & ...
    str2num(datestr(pset.time2(1,:),'mm')) == ...
    slope_release_number(slope_timestep(ii),2) & ...
    traj_slope1 >= 1);
% All but slope trajectories started in ii month
c = setdiff(a,b);
% How many particles?
c_subsample = randi([c(1) c(end)],15,1);
%     for jj = 1:length(c)
%         linem(pset.lat(1:10,c(jj)),...
%             pset.lon(1:10,c(jj)),...
%             'color',rgb('Salmon'),'LineWidth',0.3)
%     end
%     for kk = 1:length(b)
%         linem(pset.lat(1:10,b(kk)),...
%             pset.lon(1:10,b(kk)),...
%             'color',rgb('SkyBlue'),'LineWidth',0.3)
%     end
for jj = 1:length(c_subsample)
    month_this = find(str2num(datestr(pset.time2(:,c_subsample(jj)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii),1) & ...
        str2num(datestr(pset.time2(:,c_subsample(jj)),'mm')) == ...
        slope_release_number(slope_timestep(ii),2));
    linem(pset.lat(month_this,c_subsample(jj)),...
        pset.lon(month_this,c_subsample(jj)),...
        'color',rgb('DarkRed'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),c(jj)),...
    %             pset.lon(month_this(1),c(jj)),...
    %             30,rgb('DarkRed'))
    scatterm(pset.lat(month_this(end),c_subsample(jj)),...
        pset.lon(month_this(end),c_subsample(jj)),...
        30,rgb('DarkRed'),'filled')
    clear month_this
end
%clear a b c jj kk
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[450 450],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[700 700],...
%    'Color','k','LineStyle','-','LineWidth',2)
textm(46.3,-59.8,'(g)','Color',rgb('Black'),'FontSize',26)%,'Fontweight','bold')
%textm(42.5,-74.8,'19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
%geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
textm(46.3,-57.6,[datestr(time(slope_timestep(ii)),'mmm'),' ',...
    datestr(time(slope_timestep(ii)),'yyyy')],...
    'Color',rgb('Black'),'FontSize',26)
%     title([datestr(time(slope_timestep(ii)),'mmm'),' ',...
%         datestr(time(slope_timestep(ii)),'yyyy'),', ',...
%         num2str(round(slope_release_number(slope_timestep(ii),5),1)),...
%         '% to Western GB'],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

s8 = subplot('Position',plot_location(8,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[37 47],'MapLonLimit',[-59.9 -40])
hold on
pcolorm(lat,lon,squeeze(thick_pot(slope_timestep(ii)+1,:,:))' - ...
    squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(cmocean('speed'))
colormap(s8,flipud(colormap_wave5))
%     cb = colorbar('Orientation','Horizontal','Position',cbar_location(2,:),...
%         'Color',[1-eps,1,1]);%colorbar over land
%     cb.Label.String = 'H Pot anom (m)';
%     set(cb,'Fontsize',14)%,'Fontweight','bold')
caxis([-250 250])
contourm(lat,lon,bathy',[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
contourm(lat,lon,bathy',[4000 4000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
alpha 0.9

%     for jj = 1:length(c)
%         linem(pset.lat(1:10,c(jj)),...
%             pset.lon(1:10,c(jj)),...
%             'color',rgb('Salmon'),'LineWidth',0.3)
%     end
%     for kk = 1:length(b)
%         linem(pset.lat(1:10,b(kk)),...
%             pset.lon(1:10,b(kk)),...
%             'color',rgb('SkyBlue'),'LineWidth',0.3)
%     end
for jj = 1:length(c_subsample)
    month_this = find(str2num(datestr(pset.time2(:,c_subsample(jj)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii)+1,1) & ...
        str2num(datestr(pset.time2(:,c_subsample(jj)),'mm')) == ...
        slope_release_number(slope_timestep(ii)+1,2));
    linem(pset.lat(1:month_this(1)-1,c_subsample(jj)),...
        pset.lon(1:month_this(1)-1,c_subsample(jj)),...
        'color',rgb('Salmon'),'LineWidth',0.8)
    linem(pset.lat(month_this(1)-1:month_this(end),c_subsample(jj)),...
        pset.lon(month_this(1)-1:month_this(end),c_subsample(jj)),...
        'color',rgb('DarkRed'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),c(jj)),...
    %             pset.lon(month_this(1),c(jj)),...
    %             30,rgb('DarkRed'))
    scatterm(pset.lat(month_this(end),c_subsample(jj)),...
        pset.lon(month_this(end),c_subsample(jj)),...
        30,rgb('DarkRed'),'filled')
    clear month_this
end

%clear a b c jj kk
%contourm(lat,lon,squeeze(nanmean(depth(6*12-11:6*12,:,:),1))',[450 450],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(depth(6*12-11:6*12,:,:),1))',[700 700],...
%    'Color','k','LineStyle','-','LineWidth',2)
textm(46.3,-59.8,'(h)','Color',rgb('Black'),'FontSize',26)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
%geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
textm(46.3,-57.6,[datestr(time(slope_timestep(ii)+1),'mmm'),' ',...
    datestr(time(slope_timestep(ii)+1),'yyyy')],...
    'Color',rgb('Black'),'FontSize',26)
%     title([datestr(time(slope_timestep(ii)+1),'mmm'),' ',...
%         datestr(time(slope_timestep(ii)+1),'yyyy')],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

s9 = subplot('Position',plot_location(9,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[37 47],'MapLonLimit',[-59.9 -40])
hold on
pcolorm(lat,lon,squeeze(thick_pot(slope_timestep(ii)+2,:,:))' - ...
    squeeze(nanmean(thick_pot,1))')
shading interp
%colormap(cmocean('speed'))
colormap(s9,flipud(colormap_wave5))
%     cb = colorbar('Orientation','Horizontal','Position',cbar_location(3,:),...
%         'Color',[1-eps,1,1]);%colorbar over land
%     cb.Label.String = 'H Pot anom (m)';
%     set(cb,'Fontsize',14)%,'Fontweight','bold')
caxis([-250 250])
contourm(lat,lon,bathy',[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
contourm(lat,lon,bathy',[4000 4000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
alpha 0.9

%     for jj = 1:length(c)
%         linem(pset.lat(1:10,c(jj)),...
%             pset.lon(1:10,c(jj)),...
%             'color',rgb('Salmon'),'LineWidth',0.3)
%     end
%     for kk = 1:length(b)
%         linem(pset.lat(1:10,b(kk)),...
%             pset.lon(1:10,b(kk)),...
%             'color',rgb('SkyBlue'),'LineWidth',0.3)
%     end
for jj = 1:length(c_subsample)
    month_this = find(str2num(datestr(pset.time2(:,c_subsample(jj)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii)+2,1) & ...
        str2num(datestr(pset.time2(:,c_subsample(jj)),'mm')) == ...
        slope_release_number(slope_timestep(ii)+2,2));
    linem(pset.lat(1:month_this(1)-1,c_subsample(jj)),...
        pset.lon(1:month_this(1)-1,c_subsample(jj)),...
        'color',rgb('Salmon'),'LineWidth',0.8)
    linem(pset.lat(month_this(1)-1:month_this(end),c_subsample(jj)),...
        pset.lon(month_this(1)-1:month_this(end),c_subsample(jj)),...
        'color',rgb('DarkRed'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),c(jj)),...
    %             pset.lon(month_this(1),c(jj)),...
    %             30,rgb('DarkRed'))
    scatterm(pset.lat(month_this(end),c_subsample(jj)),...
        pset.lon(month_this(end),c_subsample(jj)),...
        30,rgb('DarkRed'),'filled')
    clear month_this
end
for kk = 1:length(b)
    month_this = find(str2num(datestr(pset.time2(:,b(kk)),'yyyy')) == ...
        slope_release_number(slope_timestep(ii)+2,1) & ...
        str2num(datestr(pset.time2(:,b(kk)),'mm')) == ...
        slope_release_number(slope_timestep(ii)+2,2));
    linem(pset.lat(1:month_this(1)-1,b(kk)),...
        pset.lon(1:month_this(1)-1,b(kk)),...
        'color',rgb('SkyBlue'),'LineWidth',0.8)
    linem(pset.lat(month_this(1)-1:month_this(end),b(kk)),...
        pset.lon(month_this(1)-1:month_this(end),b(kk)),...
        'color',rgb('MidnightBlue'),'LineWidth',0.8)
    %         scatterm(pset.lat(month_this(1),b(kk)),...
    %             pset.lon(month_this(1),b(kk)),...
    %             30,rgb('MidnightBlue'))
    %         scatterm(pset.lat(month_this(end),b(kk)),...
    %             pset.lon(month_this(end),b(kk)),...
    %             30,rgb('MidnightBlue'),'filled')
    clear month_this
end
clear a b c jj kk
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.05 35.05],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.1 35.1],...
%    'Color','k','LineStyle','-','LineWidth',2)
textm(46.3,-59.8,'(i)','Color',rgb('Black'),'FontSize',26)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
%geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
textm(46.3,-57.6,[datestr(time(slope_timestep(ii)+2),'mmm'),' ',...
    datestr(time(slope_timestep(ii)+2),'yyyy')],...
    'Color',rgb('Black'),'FontSize',26)
%     title([datestr(time(slope_timestep(ii)+2),'mmm'),' ',...
%         datestr(time(slope_timestep(ii)+2),'yyyy')],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')
saveas(gcf,[dir_output,'manuscript_v2/TGB2_Figure7_v2'],'png')

%% Figure 8

% Maps of Depth and Thickness (colors and contours), T or S, Speed or KE,
% averaged for 2004 and 2008.

dir_input = '/Users/afonso/Documents/Research/TGB2/Data/Eulerian/';
dir_output = '/Users/afonso/Documents/Research/TGB2/output/';

ncdisp([dir_output,'ATLg_008_layer19_monthly.nc'])

% lat and lon (high resolution)
lat_aux = ncread([dir_output,'hycom_atl_19_monthly.nc'],'lat');
lon_aux = ncread([dir_output,'hycom_atl_19_monthly.nc'],'lon');
%lat_aux = lat_aux';
%lon_aux = lon_aux' - 360;
%[lon1,lat1] = meshgrid(lon_aux,lat_aux);

% lat and lon (low resolution)
lat = ncread([dir_output,'ATLg008_eulerian_layer19.nc'],'latitude');
lon = ncread([dir_output,'ATLg008_eulerian_layer19.nc'],'longitude');
lat = lat';
lon = lon' - 360;


% Bathymetry
%ncdisp([dir_output,'depth_ATLg008_10.nc'])
bathy = ncread([dir_output,'depth_ATLg008_10.nc'],'bathymetry');
bathy(bathy >= 10^10) = NaN;% remove values on land
bathy = bathy';

% Create variable for time
year = 1993:2017;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(300,1);
for yy = 1:25
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

% Eulerian fields

d1 = dir('/Users/afonso/Documents/Research/TGB2/Data/Eulerian');
d1 = d1(4:end);
thick = [];
depth = [];
temp = [];
salt = [];
for ii=1:length(d1)
    load([dir_input,d1(ii).name]);
    thick = cat(3,thick,hh(817:1033,263:727));
    depth = cat(3,depth,zu(817:1033,263:727));
    temp = cat(3,temp,tt(817:1033,263:727));
    salt = cat(3,salt,ss(817:1033,263:727));
    clear hh ss tt zu
end
clear lat_aux lon_aux ii d1
temp(temp >= 10^10) = NaN;% remove values on land
temp(thick <= 1) = NaN;% remove values on shelf
salt(salt >= 10^10) = NaN;% remove values on land
salt(thick <= 1) = NaN;% remove values on shelf

% Potential thickness (top of layer to bottom of ocean)
f = gsw_f(lat);
f_40 = gsw_f(40);
thick_pot = NaN(size(depth));
for tt = 1:length(time)
    thick_pot(:,:,tt) = f_40 * (bathy - squeeze(depth(:,:,tt))) ./ f;
end

% U
u = ncread([dir_output,'ATLg_008_layer19_monthly.nc'],'uvel');
u = u(263:727,817:1033,:);
u = permute(u,[2,1,3]);
u(u >= 10^10) = NaN;% remove values on land
u(thick <= 1) = NaN;% remove values on shelf
u = u .* 100;% cm/s

% V
v = ncread([dir_output,'ATLg_008_layer19_monthly.nc'],'vvel');
v = v(263:727,817:1033,:);
v = permute(v,[2,1,3]);
v(v >= 10^10) = NaN;% remove values on land
v(thick <= 1) = NaN;% remove values on shelf
v = v .* 100;% cm/s

% Speed
uv = sqrt(u.^2 + v.^2);

file_traj = [dir_output,'ATLg_008_layer19_LC80'];

ncdisp([file_traj,'.nc'])

pset.trajectory = ncread([file_traj,'.nc'],'trajectory');
pset.lon = ncread([file_traj,'.nc'],'lon');
pset.lat = ncread([file_traj,'.nc'],'lat');
pset.time = ncread([file_traj,'.nc'],'time');
pset.time2 = pset.time/(60*60*24)+datenum(1993,1,1,12,0,0);

load([dir_output,'etopo1_TGB.mat'],'lat_z','lon_z','z')

% Create variable for time
year = 1993:2017;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(300,1);
for yy = 1:25
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

figure('Position', [100, 200, 800, 600])
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
hold on
c1_aux = contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
linem([36 36],[360-55 360-52],'k')
linem([46 46],[360-55 360-52],'k')
linem([36 46],[360-55 360-55],'k')
linem([36 46],[360-52 360-52],'k')
geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])

close all

c1 = NaN(2,1);
j = 1;
for i = 2:length(c1_aux)
    if c1_aux(1,i)>=360-55 & c1_aux(1,i)<=360-52 & c1_aux(2,i)>=36 & c1_aux(2,i)<=46
        c1(:,j) = c1_aux(:,i);
        if j>1
            if abs(c1_aux(1,i)-c1(1,j-1))<0.2
                j = j+1;
            else
                c1(:,j) = NaN;
            end
        else
            j = j+1;
        end
    end
end
clear i j

if isnan(c1(1,end))
    c1 = c1(:,1:end-1);% lon,lat of 1000m isobath
end

c1(3,:) = zeros;
c1_dist = gsw_distance(c1(1,:),c1(2,:));
for i = 2:length(c1)
    c1(3,i) = sum(c1_dist(1:i-1));
end
clear c1_aux c1_dist i

% Line 1 degree to the south of the 1,000-m isobath

c2 = c1;
c2(2,:) = c1(2,:) - 1;
c3 = c1;
c3(2,:) = c1(2,:) - 2;

% figure('Position', [100, 200, 800, 600])
% axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
% hold on
% contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem(c1(2,:),c1(1,:),'k','linewidth',2)
% linem(c2(2,:),c2(1,:),'--k','linewidth',1)
% linem(c3(2,:),c3(1,:),'k','linewidth',2)
% %contourm(c1(2,:),c1(2,:),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem([36 36],[360-75 360-52],'k')
% linem([46 46],[360-75 360-52],'k')
% linem([36 46],[360-75 360-75],'k')
% linem([36 46],[360-52 360-52],'k')
% geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
% 
% close all

close all

% Define run time
time_run = 1;%years
time_total = 25;%years
%Subset
time_length = floor(size(pset.time,1)*(time_run/time_total));
traj_length = floor(size(pset.time,2)*((time_total-time_run)/time_total));

traj_slope = NaN(traj_length,1);
for ii = 1:traj_length
    traj_slope(ii,1) = sum(inpolygon(pset.lon(1:90,ii),pset.lat(1:90,ii),...
        [c1(1,:)-360 fliplr(c3(1,:)-360) c1(1,1)-360],...
        [c1(2,:) fliplr(c3(2,:)) c1(2,1)]));
end
clear ii

traj_slope1 = [traj_slope;zeros(size(pset.time2,2) - length(traj_slope),1)];

% Number of particles per release time (time series)
slope_releasetime = pset.time2(1,traj_slope>=1);
slope_releasetime_yr = str2num(datestr(slope_releasetime,'yyyy'));
slope_releasetime_mo = str2num(datestr(slope_releasetime,'mm'));
perc1_slope = 100*length(slope_releasetime)/traj_length; % 15% of floats make it to slope


yr = 1993:2016;
slope_releasetime_yr1 = NaN(length(yr),1);
for ii = 1:length(yr)
    slope_releasetime_yr1(ii,1) = sum(slope_releasetime_yr==yr(ii));
end
perc_yr_slope = 100*slope_releasetime_yr1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(yr)
    perc_yr1_slope(ii,1) = 100*slope_releasetime_yr1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii));
end

mo = 1:12;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

slope_releasetime_mo1 = NaN(length(mo),1);
for ii = 1:length(mo)
    slope_releasetime_mo1(ii,1) = sum(slope_releasetime_mo==mo(ii));
end
perc_mo_slope = 100*slope_releasetime_mo1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(mo)
    perc_mo1_slope(ii,1) = 100*slope_releasetime_mo1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(ii));
end

% Area along the NWA Slope
lat_nwa = pset.lat(1:time_length,traj_slope>=1);
lon_nwa = pset.lon(1:time_length,traj_slope>=1);
time2_nwa = pset.time2(1:time_length,traj_slope>=1);

% Releases per month
slope_release_number = NaN(288,5);
for ii = 1:length(yr)
    for jj = 1:length(mo)
        slope_release_number(12*(ii-1)+jj,1) = yr(ii);
        slope_release_number(12*(ii-1)+jj,2) = mo(jj);
        slope_release_number(12*(ii-1)+jj,3) = ...
            sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii)...
            & str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(jj));
        slope_release_number(12*(ii-1)+jj,4) = sum(str2num(datestr(time2_nwa(1,:),'mm'))...
            ==mo(jj) & str2num(datestr(time2_nwa(1,:),'yyyy'))==yr(ii));
        slope_release_number(12*(ii-1)+jj,5) = 100 * slope_release_number(12*(ii-1)+jj,4)...
            /slope_release_number(12*(ii-1)+jj,3);
    end
end

slope_timestep1 = find(slope_release_number(:,5)>=20);
slope_timestep2 = find(slope_release_number(:,5)==0);
csvwrite([dir_output,'valveopen_months.csv'],slope_timestep1)
csvwrite([dir_output,'valveclosed_months.csv'],slope_timestep2)

% percentage of particles released in open valve months that made it
nansum(slope_release_number(slope_timestep1,3))
nansum(slope_release_number(slope_timestep1,4))
100*nansum(slope_release_number(slope_timestep1,4))./nansum(slope_release_number(slope_timestep1,3))
nanmean(100*slope_release_number(slope_timestep1,4)./slope_release_number(slope_timestep1,3))
nanmean(slope_release_number(slope_timestep1,5))

% split in two matrices to plot quiver in different colors
u_low = NaN(size(u));
u_high = NaN(size(u));
v_low = NaN(size(v));
v_high = NaN(size(v));

u_low(uv<=12) = u(uv<=12);
u_high(uv>12) = u(uv>12);
v_low(uv<=12) = v(uv<=12);
v_high(uv>12) = v(uv>12);

% plot_location = [0 0.64 0.5 0.32;...
%     0.5 0.64 0.5 0.32;...
%     0 0.32 0.5 0.32;...
%     0.5 0.32 0.5 0.32;...
%     0 0 0.5 0.32;...
%     0.5 0 0.5 0.32];

plot_location = [0.05 0.66 0.455 0.29;...
    0.515 0.66 0.455 0.29;...
    0.05 0.36 0.455 0.29;...
    0.515 0.36 0.455 0.29;...
    0.05 0.06 0.455 0.29;...
    0.515 0.06 0.455 0.29];

cbar_location = [0.057 0.84 0.098 0.03;...
    0.522 0.84 0.098 0.03;...
    0.057 0.54 0.098 0.03;...
    0.522 0.54 0.098 0.03;...
    0.057 0.24 0.098 0.03;...
    0.522 0.24 0.098 0.03];


colormap_aux = readtable([dir_output,'5w_BRgpb.csv'], 'HeaderLines',1);
colormap_wave = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'4w_ROTB.csv'], 'HeaderLines',1);
colormap_wave2 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'3-wave-yellow-grey-blue.csv'], 'HeaderLines',1);
colormap_wave3 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'5-step-melow-wave.csv'], 'HeaderLines',1);
colormap_wave4 = table2array(colormap_aux(:,3:5));
clear colormap_aux 

colormap_aux = readtable([dir_output,'turqoise-olive.csv'], 'HeaderLines',1);
colormap_wave5 = table2array(colormap_aux(:,3:5));
clear colormap_aux 

figure('Position', [0, 0, 1102, 740])
s1 = subplot('Position',plot_location(1,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
pcolorm(lat,lon,squeeze(nanmean(thick_pot(:,:,slope_timestep1+1),3)) - ...
    squeeze(nanmean(thick_pot,3)))
shading interp
%colormap(cmocean('speed'))
colormap(s1,flipud(colormap_wave5))
c = colorbar('Orientation','Horizontal','Position',cbar_location(1,:),...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'H Pot anom (m)';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([-90 90])
contourm(lat,lon,bathy,[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.05 35.05],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.1 35.1],...
%    'Color','k','LineStyle','-','LineWidth',2)
%textm(40,-79.4,'(c)','Color',rgb('Snow'),'FontSize',23)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
title(['Slope >= 20% (N = ',num2str(length(slope_timestep1)),')'],'FontSize',22)
plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap

s2 = subplot('Position',plot_location(2,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
pcolorm(lat,lon,squeeze(nanmean(thick_pot(:,:,slope_timestep2+1),3)) - ...
    squeeze(nanmean(thick_pot,3)))
shading interp
%colormap(cmocean('speed'))
colormap(s2,flipud(colormap_wave5))
c = colorbar('Orientation','Horizontal','Position',cbar_location(2,:),...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'H Pot anom (m)';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([-90 90])
contourm(lat,lon,bathy,[1000 1000],...
    'Color',rgb('DarkSlateGray'),'LineStyle','-','LineWidth',1)
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.05 35.05],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.1 35.1],...
%    'Color','k','LineStyle','-','LineWidth',2)
%textm(40,-79.4,'(c)','Color',rgb('Snow'),'FontSize',23)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
title(['Slope = 0% (N = ',num2str(length(slope_timestep2)),')'],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap

s3 = subplot('Position',plot_location(3,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
pcolorm(lat,lon,squeeze(nanmean(salt(:,:,slope_timestep1+1),3)))
shading interp
%colormap(cmocean('speed'))
colormap(s3,flipud(colormap_wave2))
c = colorbar('Orientation','Horizontal','Position',cbar_location(3,:),...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'Salinity';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([34.75 35.25])
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.05 35.05],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.1 35.1],...
%    'Color','k','LineStyle','-','LineWidth',2)
%textm(40,-79.4,'(c)','Color',rgb('Snow'),'FontSize',23)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%title([datestr(time(slope_timestep(ii)+2),'mmm'),' ',...
%    datestr(time(slope_timestep(ii)),'yyyy')],'FontSize',22)
plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap

s4 = subplot('Position',plot_location(4,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
pcolorm(lat,lon,squeeze(nanmean(salt(:,:,slope_timestep2+1),3)))
shading interp
%colormap(cmocean('speed'))
colormap(s4,flipud(colormap_wave2))
c = colorbar('Orientation','Horizontal','Position',cbar_location(4,:),...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'Salinity';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([34.75 35.25])
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.05 35.05],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(salt(2*12-11:2*12,:,:),1))',[35.1 35.1],...
%    'Color','k','LineStyle','-','LineWidth',2)
%textm(40,-79.4,'(c)','Color',rgb('Snow'),'FontSize',23)%,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%title([datestr(time(slope_timestep(ii)+2),'mmm'),' ',...
%    datestr(time(slope_timestep(ii)),'yyyy')],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap

u2 = squeeze(nanmean(u(:,:,slope_timestep1+1),3));
v2 = squeeze(nanmean(v(:,:,slope_timestep1+1),3));
uv2 = sqrt(u2.^2 + v2.^2);

% split in two matrices to plot quiver in different colors
u_low2 = squeeze(nanmean(u(:,:,slope_timestep1+1),3));
v_low2 = squeeze(nanmean(v(:,:,slope_timestep1+1),3));
uv_low2 = sqrt(u_low2.^2 + v_low2.^2);

u_low2(uv_low2>6) = NaN;
v_low2(uv_low2>6) = NaN;
clear uv_low*

u_high2 = squeeze(nanmean(u(:,:,slope_timestep1+1),3));
v_high2 = squeeze(nanmean(v(:,:,slope_timestep1+1),3));
uv_high2 = sqrt(u_high2.^2 + v_high2.^2);

u_high2(uv_high2<=6) = NaN;
v_high2(uv_high2<=6) = NaN;
clear uv_high*

s5 = subplot('Position',plot_location(5,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
pcolorm(lat,lon,uv2)
shading interp
%colormap(cmocean('speed'))
colormap(s5,flipud(colormap_wave4))
c = colorbar('Orientation','Horizontal','Position',cbar_location(5,:),...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'Speed (cm s^{-1})';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([0 20])
alpha 0.5
q1 = quivermc(lat,lon,...
    u_low2,...
    v_low2,...
    'reference',6,'density',100/18);
set(q1,'LineWidth',2,'Color',rgb('DarkSlateGray'))%,'MaxHeadSize',0.5)
q2 = quivermc(lat,lon,...
    u_high2,...
    v_high2,...
    'reference',15,'density',100/18);
set(q2,'LineWidth',2,'Color',rgb('DarkRed'))%,'MaxHeadSize',0.5)
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[450 450],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[700 700],...
%    'Color','k','LineStyle','-','LineWidth',2)
%textm(40,-79.4,'(a)','Color',rgb('Snow'),'FontSize',23)%,'Fontweight','bold')
%textm(42.5,-74.8,'19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%legend
f1 = fillm([37,41,41,37],[-79.7,-79.7,-73.5,-73.5],'w');
f1.FaceColor = 'w';
f1.EdgeColor = 'w';
f1.FaceAlpha = 0.75;
f1.EdgeAlpha = 1;
q4 = quivermc([40 40;39 39],[-79 -76;-79 -76],[6 0; 0 0],[0 0; 0 0],...
    'reference',6);
set(q4,'LineWidth',2,'Color',rgb('DarkSlateGray'))%,'MaxHeadSize',0.5)
q5 = quivermc([38 38;37 37],[-79 -76;-79 -76],[12 0; 0 0],[0 0; 0 0],...
    'reference',15);
set(q5,'LineWidth',2,'Color',rgb('DarkRed'))%,'MaxHeadSize',0.5)
textm(40.3,-78.3,'6 cm s^{-1}','Color',rgb('Black'),'FontSize',12)%,'Fontweight','bold')
textm(38.2,-78.3,'15 cm s^{-1}','Color',rgb('Black'),'FontSize',12)%,'Fontweight','bold')
f2 = fillm([41.5,44.5,44.5,41.5],[-52,-52,-48,-48],'w','LineWidth',3);
f2.FaceColor = 'w';
f2.EdgeColor = 'k';
f2.FaceAlpha = 0;
f2.EdgeAlpha = 1;
%title([datestr(time(slope_timestep(ii)),'mmm'),' ',...
%    datestr(time(slope_timestep(ii)),'yyyy')],'FontSize',22)
plabel('PLabelLocation',5,'Fontsize',18)
mlabel('MLabelLocation',10,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

u2 = squeeze(nanmean(u(:,:,slope_timestep2+1),3));
v2 = squeeze(nanmean(v(:,:,slope_timestep2+1),3));
uv2 = sqrt(u2.^2 + v2.^2);

% split in two matrices to plot quiver in different colors
u_low2 = squeeze(nanmean(u(:,:,slope_timestep2+1),3));
v_low2 = squeeze(nanmean(v(:,:,slope_timestep2+1),3));
uv_low2 = sqrt(u_low2.^2 + v_low2.^2);

u_low2(uv_low2>6) = NaN;
v_low2(uv_low2>6) = NaN;
clear uv_low*

u_high2 = squeeze(nanmean(u(:,:,slope_timestep2+1),3));
v_high2 = squeeze(nanmean(v(:,:,slope_timestep2+1),3));
uv_high2 = sqrt(u_high2.^2 + v_high2.^2);

u_high2(uv_high2<=6) = NaN;
v_high2(uv_high2<=6) = NaN;
clear uv_high*

s6 = subplot('Position',plot_location(6,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
pcolorm(lat,lon,uv2)
shading interp
%colormap(cmocean('speed'))
colormap(s6,flipud(colormap_wave4))
c = colorbar('Orientation','Horizontal','Position',cbar_location(6,:),...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'Speed (cm s^{-1})';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([0 20])
alpha 0.5
q1 = quivermc(lat,lon,...
    u_low2,...
    v_low2,...
    'reference',6,'density',100/18);
set(q1,'LineWidth',2,'Color',rgb('DarkSlateGray'))%,'MaxHeadSize',0.5)
q2 = quivermc(lat,lon,...
    u_high2,...
    v_high2,...
    'reference',15,'density',100/18);
set(q2,'LineWidth',2,'Color',rgb('DarkRed'))%,'MaxHeadSize',0.5)
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[450 450],...
%    'Color','k','LineStyle','--','LineWidth',2)
%contourm(lat,lon,squeeze(nanmean(depth(2*12-11:2*12,:,:),1))',[700 700],...
%    'Color','k','LineStyle','-','LineWidth',2)
%textm(40,-79.4,'(a)','Color',rgb('Snow'),'FontSize',23)%,'Fontweight','bold')
%textm(42.5,-74.8,'19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%legend
f1 = fillm([37,41,41,37],[-79.7,-79.7,-73.5,-73.5],'w');
f1.FaceColor = 'w';
f1.EdgeColor = 'w';
f1.FaceAlpha = 0.75;
f1.EdgeAlpha = 1;
q4 = quivermc([40 40;39 39],[-79 -76;-79 -76],[6 0; 0 0],[0 0; 0 0],...
    'reference',6);
set(q4,'LineWidth',2,'Color',rgb('DarkSlateGray'))%,'MaxHeadSize',0.5)
q5 = quivermc([38 38;37 37],[-79 -76;-79 -76],[12 0; 0 0],[0 0; 0 0],...
    'reference',15);
set(q5,'LineWidth',2,'Color',rgb('DarkRed'))%,'MaxHeadSize',0.5)
textm(40.3,-78.3,'6 cm s^{-1}','Color',rgb('Black'),'FontSize',12)%,'Fontweight','bold')
textm(38.2,-78.3,'15 cm s^{-1}','Color',rgb('Black'),'FontSize',12)%,'Fontweight','bold')
f2 = fillm([41.5,44.5,44.5,41.5],[-52,-52,-48,-48],'w','LineWidth',3);
f2.FaceColor = 'w';
f2.EdgeColor = 'k';
f2.FaceAlpha = 0;
f2.EdgeAlpha = 1;
%title([datestr(time(slope_timestep(ii)),'mmm'),' ',...
%    datestr(time(slope_timestep(ii)),'yyyy')],'FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',18)
mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',18)
framem
tightmap
set(gcf, 'Color', 'w')
saveas(gcf,[dir_output,'manuscript_v3/TGB2_Figure8_v4'],'png')

%% Figure 9

% depth, thick, pot_thick, SSH (2008 minus 2004)

dir_input = '/Users/afonso/Documents/Research/TGB2/Data/Eulerian/';
dir_output = '/Users/afonso/Documents/Research/TGB2/output/';

ncdisp([dir_output,'ATLg_008_layer19_monthly.nc'])

% lat and lon (high resolution)
lat_aux = ncread([dir_output,'hycom_atl_19_monthly.nc'],'lat');
lon_aux = ncread([dir_output,'hycom_atl_19_monthly.nc'],'lon');
%lat_aux = lat_aux';
%lon_aux = lon_aux' - 360;
%[lon1,lat1] = meshgrid(lon_aux,lat_aux);

% lat and lon (low resolution)
lat = ncread([dir_output,'ATLg008_eulerian_layer19.nc'],'latitude');
lon = ncread([dir_output,'ATLg008_eulerian_layer19.nc'],'longitude');
lat = lat';
lon = lon' - 360;


% Bathymetry
%ncdisp([dir_output,'depth_ATLg008_10.nc'])
bathy = ncread([dir_output,'depth_ATLg008_10.nc'],'bathymetry');
bathy(bathy >= 10^10) = NaN;% remove values on land
bathy = bathy';

% Create variable for time
year = 1993:2017;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(300,1);
for yy = 1:25
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

% Eulerian fields

d1 = dir('/Users/afonso/Documents/Research/TGB2/Data/Eulerian');
d1 = d1(4:end);
thick = [];
depth = [];
temp = [];
salt = [];
for ii=1:length(d1)
    load([dir_input,d1(ii).name]);
    thick = cat(3,thick,hh(817:1033,263:727));
    depth = cat(3,depth,zu(817:1033,263:727));
    clear hh ss tt zu
end
clear lat_aux lon_aux ii d1

% Potential thickness (top of layer to bottom of ocean)
f = gsw_f(lat);
f_40 = gsw_f(40);
thick_pot = NaN(size(depth));
for tt = 1:length(time)
    thick_pot(:,:,tt) = f_40 * (bathy - squeeze(depth(:,:,tt))) ./ f;
end

% SSH
load('/Users/afonso/Documents/Research/TGB2/Data/map_GSS_ATLg08.mat')
dir_list = dir('/Users/afonso/Documents/Research/TGB2/Data/Hycom_ATLg08/*.mat');
ssh = NaN(size(x,1),size(x,2),length(dir_list));
for ii = 1:length(dir_list)
    a = load(['/Users/afonso/Documents/Research/TGB2/Data/Hycom_ATLg08/',dir_list(ii).name]);
    ssh(:,:,ii) = a.ssh;
    clear a
end
clear ii x y dir_list pdep
ssh = ssh(1:215,:,:);

file_traj = [dir_output,'ATLg_008_layer19_LC80'];

ncdisp([file_traj,'.nc'])

pset.trajectory = ncread([file_traj,'.nc'],'trajectory');
pset.lon = ncread([file_traj,'.nc'],'lon');
pset.lat = ncread([file_traj,'.nc'],'lat');
pset.time = ncread([file_traj,'.nc'],'time');
pset.time2 = pset.time/(60*60*24)+datenum(1993,1,1,12,0,0);

load([dir_output,'etopo1_TGB.mat'],'lat_z','lon_z','z')

% Create variable for time
year = 1993:2017;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
kk = 0;
time = NaN(300,1);
for yy = 1:25
    for mm = 1:12
        kk = kk + 1;
        if mm < 10
            time(kk,1) = datenum([num2str(year(yy)),'0',num2str(mm),'15'],...
                'yyyymmdd');
        else
            time(kk,1) = datenum([num2str(year(yy)),num2str(mm),'15'],...
                'yyyymmdd');
        end
    end
end
clear kk yy mm

figure('Position', [100, 200, 800, 600])
axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
hold on
c1_aux = contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
linem([36 36],[360-55 360-52],'k')
linem([46 46],[360-55 360-52],'k')
linem([36 46],[360-55 360-55],'k')
linem([36 46],[360-52 360-52],'k')
geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])

close all

c1 = NaN(2,1);
j = 1;
for i = 2:length(c1_aux)
    if c1_aux(1,i)>=360-55 & c1_aux(1,i)<=360-52 & c1_aux(2,i)>=36 & c1_aux(2,i)<=46
        c1(:,j) = c1_aux(:,i);
        if j>1
            if abs(c1_aux(1,i)-c1(1,j-1))<0.2
                j = j+1;
            else
                c1(:,j) = NaN;
            end
        else
            j = j+1;
        end
    end
end
clear i j

if isnan(c1(1,end))
    c1 = c1(:,1:end-1);% lon,lat of 1000m isobath
end

c1(3,:) = zeros;
c1_dist = gsw_distance(c1(1,:),c1(2,:));
for i = 2:length(c1)
    c1(3,i) = sum(c1_dist(1:i-1));
end
clear c1_aux c1_dist i

% Line 1 degree to the south of the 1,000-m isobath

c2 = c1;
c2(2,:) = c1(2,:) - 1;
c3 = c1;
c3(2,:) = c1(2,:) - 2;

% figure('Position', [100, 200, 800, 600])
% axesm('MapProjection','eqdcylin','MapLatLimit',[30 60],'MapLonLimit',[-80 -40])
% hold on
% contourm(lat_z(1:10:end),lon_z(1:10:end),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem(c1(2,:),c1(1,:),'k','linewidth',2)
% linem(c2(2,:),c2(1,:),'--k','linewidth',1)
% linem(c3(2,:),c3(1,:),'k','linewidth',2)
% %contourm(c1(2,:),c1(2,:),z(1:10:end,1:10:end),[-1000 -1000],'k','linewidth',1);
% linem([36 36],[360-75 360-52],'k')
% linem([46 46],[360-75 360-52],'k')
% linem([36 46],[360-75 360-75],'k')
% linem([36 46],[360-52 360-52],'k')
% geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
% 
% close all

close all

% Define run time
time_run = 1;%years
time_total = 25;%years
%Subset
time_length = floor(size(pset.time,1)*(time_run/time_total));
traj_length = floor(size(pset.time,2)*((time_total-time_run)/time_total));

traj_slope = NaN(traj_length,1);
for ii = 1:traj_length
    traj_slope(ii,1) = sum(inpolygon(pset.lon(1:90,ii),pset.lat(1:90,ii),...
        [c1(1,:)-360 fliplr(c3(1,:)-360) c1(1,1)-360],...
        [c1(2,:) fliplr(c3(2,:)) c1(2,1)]));
end
clear ii

traj_slope1 = [traj_slope;zeros(size(pset.time2,2) - length(traj_slope),1)];

% Number of particles per release time (time series)
slope_releasetime = pset.time2(1,traj_slope>=1);
slope_releasetime_yr = str2num(datestr(slope_releasetime,'yyyy'));
slope_releasetime_mo = str2num(datestr(slope_releasetime,'mm'));
perc1_slope = 100*length(slope_releasetime)/traj_length; % 15% of floats make it to slope


yr = 1993:2016;
slope_releasetime_yr1 = NaN(length(yr),1);
for ii = 1:length(yr)
    slope_releasetime_yr1(ii,1) = sum(slope_releasetime_yr==yr(ii));
end
perc_yr_slope = 100*slope_releasetime_yr1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(yr)
    perc_yr1_slope(ii,1) = 100*slope_releasetime_yr1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii));
end

mo = 1:12;
month_str = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';...
    'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

slope_releasetime_mo1 = NaN(length(mo),1);
for ii = 1:length(mo)
    slope_releasetime_mo1(ii,1) = sum(slope_releasetime_mo==mo(ii));
end
perc_mo_slope = 100*slope_releasetime_mo1./length(slope_releasetime); % 15% of floats make it to slope

for ii = 1:length(mo)
    perc_mo1_slope(ii,1) = 100*slope_releasetime_mo1(ii,1)/...
        sum(str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(ii));
end

% Area along the NWA Slope
lat_nwa = pset.lat(1:time_length,traj_slope>=1);
lon_nwa = pset.lon(1:time_length,traj_slope>=1);
time2_nwa = pset.time2(1:time_length,traj_slope>=1);

% Releases per month
slope_release_number = NaN(288,5);
for ii = 1:length(yr)
    for jj = 1:length(mo)
        slope_release_number(12*(ii-1)+jj,1) = yr(ii);
        slope_release_number(12*(ii-1)+jj,2) = mo(jj);
        slope_release_number(12*(ii-1)+jj,3) = ...
            sum(str2num(datestr(pset.time2(1,1:traj_length),'yyyy'))==yr(ii)...
            & str2num(datestr(pset.time2(1,1:traj_length),'mm'))==mo(jj));
        slope_release_number(12*(ii-1)+jj,4) = sum(str2num(datestr(time2_nwa(1,:),'mm'))...
            ==mo(jj) & str2num(datestr(time2_nwa(1,:),'yyyy'))==yr(ii));
        slope_release_number(12*(ii-1)+jj,5) = 100 * slope_release_number(12*(ii-1)+jj,4)...
            /slope_release_number(12*(ii-1)+jj,3);
    end
end

slope_timestep1 = find(slope_release_number(:,5)>=20);
slope_timestep2 = find(slope_release_number(:,5)==0);

plot_location = [0.05 0.53 0.455 0.46;...
    0.515 0.53 0.455 0.46;...
    0.05 0.07 0.455 0.46;...
    0.515 0.07 0.455 0.46];

cbar_location = [0.06 0.83 0.083 0.02;...
    0.53 0.83 0.083 0.02;...
    0.06 0.37 0.083 0.02;...
    0.53 0.37 0.083 0.02];

colormap_aux = readtable([dir_output,'5w_BRgpb.csv'], 'HeaderLines',1);
colormap_wave = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'4w_ROTB.csv'], 'HeaderLines',1);
colormap_wave2 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'3-wave-yellow-grey-blue.csv'], 'HeaderLines',1);
colormap_wave3 = table2array(colormap_aux(:,3:5));
clear colormap_aux

colormap_aux = readtable([dir_output,'5-step-melow-wave.csv'], 'HeaderLines',1);
colormap_wave4 = table2array(colormap_aux(:,3:5));
clear colormap_aux 

figure('Position', [0, 0, 1160, 480])
s1 = subplot('Position',plot_location(1,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
%pcolorm(lat_hycom_atl,lon_hycom_atl,squeeze(nanmean(ssh_hycom_atl(:,:,16*12-11:16*12),3)) - ...
%    squeeze(nanmean(ssh_hycom_atl(:,:,12*12-11:12*12),3)))
%pcolorm(lat,lon,10*squeeze(nanmean(ssh(6*12-11:6*12,:,:),1))' - ...
%    10*squeeze(nanmean(ssh(2*12-11:2*12,:,:),1))')
pcolorm(lat(2:end-1,2:end-1),lon(2:end-1,2:end-1),squeeze(nanmean(ssh(:,:,slope_timestep2),3)) - ...
    squeeze(nanmean(ssh(:,:,slope_timestep1),3)))
shading interp
%colormap(cmocean('speed'))
colormap(s1,flipud(colormap_wave3(3:43,:)))
c = colorbar('Orientation','Horizontal','Position',cbar_location(1,:),...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'SSH (cm)';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([-30 30])
contourm(lat,lon,bathy,[200 200],'Color',rgb('DarkSlateGray'),...
    'LineWidth',1.5)
textm(40,-79,'(a)','Color',rgb('Snow'),'FontSize',20,'Fontweight','bold')
f1 = fillm([41,43,43,41],[-53,-53,-50,-50],'w');
f1.FaceColor = 'w';
f1.EdgeColor = 'k';
f1.FaceAlpha = 0;
f1.EdgeAlpha = 1;
f1.LineWidth = 2.5;
%textm(42.5,-74.8,'19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%title('2004','FontSize',22)
plabel('PLabelLocation',5,'Fontsize',18)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

s2 = subplot('Position',plot_location(2,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
pcolorm(lat,lon,squeeze(nanmean(depth(:,:,slope_timestep2),3)) - ...
    squeeze(nanmean(depth(:,:,slope_timestep1),3)))
shading interp
%colormap(cmocean('speed'))
colormap(s2,flipud(colormap_wave3(3:43,:)))
c = colorbar('Orientation','Horizontal','Position',cbar_location(2,:),...
    'Color',[1-eps,1,1],'Ticks',[-100;0;100]);%colorbar over land
c.Label.String = 'Depth (m)';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([-150 150])
contourm(lat,lon,bathy,[200 200],'Color',rgb('DarkSlateGray'),...
    'LineWidth',1.5)
textm(40,-79,'(b)','Color',rgb('Snow'),'FontSize',20,'Fontweight','bold')
f1 = fillm([41,43,43,41],[-53,-53,-50,-50],'w');
f1.FaceColor = 'w';
f1.EdgeColor = 'k';
f1.FaceAlpha = 0;
f1.EdgeAlpha = 1;
f1.LineWidth = 2.5;
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%title('2008','FontSize',22)
%plabel('PLabelLocation',5,'Fontsize',14)
%mlabel('MLabelLocation',5,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

s3 = subplot('Position',plot_location(3,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
pcolorm(lat,lon,squeeze(nanmean(thick_pot(:,:,slope_timestep2),3)) - ...
    squeeze(nanmean(thick_pot(:,:,slope_timestep1),3)))
shading interp
%colormap(cmocean('speed'))
colormap(s3,flipud(colormap_wave3(3:43,:)))
c = colorbar('Orientation','Horizontal','Position',cbar_location(3,:),...
    'Color',[1-eps,1,1],'Ticks',[-100;0;100]);%colorbar over land
c.Label.String = 'Pot Thick (m)';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([-150 150])
contourm(lat,lon,bathy,[200 200],'Color',rgb('DarkSlateGray'),...
    'LineWidth',1.5)
textm(40,-79,'(c)','Color',rgb('Snow'),'FontSize',20,'Fontweight','bold')
f1 = fillm([41,43,43,41],[-53,-53,-50,-50],'w');
f1.FaceColor = 'w';
f1.EdgeColor = 'k';
f1.FaceAlpha = 0;
f1.EdgeAlpha = 1;
f1.LineWidth = 2.5;
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
plabel('PLabelLocation',5,'Fontsize',18)
mlabel('MLabelLocation',10,'MLabelParallel','south')
framem
tightmap
set(gcf, 'Color', 'w')

s4 = subplot('Position',plot_location(4,:));
axesm('MapProjection','eqdcylin','MapLatLimit',[34 47],'MapLonLimit',[-79.9 -40])
hold on
pcolorm(lat,lon,squeeze(nanmean(thick(:,:,slope_timestep2),3)) - ...
    squeeze(nanmean(thick(:,:,slope_timestep1),3)))
shading interp
%colormap(cmocean('speed'))
colormap(s4,flipud(colormap_wave3(3:43,:)))
c = colorbar('Orientation','Horizontal','Position',cbar_location(4,:),...
    'Color',[1-eps,1,1]);%colorbar over land
c.Label.String = 'Thickness (m)';
set(c,'Fontsize',18)%,'Fontweight','bold')
caxis([-20 20])
contourm(lat,lon,bathy,[200 200],'Color',rgb('DarkSlateGray'),...
    'LineWidth',1.5)
textm(40,-79,'(d)','Color',rgb('Snow'),'FontSize',20,'Fontweight','bold')
%textm(43.5,-74.8,'Layer 19','Color','w','Fontsize',20,'Fontweight','bold')
%alpha 0.5
%quiverm(lat_mesh(1:6:end,1:6:end),lon_mesh(1:6:end,1:6:end),...
%    nanmean(v_19(1:6:end,1:6:end,:),3)',...
%    nanmean(u_19(1:6:end,1:6:end,:),3)',2)
geoshow('GSHHS_l_L1.shp','FaceColor','k')%[.7 .7 .7])
%geoshow('landareas.shp','FaceColor','k')%[.7 .7 .7])
%title('HYCOM Atl Layer 20 - Mean Speed (2003-2012)')
%plabel('PLabelLocation',5,'Fontsize',14)
mlabel('MLabelLocation',10,'MLabelParallel','south','Fontsize',18)
framem
tightmap
set(gcf, 'Color', 'w')
saveas(gcf,[dir_output,'manuscript_v3/TGB2_Figure9_v4'],'png')
