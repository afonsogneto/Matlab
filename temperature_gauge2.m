%% temperature_gauge2

% Afonso Gonçalves Neto
% January 8, 2022

%% Global 1880-present

% Input
% [min]: minimum value (scalar)
% [max]: maximum value (scalar)
% [datapoint]: datapoint value (scalar)

% Output
% Plot of gauge indicating datapoint

dir_input = '/Users/afonso/Documents/Research/Birthday/';

% Time 
time_aux = ncread([dir_input,'gistemp250_GHCNv4.nc'],'time');
time = double(time_aux);
time_str = string(datestr(time + datenum('18000101','yyyymmdd'),'yyyymm'));
clear time_aux

time1 = ncread([dir_input,'gistemp250_GHCNv4.nc'],'time');
time = double(time1) + datenum('18000101','yyyymmdd');
clear time1

% Lat and Lon
lat = ncread([dir_input,'gistemp250_GHCNv4.nc'],'lat');
lon = ncread([dir_input,'gistemp250_GHCNv4.nc'],'lon');

% Temperature
T_aux = ncread([dir_input,'gistemp250_GHCNv4.nc'],'tempanomaly');
T = squeeze(nanmean(nanmean(T_aux,1),2));
clear T_aux

% Data limits
%min_value = min(T);
%max_value = max(T);
min_value = -2.2;
max_value = 2.2;

% Colormap
dir_colormap = '/Users/afonso/Documents/StudentSeminar 2/';
colormap_aux = readtable([dir_colormap,'w5m4.csv'], 'HeaderLines',1);
cbar_diverge2 = table2array(colormap_aux(:,3:5));
clear colormap_aux dir_colormap
cbar_stripes = NaN(450,3);
for ii = 1:3
    cbar_stripes(:,ii) = interp1((1:length(cbar_diverge2(2:end,ii)))',...
        cbar_diverge2(2:end,ii),(1:1/56.2:length(cbar_diverge2(2:end,ii)))');
end
cbar_stripes = flipud(cbar_stripes);
clear ii cbar_diverge2

% Info
radius = 1.04;
radius_out = 1.08;
radius_in = 0.85;
radius_ticks1 = 0.96;
radius_ticks2 = 1.01;
radius_base = 0.05;
CenterX = 0;
CenterY = 0;

% Gauge
angles1 = linspace(-pi/4, 5*pi/4, 450);
x1 = radius * cos(angles1) + CenterX;
y1 = radius * sin(angles1) + CenterY;

% Panel
angles2 = linspace(0, 2*pi, 450);
x2 = radius * cos(angles2) + CenterX;
y2 = radius * sin(angles2) + CenterY;

% Out
angles3 = linspace(0, 2*pi, 450);
x3 = radius_out * cos(angles3) + CenterX;
y3 = radius_out * sin(angles3) + CenterY;

% Stripes
angles_stripes = NaN(450,3);
x_stripes = NaN(450,3);
y_stripes = NaN(450,3);
for ii = 1:450
    angles_stripes(ii,1:3) = linspace((ii-1)*pi/360 - pi/4,(ii)*pi/360,3);
    x_stripes(ii,1:3) = radius * cos(angles_stripes(ii,:)) + CenterX;
    y_stripes(ii,1:3) = radius * sin(angles_stripes(ii,:)) + CenterY;
end

% Base
angles_base = linspace(0, 2*pi, 450);
x_base = radius_base * cos(angles_base) + CenterX;
y_base = radius_base * sin(angles_base) + CenterY;

% Ticks1
ticks1 = [-2:0.5:2];
ticks1_str = {'-2.0';'-1.5';'-1.0';'-0.5';'0.0';'0.5';'1.0';'1.5';'2.0'};
angles_ticks1 = NaN(length(ticks1),1);
x_ticks1_out = NaN(length(ticks1),1);
y_ticks1_out = NaN(length(ticks1),1);
x_ticks1_in = NaN(length(ticks1),1);
y_ticks1_in = NaN(length(ticks1),1);
x_in = NaN(length(ticks1),1);
y_in = NaN(length(ticks1),1);
for ii = 1:length(ticks1)
    angles_ticks1(ii,1) = 1.5*pi*(1-(ticks1(ii)-min_value)/(max_value-min_value)) - 0.25*pi;
    x_ticks1_out(ii,1) = radius * cos(angles_ticks1(ii)) + CenterX;
    y_ticks1_out(ii,1) = radius * sin(angles_ticks1(ii)) + CenterY;
    x_ticks1_in(ii,1) = radius_ticks1 * cos(angles_ticks1(ii)) + CenterX;
    y_ticks1_in(ii,1) = radius_ticks1 * sin(angles_ticks1(ii)) + CenterY;
    x_in(ii,1) = radius_in * cos(angles_ticks1(ii)) + CenterX;
    y_in(ii,1) = radius_in * sin(angles_ticks1(ii)) + CenterY;
end

% Ticks2
ticks2 = [-2.2:0.1:2.2];
angles_ticks2 = NaN(length(ticks2),1);
x_ticks2_out = NaN(length(ticks2),1);
y_ticks2_out = NaN(length(ticks2),1);
x_ticks2_in = NaN(length(ticks2),1);
y_ticks2_in = NaN(length(ticks2),1);
for ii = 1:length(ticks2)
    angles_ticks2(ii,1) = 1.5*pi*(1-(ticks2(ii)-min_value)/(max_value-min_value)) - 0.25*pi;
    x_ticks2_out(ii,1) = radius * cos(angles_ticks2(ii)) + CenterX;
    y_ticks2_out(ii,1) = radius * sin(angles_ticks2(ii)) + CenterY;
    x_ticks2_in(ii,1) = radius_ticks2 * cos(angles_ticks2(ii)) + CenterX;
    y_ticks2_in(ii,1) = radius_ticks2 * sin(angles_ticks2(ii)) + CenterY;
end

Vid = VideoWriter([dir_input,'temperature_gauge_Global']);
%Vid.VideoCompressionMethod = ‘H.294’;
Vid.FrameRate = 30;
open(Vid)
for tt = length(T)+1:length(T)+90
    % Datapoint
    if tt <= length(T)
        datapoint = T(tt);
    else
        datapoint = T(end);
    end
    
    angles_dp = 1.5*pi*(1-(datapoint-min_value)/(max_value-min_value)) - 0.25*pi;
    x_dp = radius * cos(angles_dp) + CenterX;
    y_dp = radius * sin(angles_dp) + CenterY;
    
    % Base
    angles_base_left = 1.5*pi*(1-(datapoint-min_value)/(max_value-min_value)) - 0.25*pi - pi/2;
    x_base_left = radius_base * cos(angles_base_left) + CenterX;
    y_base_left = radius_base * sin(angles_base_left) + CenterY;
    angles_base_right = 1.5*pi*(1-(datapoint-min_value)/(max_value-min_value)) - 0.25*pi + pi/2;
    x_base_right = radius_base * cos(angles_base_right) + CenterX;
    y_base_right = radius_base * sin(angles_base_right) + CenterY;
    
    figure('Position', [200, 200, 700, 700])
    subplot('Position',[0 0 1 1]);
    hold on
    
    % Panel
    patch(x2,y2, rgb('lightgray'), 'EdgeColor', rgb('lightgray'))
    
    % Out
    plot(x3, y3, 'k-', 'LineWidth', 16)
    %plot([CenterX x3(1)],[CenterY y3(1)], 'k-', 'LineWidth', 5)
    %plot([CenterX x3(end)],[CenterY y3(end)], 'k-', 'LineWidth', 5)
    
    % Stripes
    for ii = 1:450
        patch([CenterX,x_stripes(ii,:),CenterX],...
            [CenterY,y_stripes(ii,:),CenterY],...
            cbar_stripes(ii,:), 'EdgeColor', cbar_stripes(ii,:))
    end
    plot(x1, y1, 'k-', 'LineWidth', 2)
    plot([CenterX x1(1)],[CenterY y1(1)], 'k-', 'LineWidth', 5)
    plot([CenterX x1(end)],[CenterY y1(end)], 'k-', 'LineWidth', 5)
    
    % Base
    %plot([x1(1) x1(end)], [y1(1) y1(end)], 'k-', 'LineWidth', 8)
    patch(x_base,y_base, rgb('lightgray'), 'EdgeColor', rgb('lightgray'))
    plot(x_base, y_base, 'k-', 'LineWidth', 5)
    
    % Ticks2
    for ii = 1:length(ticks2)
        plot([x_ticks2_in(ii) ,x_ticks2_out(ii)],...
            [y_ticks2_in(ii) ,y_ticks2_out(ii)],'k','LineWidth',5)
    end
    
    % Ticks1
    for ii = 1:length(ticks1)
        plot([x_ticks1_in(ii) ,x_ticks1_out(ii)],...
            [y_ticks1_in(ii) ,y_ticks1_out(ii)],'k','LineWidth',5)
    end
    
    % Tick Labels
    for ii = 1:length(ticks1)
        text(x_in(ii),y_in(ii),ticks1_str{ii},'Color','k','FontSize',24,'FontWeight','Bold',...
        'FontName','Monospaced','HorizontalAlignment','Center')
    end
    
    % Datapoint
    patch([x_base_left x_base_right x_dp x_base_left],...
        [y_base_left y_base_right y_dp y_base_left],...
        rgb('lightgray'), 'EdgeColor',rgb('lightgray'))
    plot([x_base_left x_dp], [y_base_left y_dp], 'k-', 'LineWidth', 5)
    plot([x_base_right x_dp], [y_base_right y_dp], 'k-', 'LineWidth', 5)
    
    % Panel
    plot(x2, y2, 'k-', 'LineWidth', 5)
    
    % LCD
    rectangle('Position',[-0.35,-0.75,0.7,0.3],'Curvature',[0.1,0.1],...
        'EdgeColor',rgb('Black'),'FaceColor',rgb('mediumspringgreen'),...
        'LineWidth',5)
    if tt <= length(T)
        text(0,-0.53,[datestr(time(tt),'mmm'),' ',datestr(time(tt),'yyyy')],'FontSize',30,...
            'HorizontalAlignment','Center','FontWeight','Bold',...
            'FontName','Digital-7')
        text(0,-0.65,[num2str(round(T(tt),2)),'°C'],'FontSize',40,...
            'HorizontalAlignment','Center','FontWeight','Bold',...
            'FontName','Digital-7')
    else
        text(0,-0.53,[datestr(time(end),'mmm'),' ',datestr(time(end),'yyyy')],'FontSize',30,...
            'HorizontalAlignment','Center','FontWeight','Bold',...
            'FontName','Digital-7')
        text(0,-0.65,[num2str(round(T(end),2)),'°C'],'FontSize',40,...
            'HorizontalAlignment','Center','FontWeight','Bold',...
            'FontName','Digital-7')
    end
    text(0,-0.83,'Data: GISTEMP.v4 (NASA)','FontSize',16,...
        'HorizontalAlignment','Center','FontWeight','Bold',...
        'FontName','Monospaced')
    text(0,-0.88,'Viz: @afonsogneto','FontSize',16,...
        'HorizontalAlignment','Center','FontWeight','Bold',...
        'FontName','Monospaced')
    
    grid off;
    %xlabel('X', 'FontSize', 14);
    %ylabel('Y', 'FontSize', 14);
    xlim([-1.1 1.1])
    ylim([-1.1 1.1])
    set(gca,'yticklabel',[],'xticklabel',[],'TickLength',[0 0])
    
    clear *_dp datapoint
    currFrame = getframe(gcf);
    writeVideo(Vid,currFrame);
    close all
end
close(Vid)

% Static image

% Datapoint
tt = 1691;
datapoint = T(tt);

angles_dp = 1.5*pi*(1-(datapoint-min_value)/(max_value-min_value)) - 0.25*pi;
x_dp = radius * cos(angles_dp) + CenterX;
y_dp = radius * sin(angles_dp) + CenterY;

% Base
angles_base_left = 1.5*pi*(1-(datapoint-min_value)/(max_value-min_value)) - 0.25*pi - pi/2;
x_base_left = radius_base * cos(angles_base_left) + CenterX;
y_base_left = radius_base * sin(angles_base_left) + CenterY;
angles_base_right = 1.5*pi*(1-(datapoint-min_value)/(max_value-min_value)) - 0.25*pi + pi/2;
x_base_right = radius_base * cos(angles_base_right) + CenterX;
y_base_right = radius_base * sin(angles_base_right) + CenterY;

figure('Position', [200, 200, 700, 700])
subplot('Position',[0 0 1 1]);
hold on

% Panel
patch(x2,y2, rgb('lightgray'), 'EdgeColor', rgb('lightgray'))

% Out
plot(x3, y3, 'k-', 'LineWidth', 16)
%plot([CenterX x3(1)],[CenterY y3(1)], 'k-', 'LineWidth', 5)
%plot([CenterX x3(end)],[CenterY y3(end)], 'k-', 'LineWidth', 5)

% Stripes
for ii = 1:450
    patch([CenterX,x_stripes(ii,:),CenterX],...
        [CenterY,y_stripes(ii,:),CenterY],...
        cbar_stripes(ii,:), 'EdgeColor', cbar_stripes(ii,:))
end
plot(x1, y1, 'k-', 'LineWidth', 2)
plot([CenterX x1(1)],[CenterY y1(1)], 'k-', 'LineWidth', 5)
plot([CenterX x1(end)],[CenterY y1(end)], 'k-', 'LineWidth', 5)

% Base
%plot([x1(1) x1(end)], [y1(1) y1(end)], 'k-', 'LineWidth', 8)
patch(x_base,y_base, rgb('lightgray'), 'EdgeColor', rgb('lightgray'))
plot(x_base, y_base, 'k-', 'LineWidth', 5)

% Ticks2
for ii = 1:length(ticks2)
    plot([x_ticks2_in(ii) ,x_ticks2_out(ii)],...
        [y_ticks2_in(ii) ,y_ticks2_out(ii)],'k','LineWidth',5)
end

% Ticks1
for ii = 1:length(ticks1)
    plot([x_ticks1_in(ii) ,x_ticks1_out(ii)],...
        [y_ticks1_in(ii) ,y_ticks1_out(ii)],'k','LineWidth',5)
end

% Tick Labels
for ii = 1:length(ticks1)
    text(x_in(ii),y_in(ii),ticks1_str{ii},'Color','k','FontSize',24,'FontWeight','Bold',...
        'FontName','Monospaced','HorizontalAlignment','Center')
end

% Datapoint
patch([x_base_left x_base_right x_dp x_base_left],...
    [y_base_left y_base_right y_dp y_base_left],...
    rgb('lightgray'), 'EdgeColor',rgb('lightgray'))
plot([x_base_left x_dp], [y_base_left y_dp], 'k-', 'LineWidth', 5)
plot([x_base_right x_dp], [y_base_right y_dp], 'k-', 'LineWidth', 5)

% Panel
plot(x2, y2, 'k-', 'LineWidth', 5)

% LCD
rectangle('Position',[-0.35,-0.75,0.7,0.3],'Curvature',[0.1,0.1],...
    'EdgeColor',rgb('Black'),'FaceColor',rgb('mediumspringgreen'),...
    'LineWidth',5)
if tt <= length(T)
    text(0,-0.53,[datestr(time(tt),'mmm'),' ',datestr(time(tt),'yyyy')],'FontSize',30,...
        'HorizontalAlignment','Center','FontWeight','Bold',...
        'FontName','Digital-7')
    text(0,-0.65,[num2str(round(T(tt),2)),'°C'],'FontSize',40,...
        'HorizontalAlignment','Center','FontWeight','Bold',...
        'FontName','Digital-7')
else
    text(0,-0.53,[datestr(time(end),'mmm'),' ',datestr(time(end),'yyyy')],'FontSize',30,...
        'HorizontalAlignment','Center','FontWeight','Bold',...
        'FontName','Digital-7')
    text(0,-0.65,[num2str(round(T(end),2)),'°C'],'FontSize',40,...
        'HorizontalAlignment','Center','FontWeight','Bold',...
        'FontName','Digital-7')
end
text(0,-0.83,'Data: GISTEMP.v4 (NASA)','FontSize',16,...
    'HorizontalAlignment','Center','FontWeight','Bold',...
    'FontName','Monospaced')
text(0,-0.88,'Viz: @afonsogneto','FontSize',16,...
    'HorizontalAlignment','Center','FontWeight','Bold',...
    'FontName','Monospaced')

grid off;
%xlabel('X', 'FontSize', 14);
%ylabel('Y', 'FontSize', 14);
xlim([-1.1 1.1])
ylim([-1.1 1.1])
set(gca,'yticklabel',[],'xticklabel',[],'TickLength',[0 0])
saveas(gcf,[dir_input,'speedometer.png'])

%% BRA Afonso

% Input
% [min]: minimum value (scalar)
% [max]: maximum value (scalar)
% [datapoint]: datapoint value (scalar)

% Output
% Plot of gauge indicating datapoint

% Name
name = "Afonso";

% Date of Birth
MOB = "198904";

% Country of Birth
country = "BRA";

%ncdisp("gistemp250_GHCNv4.nc")

% Time 
time1 = ncread('gistemp250_GHCNv4.nc','time');
time2 = double(time1) + datenum('18000101','yyyymmdd');
time3 = string(datestr(time2,'yyyymm'));
time = time2(find(strcmp(time3,MOB)):end);

% Lat and Lon
lat1 = ncread('gistemp250_GHCNv4.nc','lat');
lon1 = ncread('gistemp250_GHCNv4.nc','lon');

country_shape = shaperead(strcat('shp/',country,'/gadm36_',country,'_0.shp'));
country_limits = country_shape.BoundingBox;

lat = lat1(lat1 > country_limits(1,2) & lat1 < country_limits(2,2));
lon = lon1(lon1 > country_limits(1,1) & lon1 < country_limits(2,1));

T_start = [find(lon1 > country_limits(1,1) & lon1 < country_limits(2,1),1,'first'),...
    find(lat1 > country_limits(1,2) & lat1 < country_limits(2,2),1,'first'),...
    find(strcmp(time3,MOB))];
T_count = [length(lon),length(lat),length(time)];

% Temperature
T1 = ncread('gistemp250_GHCNv4.nc','tempanomaly',T_start,T_count);
T2 = squeeze(nanmean(nanmean(T1,1),2));
T = repmat(T2,1,2);
clear T1 T2 T_start T_count lat lon time1 time2 time3

% Data limits
min_value = min(T);
max_value = max(T);

% Colormap
dir_colormap = '/Users/afonso/Documents/StudentSeminar 2/';
colormap_aux = readtable([dir_colormap,'w5m4.csv'], 'HeaderLines',1);
cbar_diverge2 = table2array(colormap_aux(:,3:5));
clear colormap_aux dir_colormap
cbar_stripes = NaN(360,3);
for ii = 1:3
    cbar_stripes(:,ii) = interp1((1:length(cbar_diverge2(:,ii)))',...
        cbar_diverge2(:,ii),(1:1/39.9:length(cbar_diverge2(:,ii)))');
end
cbar_stripes = flipud(cbar_stripes);
clear ii cbar_diverge2

% Info
radius = 1;
radius_base = 0.05;
CenterX = 0;
CenterY = 0;

% Gauge
angles1 = linspace(-pi/4, 5*pi/4, 360);
x1 = radius * cos(angles1) + CenterX;
y1 = radius * sin(angles1) + CenterY;

% Panel
angles2 = linspace(0, 2*pi, 360);
x2 = radius * cos(angles2) + CenterX;
y2 = radius * sin(angles2) + CenterY;

% Stripes
angles_stripes = NaN(360,3);
x_stripes = NaN(360,3);
y_stripes = NaN(360,3);
for ii = 1:360
    angles_stripes(ii,1:3) = linspace((ii-1)*pi/360 - pi/4,(ii)*pi/360 + pi/4,3);
    x_stripes(ii,1:3) = radius * cos(angles_stripes(ii,:)) + CenterX;
    y_stripes(ii,1:3) = radius * sin(angles_stripes(ii,:)) + CenterY;
end

% Base
angles_base = linspace(0, 2*pi, 360);
x_base = radius_base * cos(angles_base) + CenterX;
y_base = radius_base * sin(angles_base) + CenterY;

Vid = VideoWriter(['temperature_gauge_',name]);
%Vid.VideoCompressionMethod = ‘H.294’;
Vid.FrameRate = 50;
open(Vid)
for tt = 1:length(T)
    % Datapoint
    datapoint = T(tt);
    angles_dp = pi*(1-(datapoint-min_value)/(max_value-min_value));
    x_dp = radius * cos(angles_dp) + CenterX;
    y_dp = radius * sin(angles_dp) + CenterY;
    % Base
    angles_base_left = pi*(1-(datapoint-min_value)/(max_value-min_value)) - pi/2;
    x_base_left = radius_base * cos(angles_base_left) + CenterX;
    y_base_left = radius_base * sin(angles_base_left) + CenterY;
    angles_base_right = pi*(1-(datapoint-min_value)/(max_value-min_value)) + pi/2;
    x_base_right = radius_base * cos(angles_base_right) + CenterX;
    y_base_right = radius_base * sin(angles_base_right) + CenterY;

    figure('Position', [200, 200, 700, 700])
    subplot('Position',[0 0 1 1]);
    hold on
    
    % Panel
    patch(x2,y2, rgb('gainsboro'), 'EdgeColor', rgb('gainsboro'))
    
    for ii = 1:360
        patch([CenterX,x_stripes(ii,:),CenterX],...
            [CenterY,y_stripes(ii,:),CenterY],...
            cbar_stripes(ii,:), 'EdgeColor', cbar_stripes(ii,:))
    end
    plot(x1, y1, 'k-', 'LineWidth', 2)
    plot([CenterX x1(1)],[CenterY y1(1)], 'k-', 'LineWidth', 5)
    plot([CenterX x1(end)],[CenterY y1(end)], 'k-', 'LineWidth', 5)
    
    % Base
    %plot([x1(1) x1(end)], [y1(1) y1(end)], 'k-', 'LineWidth', 8)
    patch(x_base,y_base, 'w', 'EdgeColor', 'w')
    plot(x_base, y_base, 'k-', 'LineWidth', 5)
    
    % Datapoint
    patch([x_base_left x_base_right x_dp x_base_left],...
        [y_base_left y_base_right y_dp y_base_left],...
        'w', 'EdgeColor', 'w')
    plot([x_base_left x_dp], [y_base_left y_dp], 'k-', 'LineWidth', 5)
    plot([x_base_right x_dp], [y_base_right y_dp], 'k-', 'LineWidth', 5)
    
    % Panel
    plot(x2, y2, 'k-', 'LineWidth', 5)
    
    % LCD
    rectangle('Position',[-0.35,-0.75,0.7,0.3],'Curvature',[0.2,0.2],...
        'EdgeColor',rgb('DarkSlateGray'),'FaceColor',rgb('palegreen'),...
        'LineWidth',3)
    text(0,-0.53,[datestr(time(tt),'mmm'),' ',datestr(time(tt),'yyyy')],'FontSize',30,...
        'HorizontalAlignment','Center','FontWeight','Bold',...
        'FontName','Monaco')
    text(0,-0.67,[num2str(round(T(tt),2)),'°C'],'FontSize',36,...
        'HorizontalAlignment','Center','FontWeight','Bold',...
        'FontName','Monaco')
    
    grid off;
    xlabel('X', 'FontSize', 14);
    ylabel('Y', 'FontSize', 14);
    xlim([-1.1 1.1])
    ylim([-1.1 1.1])
    set(gca,'xticklabel',[],'TickLength',[0 0])
    
    clear *_dp datapoint
    currFrame = getframe(gcf);
    writeVideo(Vid,currFrame);
    close all
end
close(Vid)
