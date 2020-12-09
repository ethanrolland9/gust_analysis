function [max_cos_amp,max_hcos_amp,cos_amplitudes,hcos_amplitudes] = mms_gusts(date,phase)
%------------------------------------------------------------------------%
% MMS instrument information: 
% https://www.nasa.gov/centers/ames/earthscience/programs/MMS/instruments
% List of MMS missions: 
% https://airbornescience.nasa.gov/mms
% Flight data: 
% https://espoarchive.nasa.gov/archive/browse
%
% Want to look for data that was recorded by high-altitude platforms
% (ER-2/WB-57/Global Hawk) at 20Hz. Missions before the early 2000's
% were only recorded at 5Hz. Almost all missions have MMS data that is
% low-passed to 1Hz, which is useful for quickly checking altitudes
% locations (see mms_data_scan.m script). For some missions, the 20Hz 
% data is only available upon request from the NASA POC (T. Paul Bui: 
% thaopaul.v.bui@nasa.gov)
%
% INPUTS:
% date = character array corresponding to flight date (e.g. '1021' for 
%        October 21st. *Have to change line 28 for each mission*
% phase = optional string specifying desired flight phase (default is
%         entire flight)
%         "ASCENT" = takeoff through 13000 m altitude
%         "DESCENT" = 13000 m altitude through landing
%         "STRATOSPHERE" = altitude is above specified tropopause
%         "TROPOSPHERE" = altitude between 13000 m and tropopause
% OUTPUTS:
% max_cos_amp = maximum cosine gust amplitude
% max_hcos_amp = maximun half-cosine gust amplitude
% cos_amplitudes = all cosine gust amplitudes during flight
% hcos_amplitudes = all half-cosine gust amplitudes during flight
%
%------------------------------------------------------------------------%
%% Variables
% clear all;
% close all;

filename = ['MMS-20HZ_WB57_2016',date,'.ict'];
data_freq = 20;       % (Hz) MMS sampling frequency, if changed make sure to double-check variable assignments in the next section
interval = 1;         % (m) Interval used to sample data
wingspan = 40;        % (m) Dawn aircraft wingspan
tropopause = 17000;   % (m) Estimated tropopause altitude, determined by looking at temperature profile
TEMP_PLOT_FLAG = 0;   % Temperature profile plot ON (= 1) or OFF (= 0)
EAS = 10;             % (m/s) Equivalent (sea-level) airspeed of Dawn aircraft
origin = 'Guam';      % Operating base for the mission
utc_to_local = 10;    % (hrs) Convert from UTC time to local time to determine time of day for ascent and descent (+ is ahead of UTC time, - is behind UTC time)
SHADING_FLAG = 0;     % Shading of stratosphere in plots ON (= 1) or OFF (= 0)
MARKER_FLAG = 0;      % Interval markers in plots ON (= 1) or OFF (= 0)
GUST_FLAG = 1;        % Solve for worst-case cosine gust ON (= 1) or OFF (= 0)
GUST_PLOT_FLAG = 0; % Wingtip and fuselage markers in plots ON (= 1) or OFF (= 0)
plot_title = ['WB-57 Flight: ',date(1:2),'/',date(3:4),'/2016'];   %Title for all plots

%% Load Data
% disp('Loading data...');
delimiterIn = ',';
headerlinesIn = 70;   %check this number in first row of data file
A = importdata(filename,delimiterIn,headerlinesIn);

W = A.data(:,7)/(1e3);     
index = find(W < -50);    %find where vertical wind speed data starts recording
A.data(index,:) = [];

%filter data based on specified phase of flight
if (~exist('phase','var'))
    phase = "ALL";
end
if phase == "ASCENT"
    altitude = A.data(:,12)/(1e1);
    index = find(altitude > 13000);
    A.data(index,:) = [];
    altitude = A.data(:,12)/(1e1);
    [~,I] = max(altitude);
    A.data(I:end,:) = [];
elseif phase == "DESCENT"
    altitude = A.data(:,12)/(1e1);
    index = find(altitude > 13000);
    A.data(index,:) = [];
    altitude = A.data(:,12)/(1e1);
    [~,I] = max(altitude);
    A.data(1:I,:) = [];
elseif phase == "STRATOSPHERE"
    altitude = A.data(:,12)/(1e1);
    if max(altitude)>tropopause
        index = find(altitude < tropopause);
        A.data(index,:) = [];
    else
        max_cos_amp = NaN;
        max_hcos_amp = NaN;
        return
    end
elseif phase == "TROPOSPHERE"
    altitude = A.data(:,12)/(1e1);
    index1 = altitude > tropopause;
    index2 = altitude < 13000;
    total = index1+index2;
    index = find(total>0);
    A.data(index,:) = [];
end

%load data
time = A.data(:,1);     %time from UTC midnight (s)
temp = A.data(:,3)/(1e3);     %static temperature (K)
TAS = A.data(:,4)/(1e2);     %true airspeed (m/s)
W = A.data(:,7)/(1e3);     %vertical wind speed (m/s)
lat = A.data(:,10)/(1e5);     %latitude from LN100 INS (deg +N)
lon = A.data(:,11)/(1e5);     %longitude from LN100 INS (deg +E)
altitude = A.data(:,12)/(1e1);     %altitude from LN100 INS (m)

%% Find Distance Traveled by Aircraft
% disp('Converting to distance traveled...');
freq = 1/data_freq;
distance = zeros(length(time),1);
distance(1) = time(1)*TAS(1);
for i = 2:length(time)
   distance(i) = distance(i-1) + (TAS(i)*freq);
end
distance = distance - distance(1);


%% Find Delta AoA Perturbations
TAS_dawn = zeros(length(TAS),1);
[~,~,~,rho] = atmoscoesa(altitude);
for i = 1:length(altitude)      %convert Dawn EAS to TAS based on density profile
    TAS_dawn(i) = EAS/sqrt(rho(i)/1.225);     %https://en.wikipedia.org/wiki/Equivalent_airspeed
end
delta_alpha = W./TAS_dawn*(180/pi);     %delta AoA = vertical gust/aircraft cruise velocity

%% Interpolate Data at Desired Intervals
% disp(['Interpolating data at ',num2str(interval),' m intervals...']);
max_dist = max(distance);
array_length = floor(max_dist/interval) + 10;
markers = zeros(1,array_length);
max_dist = max(distance);
for i = 2:array_length
   dist = interval*i;
   if dist <= max_dist
       markers(i) = dist;
   end
end
markers = nonzeros(markers);
w_sampled = interp1(distance,W,markers,'linear');    %vertical wind velocity
alt_sampled = interp1(distance,altitude,markers,'linear');     %altitude
alpha_sampled = interp1(distance,delta_alpha,markers,'linear');   %delta AoA perturbations

%% Find Tropopause
if TEMP_PLOT_FLAG == 1
    figure;
    plot(temp,altitude/1000)
    xlabel('Static Temperature (K)');ylabel('Altitude (km)');
    title(plot_title);
    grid on;
end

%% Stratosphere Shading for Plots
if SHADING_FLAG == 1
    filter = altitude>tropopause;
    indices = find(filter==1);
    if isempty(indices)
        warning('Aircraft does not reach stratosphere, turn shading off or adjust tropopause height')
    end
    points = [0;0];
    block_start = indices(1);
    for i = 2:length(indices)
        if indices(i) ~= indices(i-1)+1
           block_end = indices(i-1);
           point = [block_start;block_end];
           points = [points,point];
           block_start = indices(i);
        end
    end
    point = [block_start;indices(i)];
    points = [points,point];
    points = points(:,2:end);
end

%% Solve for Worst-case Gust
if GUST_FLAG == 1
%     disp('Solving for worst-case gusts...')
    max_cos_amp = 0;
    max_amplitude_start_index = 0;
    max_hcos_amp = 0;
    max_deflect_start_index = 0;
    index_gap = wingspan/interval;
    if rem(index_gap,2) ~= 0
        warning('Choose interval that divides evenly into aircraft wingspan')
    end
    cos_amplitudes = zeros(length(markers)-index_gap,1);
    hcos_amplitudes = zeros(length(markers)-index_gap,1);
    for index = 1:length(markers)-index_gap
        left_wing_gust = w_sampled(index);     %windspeed at left wingtip
        right_wing_gust = w_sampled(index+index_gap);     %windspeed at right wingtip
        fuselage_gust = w_sampled(index+(index_gap/2));     %windspeed at fuselage
        
        left_halfspan_gust = fuselage_gust-left_wing_gust;    %differential between left wingtip and fuselage
        right_halfspan_gust = right_wing_gust-fuselage_gust;     %differential between right wingtip and fuselage
        
        %cosine gust
        if left_halfspan_gust*right_halfspan_gust < 0  %check that wingtip gusts are opposite sign from fuselage gust (cosine condition)
           c_amplitude = (abs(left_halfspan_gust) + abs(right_halfspan_gust))/4;   %amplitude defined as average difference over half-spans
           if c_amplitude > max_cos_amp
              max_cos_amp = c_amplitude;
              max_amplitude_start_index = index;
           end
           cos_amplitudes(index) = c_amplitude;
        end
        
        %span-wise gust
        hc_amplitude = abs(left_wing_gust - right_wing_gust)/2;    %amplitude defined as half the difference between wingtips
        if hc_amplitude > max_hcos_amp
           max_hcos_amp = hc_amplitude;
           max_deflect_start_index = index;
        end
        hcos_amplitudes(index) = hc_amplitude;
        
    end
    
%     disp(['Worst-case cosine gust amplitude: ',num2str(max_amplitude),' m/s'])
%     disp(['Distance location of worst-case amplitude: ',num2str(markers(max_amplitude_start_index)),' m'])
%     disp(['Worst-case span-wise gust difference: ',num2str(max_deflect),' m/s'])
%     disp(['Distance location of worst-case difference: ',num2str(markers(max_deflect_start_index)),' m'])
end

%points for plotting
cos_pts = [markers(max_amplitude_start_index) w_sampled(max_amplitude_start_index);
              markers(max_amplitude_start_index+index_gap/2) w_sampled(max_amplitude_start_index+index_gap/2);
              markers(max_amplitude_start_index+index_gap) w_sampled(max_amplitude_start_index+index_gap)];
hcos_pts = [markers(max_deflect_start_index) w_sampled(max_deflect_start_index);
               markers(max_deflect_start_index+index_gap),w_sampled(max_deflect_start_index+index_gap)];

%% Plotting
% disp('Plotting...');

%FLIGHT TRACK
% figure
% geoplot(lat,lon,'LineWidth',1.5)
% geobasemap bluegreen
% title(plot_title)

%ALTITUDE AND VERTICAL WIND SPEEDS VS. TIME
% figure;
% yyaxis right
% plot(time/3600+utc_to_local,altitude/1000)
% ylabel('Altitude (km)');
% yl = ylim;
% if SHADING_FLAG == 1
%     for i = 1:length(points)
%         patch([time(points(1,i))/3600 + utc_to_local time(points(2,i))/3600 + utc_to_local time(points(2,i))/3600 + utc_to_local time(points(1,i))/3600 + utc_to_local],[yl(1) yl(1) yl(2) yl(2)],'green','EdgeColor','none','FaceAlpha',0.2)
%     end
% end
% yyaxis left
% plot(time/3600+utc_to_local,W)
% xlabel([origin, ' Local Time']);ylabel('Vertical Wind Speed (m/s)');
% title(plot_title);
% grid on;

%ALTITUDE AND VERTICAL WIND SPEEDS VS. DISTANCE
% figure
% yyaxis right
% hold on
% plot(distance/1000,altitude/1000,'LineWidth',1.5,'HandleVisibility','off')
% grid on
% xlabel('Distance (km)');ylabel('Altitude (km)');
% yl = ylim;
% if SHADING_FLAG == 1
%     patch([distance(points(1,1)) distance(points(2,1)) distance(points(2,1)) distance(points(1,1))],[yl(1) yl(1) yl(2) yl(2)],'green','EdgeColor','none','FaceAlpha',0.2,'DisplayName','Aircraft in Stratosphere')
%     if length(points(1,:)) > 1
%         for i = 2:length(points(1,:))
%             patch([distance(points(1,i)) distance(points(2,i)) distance(points(2,i)) distance(points(1,i))],[yl(1) yl(1) yl(2) yl(2)],'green','EdgeColor','none','FaceAlpha',0.2,'HandleVisibility','off')
%         end
%     end
% end
% ylim(yl);
% yyaxis left
% hold on
% plot(distance/1000,W,'LineWidth',1.5,'HandleVisibility','off')
% if MARKER_FLAG == 1
%     scatter(markers,w_sampled,150,[0, 0.4470, 0.7410],'.','HandleVisibility','off');
% end
% if GUST_FLAG == 1 && GUST_PLOT_FLAG == 1
%     plot(cos_pts(:,1),cos_pts(:,2),':rs','LineWidth',1.5,'MarkerSize',10,'DisplayName','Worst Cosine Gust')
%     plot(hcos_pts(:,1),hcos_pts(:,2),'--go','LineWidth',1.5,'DisplayName','Worst Span-wise Gust Difference ')
%     legend
% end
% ylabel('Vertical Wind Speed (m/s)');
% title(plot_title);

%ALTITUDE AND DELTA AOA VS. DISTANCE
% figure
% yyaxis right
% hold on
% plot(distance,altitude/1000)
% grid on
% xlabel('Distance (m)');ylabel('Altitude (km)');
% yl = ylim;
% if SHADING_FLAG == 1
%     for i = 1:length(points)
%         patch([distance(points(1,i)) distance(points(2,i)) distance(points(2,i)) distance(points(1,i))],[yl(1) yl(1) yl(2) yl(2)],'green','EdgeColor','none','FaceAlpha',0.2)
%     end
% end
% yyaxis left
% hold on
% plot(distance,delta_alpha)
% if MARKER_FLAG == 1
%     scatter(markers,alpha_sampled,150,[0, 0.4470, 0.7410],'.','HandleVisibility','off');
% end
% ylabel('\Delta\alpha (deg)');
% title(plot_title);

%APPROXIMATE DERIVATIVE OF VERTICAL GUSTS OVER TIME
% figure
% yyaxis right
% hold on
% plot(time/3600 + utc_to_local,altitude/1000)
% grid on
% xlabel([origin, ' Local Time']);ylabel('Altitude (km)');
% yl = ylim;
% if SHADING_FLAG == 1
%     for i = 1:length(points)
%         patch([time(points(1,i))/3600 + utc_to_local time(points(2,i))/3600 + utc_to_local time(points(2,i))/3600 + utc_to_local time(points(1,i))/3600 + utc_to_local],[yl(1) yl(1) yl(2) yl(2)],'green','EdgeColor','none','FaceAlpha',0.2)
%     end
% end
% yyaxis left
% hold on
% plot(time(2:end)/3600 + utc_to_local,diff(W))
% ylabel('Difference Between Data Points (m/s)');
% grid on
% title(plot_title);
end