%------------------------------------------------------------------------%
% MMS instrument information: 
% https://www.nasa.gov/centers/ames/earthscience/programs/MMS/instruments
% List of MMS missions: 
% https://airbornescience.nasa.gov/mms
% Flight data: 
% https://espoarchive.nasa.gov/archive/browse
%
% This script can be used to scan through a list of 1Hz data files (the
% 1Hz files are much faster to download than the 20Hz files) to
% double-check altitude profiles and flight tracks/locations. The file
% name template and mission date list will have to be changed for each
% different mission.
%
%------------------------------------------------------------------------%
clear all;
close all;

missions = ["0913","0915","1001","1002","1003","1008","1012","1014","1015","1018","1019","1021","1024","1025","1026","1028","1030","1031","1101","1102"];  %CHANGE FOR EACH MISSION
plot_title = 'WB-57 Flight';   %changes title for all plots
origin = 'Guam';   % operating base for the mission
utc_to_local = 10;   % (hrs) convert from UTC time to local time to see what time of day aircraft is ascending and descending (+ is ahead of UTC time, - is behind UTC time)

for i = 1:length(missions)
    filename = strcat('MMS-1HZ_WB57_2016',num2str(missions(i)),'.ict');   %CHANGE FOR EACH MISSION

    disp('Loading data...');
    delimiterIn = ',';
    headerlinesIn = 65;   %check this number in first row of data file
    A = importdata(filename,delimiterIn,headerlinesIn);

    W = A.data(:,7)/(1e3);     
    index = find(W < -50);    %find where vertical wind speed data starts recording
    
    A.data(index,:) = [];

    time = A.data(:,1);     %time from UTC midnight (s)
    temp = A.data(:,3)/(1e3);     %static temperature (K)
    lat = A.data(:,10)/(1e5);     %latitude from LN100 INS (deg +N)
    lon = A.data(:,11)/(1e5);     %longitude from LN100 INS (deg +E)
    altitude = A.data(:,12)/(1e1);     %altitude from LN100 INS (m)
    
    
    i_lat = find(lat<-90);   %get rid of points where lat/lon was not recorded
    lat(i_lat) = [];
    i_lon = find(lon<-900);
    lon(i_lon) = [];
    
    %Temperature profile for tropopause estimation
%     figure
%     plot(temp,altitude/1000)
%     xlabel('Static Temperature (K)');ylabel('Altitude (km)');
%     title(plot_title);
%     grid on;
    
    %Flight track
%     figure
%     geoplot(lat,lon,'LineWidth',1.5);
%     geobasemap bluegreen
%     title(plot_title);
    
    %Altitude profile
%     figure
%     plot(time/3600+utc_to_local,altitude/1000)
%     xlabel([origin, ' Local Time']);ylabel('Altitude (km)');
%     title(plot_title);
%     grid on;
    
end
