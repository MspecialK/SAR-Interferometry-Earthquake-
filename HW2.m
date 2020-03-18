close all;
clear
clc

%% Parameters
f=5405000000;  %Hz
c=299792458;   %m/s


%% TASK 1 -------------------------------------------------------------------------------------
% 1.1 Loading the files and subsampling
subs = 3;

infosFileMaster04 = h5info('C:\Users\mykha\Desktop\Radar Systems\H2\sardata\20190704.h5');
Data04_Re = h5read(infosFileMaster04.Filename,"/i_VV");
Data04_Im = h5read(infosFileMaster04.Filename, "/q_VV");
Distance_04 = h5read(infosFileMaster04.Filename,"/topoDistance");

Data04_Re = Data04_Re(1:subs:end, 1:subs:end);
Data04_Im = Data04_Im(1:subs:end, 1:subs:end);
Distance_04 = Distance_04(1:subs:end, 1:subs:end);

master = single(Data04_Re)+1j.*single(Data04_Im); 
master_amp = abs(master);
master_phase = angle(master);

clear Data04_Re;  % Saves memory
clear Data04_Im;  % Saves memory

%% Plotting
figure

subplot(1,3,1)
imagesc(master_amp);                      % Amplitude before removing sparks
title('Original data Amplitude');

subplot(1,3,3)
imagesc(master_phase);                    % Phase plot
title('Original data phase');

%%Sparks removal
lim=quantile(master_amp,0.98,'all');      % Limiting value to .98 percentile
master_amp(master_amp>=lim)=lim;          % Limiting all values above limit to the limit
clear lim;  % Saves memory

subplot(1,3,2)
imagesc(master_amp);                      % Amplitude limited
title('Spark-limited to 98 percentile');

%% Averaging 
a=ones(5)/(5^2);                          % 5x5 Window
b=ones(10)/(10^2);                        % 10x10 Window

master_abs_1 = conv2(master_amp,a,'same');
master_abs_2 = conv2(master_amp,b,'same');

figure  % We are going to plot different filter size effects on the Amplitude 

% No Filter aplied plot and histogram of the values
subplot(2,3,1)
imagesc(master_amp);
title ('Amplitude');       
subplot(2,3,4)
histogram(master_amp);         
title ('Amplitude - Histogram');

% 5x5 Filter Applied
subplot(2,3,2)
imagesc(master_abs_1);         
title ('Amplitude 5');
subplot(2,3,5)
histogram(master_abs_1);       
title ('Amplitude 5 - Histogram');

% 10x10 Filter Applied
subplot(2,3,3)
imagesc(master_abs_2);
title ('Amplitude 10');
subplot(2,3,6)
histogram(master_abs_2);
title ('Amplitude 10 - Histogram');


%% SAVE for QGIS
load('C:\Users\mykha\Desktop\Radar Systems\H2\sardata\geocoding_infos.mat'); 
geocodingInfos = GecSubs(geocodingInfos, [subs subs]); 
geocodedImage = geocodeRadarImage(master_abs_1, geocodingInfos);
geotiffwrite("C:\Users\mykha\Desktop\ima.tif",geocodedImage,geocodingInfos.xref);

clear master_abs_1;  % Saves memory
clear master_abs_2;  % Saves memory

%% TASK 1.2 -------------------------------------------------------------------------------------

infosFileMaster16 = h5info('C:\Users\mykha\Desktop\Radar Systems\H2\sardata\20190716.h5');
Data16_Re = h5read(infosFileMaster16.Filename,"/i_VV");
Data16_Im = h5read(infosFileMaster16.Filename, "/q_VV");
Distance_16 = h5read(infosFileMaster16.Filename,"/topoDistance");

Data16_Re = Data16_Re(1:subs:end, 1:subs:end);
Data16_Im = Data16_Im(1:subs:end, 1:subs:end);
Distance_16 = Distance_16(1:subs:end, 1:subs:end);

slave = single(Data16_Re)+1j.*single(Data16_Im);
clear Data16_Re;  % Saves memory
clear Data16_Im;  % Saves memory

interferogram = master.*conj(slave);
figure
subplot(1,2,1)
imagesc(angle(interferogram));
colorbar

slave_d = slave.*exp(1j.*(4*pi*f/c)*Distance_16);
master_d = master.*exp(1j.*(4*pi*f/c)*Distance_04);

clear Distance_04;  % Saves memory
clear Distance_16;  % Saves memory

interferogram_d=master_d.*conj(slave_d);
subplot(1,2,2)
imagesc(angle(interferogram_d));
colorbar

%%Geocoding
geocodedImage = geocodeRadarImage(angle(interferogram_d), geocodingInfos);
geotiffwrite("C:\Users\mykha\Desktop\ima_t2.tif",geocodedImage,geocodingInfos.xref);

%% TASK 1.3

M=ones(43,11)*(1/473);
interferogram_clean=conv2(angle(interferogram_d),M,'same');

figure 
subplot(1,2,1)
imagesc(angle(interferogram_d));
colorbar

subplot(1,2,2)
imagesc(interferogram_clean);
colorbar

geocodedImage = geocodeRadarImage(interferogram_clean, geocodingInfos);
geotiffwrite("C:\Users\mykha\Desktop\ima_bello.tif",geocodedImage,geocodingInfos.xref);

%% TASK 2.1 --------------------------------------------------------------------------------------

exy=conv2(interferogram_d,a,'same');
ex1=conv2((abs(master_d).^2),a,'same');
ex2=conv2((abs(slave_d).^2),a,'same');
coherence=exy./sqrt(ex1.*ex2);

clear exy;  % Saves memory
clear ex1;  % Saves memory
clear ex2;  % Saves memory


figure
subplot(1,2,1)
imagesc(angle(coherence));
colorbar

subplot(1,2,2)
imagesc(abs(coherence));
colorbar

geocodedImage = geocodeRadarImage(abs(coherence), geocodingInfos);
geotiffwrite("C:\Users\mykha\Desktop\ima_coherence.tif",geocodedImage,geocodingInfos.xref);

%% TASK 2.2 ------------------------------------------------------------------------------------

coh_sub=coherence(1:20:end, 1:5:end); 
interferogram_unw = unwrap_IRLS(double(angle(coh_sub)), abs(coh_sub), 100, [], 1);

figure
imagesc(interferogram_unw);
colorbar

%%Interpolation to save the right image

newNumberOfRows = 8425; 
newNumberOfCols = 2292;  
[x, y] = meshgrid(1:size(interferogram_unw,2), 1:size(interferogram_unw,1));
[xq, yq] = meshgrid(linspace(1, size(interferogram_unw, 2), newNumberOfCols), linspace(1, size(interferogram_unw, 1), newNumberOfRows));
Unwraped_phase = interp2(x, y, interferogram_unw, xq, yq);

% Saving the image
geocodedImage = geocodeRadarImage(Unwraped_phase, geocodingInfos);
geotiffwrite("C:\Users\mykha\Desktop\ima_coherence_UNW.tif",geocodedImage,geocodingInfos.xref);

%% TASK 3.1 -----------------------------------------------------------------------------------------------

[FX,FY] = gradient(Unwraped_phase); 
figure
subplot(1,2,1)
imagesc(FX);
title('Gradient along the X axis');

subplot(1,2,2)
imagesc(FY);
title('Gradient along the Y axis');

geocodedImage = geocodeRadarImage(FX, geocodingInfos);
geotiffwrite("C:\Users\mykha\Desktop\ima_XX_UNW.tif",geocodedImage,geocodingInfos.xref);

geocodedImage = geocodeRadarImage(FY, geocodingInfos);
geotiffwrite("C:\Users\mykha\Desktop\ima_YY_UNW.tif",geocodedImage,geocodingInfos.xref);

%% TASK 3.2 ----------------------------------------------------------------------------------------------

line=movmean(FX(4000,:),100);

figure
subplot(2,1,2)
plot(line);
grid on;
title('Gradient along X-axe');

subplot(2,1,1)
plot(Unwraped_phase(4000,:));
grid on;
title('Unwraped phase along X-axe');

lambda=c/f;
finale=(lambda/(4*pi).*Unwraped_phase); %displacement in centimeters
figure
contour(finale,'ShowText','on');
title('Displacement in m');

geocodedImage = geocodeRadarImage(finale, geocodingInfos);
geotiffwrite("C:\Users\mykha\Desktop\ima_cont.tif",geocodedImage,geocodingInfos.xref);

% Instead of using improfile that is boring and just some lines we can plot
% the FX and FY complete with a little adjustment by avareging. 
z=ones(100)/10000;
FXXX=conv2(FX,z,'same');
FXXX=FXXX(1:10:end,1:10:end);
figure;
surf(FXXX);




