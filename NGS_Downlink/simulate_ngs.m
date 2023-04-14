clear all();

%% Parameters
nRes = 256; % [pix/m] 最終的な分解能 
elevation = 55; % [degree] Elevation
zenithAngle = 90 - elevation; % [degree] zenith angle 
frameRate = 400; % [Hz] シミュレーションのフレームレート
duration = 10; % [seconds] 計算する大気揺らぎの時間方向の長さ
frameNumber = frameRate * duration; % 計算するフレーム数
atmModel = 'Oya2014';
%atmProfile = 'good';
%atmProfile = 'median';
atmProfile = 'bad';

%% 3 graphs (propagation, cam, wfs_image)
fig = figure(3);
fig.Position

%% telescope class
tel = telescope(...
    8,... % telescope diameter
    'resolution', nRes,... % telescope resolution
    'obstructionRatio', 0.2777,...
    'fieldOfViewInArcsec', 1,...
    'samplingTime', 1/frameRate);

%% Atmosphere
altitude = [0,0.5,1,2,4,8,16]*1e3; % [m] 大気層の高さ
%r0 = 0.1813; % [m] Fried parameter 589nm good
%r0 = 0.1436; % [m] Fried parameter 589nm median
r0 = 0.1107; % [m] Fried parameter 589nm bad
%r0 = 0.1180; % [m] Fried parameter 500nm
L0 = 30; % [m] Outer scale
%fractionalR0 = [0.7546,0.0497,0.0141,0.0133,0.0545,0.0467,0.0671]; % 大気揺らぎの強度比 good
%fractionalR0 = [0.7316,0.0650,0.0193,0.0252,0.0574,0.0500,0.0515]; % 大気揺らぎの強度比 median
fractionalR0 = [0.6882,0.0798,0.0398,0.0395,0.0551,0.0548,0.0428]; % 大気揺らぎの強度比 bad
windSpeed = [5.6,5.77,6.25,7.57,13.31,19.06,12.14]; % [m/s] 風速
windDirection = [0,50,100,150,200,250,300,350]*pi/180; % [radian] 風向​

% Elevationによるスケール
altitude = altitude/cosd(zenithAngle);
r0 = r0*cosd(zenithAngle)^(3/5);

% generate atmosphere class
atm = atmosphere(photometry.Na,...
    r0,...
    L0,...
    'fractionnalR0',fractionalR0,...
    'altitude',altitude,...
    'windSpeed',windSpeed,...
    'windDirection',windDirection,...
    'randStream',RandStream.create('mt19937ar','Seed', 'shuffle')); % 大気揺らぎ生成の乱数シードの変更
%atm.wavelength = photometry.Na; % LGS (NGSの場合は無視)

%% source class
% 4 source, 40 arcsec off-axis
%ast = {[1,arcsec(40),0]}; % 40arcsecで1つのLGS
%src = source(...
%    'asterism', ast,...
%    'wavelength', photometry.Na);
src = source('wavelength', photometry.Na); % single on-axis NGS

%src.nPhoton = 8.043*10^6; % LGS, EL90 40arcsec max
%src.nPhoton = 3.973*10^6; % LGS, EL90 40arcsec min

%% Imager
cam = imager(...
    'fieldStopSize', 128,...
    'nyquistSampling', 1,...
    'exposureTime', 1); % 400 Hz(tel) * 1 time

%% Propagation
src = src.*tel*cam; % propagation
%ax1 = subplot(1,3,1);
%imagesc(cam); % show imager image
%imagesc(src.catMeanRmOpd()); % tip/tilt from wavefront

%% Shack Hartmann
% 32x32 sub-aperture
% intensity threshold 0.6
wfs = shackHartmann(16*2, 256*4, 0.6); % tmp!!!!!!!!!!!!!!!
wfs.lenslets.nyquistSampling = 2; % FWHMを4pixelでサンプリング
wfs.lenslets.fieldStopSize = 8; % 256/32=8で最大の視野を使用
%wfs.lenslets.throughput = 0.45; % overall system throughput
wfs.camera.exposureTime = 1; % 400 Hz(tel) * 1 time
%wfs.camera.photonNoise = true;
wfs.camera.nPhotonBackground = 0; % tmp!!!!!!!!!!!!!!!月次第(レーザー想定で0)
%wfs.camera.readOutNoise = 1.6; % tmp!!!!!!!!!!!!!!!
%wfs.framePixelThreshold = 5e4; % 高周波数に効く

src = src.*tel*wfs; % first propagation
wfs.INIT() % 60%以上照らされているsubapertureだけで重心計算

%% Shack Hartmann (1x1)
% intensity threshold 0.6
wfs1 = shackHartmann(1, 64, 1); % tmp!!!!!!!!!!!!!!!
wfs1.lenslets.nyquistSampling = 1; % FWHMを4pixelでサンプリング
wfs1.lenslets.fieldStopSize = 256; % 256/1=256で最大の視野を使用
%wfs1.lenslets.throughput = 0.45; % overall system throughput
wfs1.camera.exposureTime = 1; % 400 Hz(tel) * 1 time
%wfs1.camera.photonNoise = true;
wfs1.camera.nPhotonBackground = 0; % tmp!!!!!!!!!!!!!!!月次第(レーザー想定で0)
%wfs1.camera.readOutNoise = 1.6; % tmp!!!!!!!!!!!!!!!
%wfs1.framePixelThreshold = 1e5; % 高周波数に効く

src = src.*tel*wfs1; % first propagation
wfs1.INIT() % initialize

%% plot wfs image
ax2 = subplot(1,3,2);
imagesc(wfs1.camera); % show wfs image
ax3 = subplot(1,3,3);
imagesc(wfs.camera);

%% 配列の確保
slopes_ave = zeros(frameNumber, 2);
slopes_ave1wfs = zeros(frameNumber, 2);
wave_tt = zeros(frameNumber, 2);

%% zernike class for tip/tilt decoposition

% generate zernike class for 2,3 modes (tip/tilt)
zer = zernike(2:3,'resolution',tel.resolution(1),'pupil',tel.pupilLogical);

% get zernike wavefront within telescope pupil
zer = zer.modes(tel.pupilLogical(:),:)/4; % PV = 1m

% invert zernike mode
izer = pinv(zer);

%% generate atmosphere layers
tel = tel + atm;
%ax1 = subplot(2,3,[1,2,3]);
%imagesc(tel,src);
ax1 = subplot(1,3,1);
zer_image = src.catMeanRmOpd();
imagesc(zer_image)

%% time evolution
for k = 1:frameNumber
    fprintf('%d\n',k);
    +tel;
    src = src.*tel*cam*wfs*wfs1;

    %% tip/tilt from wavefront in arcsec
    wave_tt(k,:)=(izer*src.catMeanRmOpd(tel.pupilLogical)/tel.D/arcsec(1))';

    %% tip/tilt from 1x1 SH-WFS in arcsec
    slopes_ave1wfs(k,:) = wfs1.slopes'*589e-9/tel.D/2*8/arcsec(1);
    
    %% xslope_ave, yslope_ave
    slopes = wfs.slopes;
    %disp(slopes) % (32*32)*2 = 2048*1
    % INIT() -> 1488
    xslopes = slopes(1:744);
    yslopes = slopes(745:1488);
    xslope_ave = mean(xslopes, 'omitnan'); % ignore NaN
    yslope_ave = mean(yslopes, 'omitnan'); % ignore NaN
    disp(xslope_ave)
    disp(yslope_ave)
    slopes_ave(k,1) = xslope_ave;
    slopes_ave(k,2) = yslope_ave;
    
    %% view
    %src = src.*tel*cam;
    %imagesc(tel); % propagation
    %imagesc(cam); % cam
    imagesc(src.catMeanRmOpd()) % tip/tilt from wavefront
    %src = src.*tel*wfs*wfs1;
    imagesc(wfs1.camera); % show wfs1 image
    imagesc(wfs.camera); % show wfs image
    %imagesc(slopes)
    drawnow
    
end

%% Save
%save('slopes_ave_update20230305_bad.mat','slopes_ave');
%save('slopes_ave1wfs_update20230305_bad.mat','slopes_ave1wfs');
%save('wave_tt_update20230305_bad.mat','wave_tt');