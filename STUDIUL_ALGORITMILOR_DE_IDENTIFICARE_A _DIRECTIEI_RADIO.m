clc 
clear
close all
%% 
raw = importdata('frecventa4_6694.txt');   
data = raw.data;                         

theta = data(:,1);   
phi   = data(:,2);  
gain_dB = data(:,3); 
u_theta = unique(theta);   
u_phi   = unique(phi);     
[PhiGrid, ThetaGrid] = meshgrid(u_phi, u_theta);   
GainGrid = reshape(gain_dB, numel(u_phi), numel(u_theta));  
GainGrid = GainGrid.';                                      
az = u_phi;             
az(az > 180) = az(az > 180) - 360;  
el = 90 - u_theta;     
[az_sorted, idx_az] = sort(az, 'ascend');
magPattern = GainGrid(:, idx_az);   
[el_sorted, idx_el] = sort(el, 'ascend');
magPattern = magPattern(idx_el, :); 
nan_mask = isnan(magPattern);
if any(nan_mask, 'all')
    minVal = min(magPattern(~nan_mask));
    magPattern(nan_mask) = minVal;
end

phasePattern = zeros(size(magPattern)); 
fc_cst = 4.6694e9;    
freqvec = [fc_cst-1e6  fc_cst+1e6];
size(magPattern) ;   
size(phasePattern);  
magPattern3D   = cat(3, magPattern, magPattern);      
phasePattern3D = cat(3, phasePattern, phasePattern);  
elem = phased.CustomAntennaElement( ...
    'AzimuthAngles',    az_sorted, ...
    'ElevationAngles',  el_sorted, ...
    'MagnitudePattern', magPattern3D, ...
    'PhasePattern',     phasePattern3D, ...
    'FrequencyVector',  freqvec);

%% ULA
ula = phased.ULA('NumElements',8,'ElementSpacing',1.5, 'Element',elem);
ang1 = [-45; 0];         
ang2 = [45; 0];       
angs = [ang1 ang2];

c = physconst('LightSpeed');
fc = 4.6694e9;              
lambda = c/fc_cst;
pos = getElementPosition(ula)/lambda;
Nsamp = 1000;
nPower = 0.01;
rs = rng(2007);
signal = sensorsig(pos,Nsamp,angs,nPower);

%% Beamscan DOA Estimation - ULA
broadsideAngle = az2broadside(angs(1,:),angs(2,:));
spatialspectrum = phased.BeamscanEstimator('SensorArray',ula,...
            'OperatingFrequency',fc_cst,'ScanAngles',-90:90);
spatialspectrum.DOAOutputPort = true;
spatialspectrum.NumSignals = 2;
[~,ang] = spatialspectrum(signal);

%% MVDR and MUSIC-ULA
ang1 = [-45; 0];     
ang2 = [45; 0];        
signal = sensorsig(pos,Nsamp,[ang1 ang2],nPower);
[~,ang] = spatialspectrum(signal);

figure(1)
plotSpectrum(spatialspectrum);
title("ULA-Beamscan")
ax = gca; 
yl = ax.YLim;
ax.YLim = [yl(1), yl(2) + 0.1 * abs(yl(2))];  
mvdrspatialspect = phased.MVDREstimator('SensorArray',ula,...
        'OperatingFrequency',fc_cst,'ScanAngles',-90:90,...
        'DOAOutputPort',true,'NumSignals',2);
[~,ang] = mvdrspatialspect(signal);

figure(2)
plotSpectrum(mvdrspatialspect);
title("ULA-MVDR")
ax = gca;  
yl = ax.YLim;
ax.YLim = [yl(1), yl(2) + 0.1 * abs(yl(2))]; 

musicspatialspect = phased.MUSICEstimator('SensorArray',ula,...
        'OperatingFrequency',fc_cst,'ScanAngles',-90:90,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',2);
[~,ang] = musicspatialspect(signal);

figure(3)
ymvdr = mvdrspatialspect(signal);
ymusic = musicspatialspect(signal);
helperPlotDOASpectra(mvdrspatialspect.ScanAngles,...
  musicspatialspect.ScanAngles,ymvdr,ymusic,'ULA')
title("ULA-MVDR and MUSIC")

figure(4)
plotSpectrum(musicspatialspect)
title("ULA-MUSIC")
ax = gca;  
yl = ax.YLim;
ax.YLim = [yl(1), yl(2) + 0.1 * abs(yl(2))];  

ang1 = [-40; 35]; ang2 = [40; 35]; ang3 = [0;0];
signal = sensorsig(pos,Nsamp,[ang1 ang2 ang3],nPower);
[~,ang] = mvdrspatialspect(signal);
ang = broadside2az(ang,25);


%% Beamscan DOA Estimation-URA
ura = phased.URA('Size',[8 8],'ElementSpacing',[1 1], 'Element',elem);
ang1 = [-45; -25];         
ang2 = [45; 25];         
signal = sensorsig(getElementPosition(ura)/lambda,Nsamp, ...
  [ang1 ang2],nPower);
rng(rs);                 
azelspectrum = phased.BeamscanEstimator2D('SensorArray',ura,...
       'OperatingFrequency',fc_cst,...
       'AzimuthScanAngles',-90:90,'ElevationScanAngles',-45:45,...
       'DOAOutputPort',true,'NumSignals',2);
[~,ang] = azelspectrum(signal);

figure(5)
spectruTest = plotSpectrum(azelspectrum);
title("URA-Beamscan")

%% MVDR DOA Estimation-URA
mvdrazelspectrum = phased.MVDREstimator2D('SensorArray',ura,...
        'OperatingFrequency',fc_cst,...
        'AzimuthScanAngles',-90:90,'ElevationScanAngles',-45:45,...
        'DOAOutputPort',true,'NumSignals',2);
[~,ang] = mvdrazelspectrum(signal);

figure(6)
plotSpectrum(mvdrazelspectrum);
title("URA-MVDR")

%% MUSIC DOA Estimation-URA
musicazelspectrum = phased.MUSICEstimator2D('SensorArray',ura,...
        'OperatingFrequency',fc_cst,...
        'AzimuthScanAngles',-90:90,'ElevationScanAngles',-45:45,...
        'DOAOutputPort',true,'NumSignalsSource','Property','NumSignals',2);
[~,ang] = musicazelspectrum(signal);

figure(7)
plotSpectrum(musicazelspectrum);
title("URA-MUSIC")
