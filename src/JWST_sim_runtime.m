%% JWST Dashboard Core Functionality
% A more open version of the JWST exercise, showing syntax needed to access
% outputs directly.
function [OPD_RMS, Strehl, RMS_Spotsize, MTF_data]  = simulate(CorrectionMatrix)
    % close all; clear all; clc;

    %% Load Lens File
    OS = JWST(1);
    sampling = 5;


    %% the Tolerance Matrix is a matrix with 4 columns and 18 rows (1 row per segment) 
    % Column 1: Piston tolerance [micron]
    % Column 2: Tilt around x [murad]
    % Column 3: Tilt around y [murad]
    % Column 4: RMS Focus Error [nm]
    % Feel free to replace with own randomized inputs, or ToleranceMatrix = zeros(18,4) 
    % for nominal performance.

    ToleranceMatrix(:,1) = [0.414170806677252  0.300243140031183 0.439976576044261 0.0152119184255706 -0.128208842154579 0.478345378762405 0.0754076647493227 0.651252356266327 0.101404608428473 0.0391465366416755 0.548089404530392 0.266178356301017 -0.360151894856671 0.385371897751932 0.300917231394459 -0.546668246684149 -0.238745289802985 -0.585858177909011];
    ToleranceMatrix(:,2) = [0.139276147276064  -0.390138294249314 -0.411334159967717 0.338255587537226 0.350357337374621 -0.244629702055557 -0.188140054852467 -0.208542872352273 0.252227970049943 0.0246373453963109 0.411647424007853 0.210408685278169 0.472651076977497 -0.182166946273771 0.0847186192633207 -0.390257631406097 -0.230116336295599 -0.0477925462370182];
    ToleranceMatrix(:,3) = [-0.214891914341358  -0.0239202817325437 0.326573979042765 0.0709910754624060 -0.438971070807491 -0.309013559302602 0.327732173448263 -0.292396965620019 0.176871093438419 -0.181895273849737 -0.352344222848263 -0.330232933973511 0.448108735396022 -0.366189014643874 0.171462889478031 -0.0574700377971163 -0.106588493632423 0.0846413033551107];
    ToleranceMatrix(:,4) = [-4.77158933834332,-5.90414896183012,6.46915662700510,12.3803995282438,-7.83571389646577,7.02863171336911,-7.60505047433170,12.3089743294982,7.55653924410722,-1.89433225241735,5.72159670128435,-0.0411604369112084,7.72475667997612,-3.58727666297385,-10.6689141377548,2.27478638187123,10.2546957682035,-7.65585159583474];

    %% Correction Matrix
    % The correction matrix. The array has the same specifications as the 
    % your optimization algortihm should ideally converge to a matrix that is
    % equivalent to the ToleranceMatrix, though there tend to be several
    % solutions (e.g. a small but consistent piston error across the pupil has a
    % neglicable performance impact)

    %CorrectionMatrix = zeros(size(ToleranceMatrix));


    %% Apply Tolerances and Corrections
    % The following code applies the tolerances and corrections to the optical
    % system. 

    OS_out = OS;      
    for seg = 1:18
        OS_out(1).surface(2).segment(seg) = Local2Global(OS_out(1).surface(2).segment(seg), [0 0 (CorrectionMatrix(seg, 1)-ToleranceMatrix(seg, 1))/1000], ...
            [CorrectionMatrix(seg, 2)-ToleranceMatrix(seg, 2), CorrectionMatrix(seg, 3)-ToleranceMatrix(seg, 3) 0]/1e6, OS_out(1).geometry(seg));
        if (CorrectionMatrix(seg, 4)- ToleranceMatrix(seg, 4))/1000 ~= 0
            OS_out(1).surface(2).segment(seg).zernike = ...
                {1000, [0 0 0 (CorrectionMatrix(seg, 4)- ToleranceMatrix(seg, 4))/1e6], OS_out(1).geometry(seg).position(1:2)};
        end

    end


    %% Calculate Outputs

    %% Optical Layout
    settings = FORTASettings;
    settings.field = 1;
    %Layout3D(OS_out,settings, 1, true, 2);

    %% Wavefront Diagram and PSF
    settings = FORTASettings;
    settings.field = 1;
    settings.plot = false; % if you want to run this in an optimizer it may be advisable to set this to false, otherwise you end up with lots of figures;
    settings.ref = 'Chief';
    settings.sampling = 5;
    settings.remove_tilt = 0;
    settings.initialpadding = 2*2^(sampling+4); % Zeropadding applied when calculating PSF.
    settings.viewfrac = 0.1; % default zoom level for PSF (visual effect only)

    %% Calculate Wavefront and Other Geometric Results:

    % Calculate first order properties (F/number, focal length, etc):
    OP = FirstOrderProperties(OS.Update.MakeRef,1);
    % Calculate Wavefront (and other parameters).
    % the final input is just the figure number.
    OPD_data = GetOPD(OS_out, OP, settings, 4); 

    % You can access the geometrical optical data as follows. These outputs may
    % be useful in an optimization function. 
    % OPD = OPD_data(1).field(1).OPD; % OPD data retrieved from struct;
    RMS_Spotsize = OPD_data(1).field(1).SPOT_RMS; % RMS Spotsize;
    % Full_Spotsize = OPD_data(1).field(1).SPOT_GEO; % Full Spotsize;
    % Centroid_X = OPD_data(1).field(1).CENX; % X location in the image plane of the centroid of the spot;
    % Centroid_Y = OPD_data(1).field(1).CENY; % Y location in the image plane of the centroid of the spot;
    OPD_RMS = OPD_data(1).field(1).RMS; % RMS Wavefront Error;
    % OPD_PV = OPD_data(1).field(1).PV; % Peak-to-Valley Wavefront Error;

    %% Calculate PSF and MTF

    % Calculate PSF and MTF from the wavefront data. Final input are the figure numbers (can be omitted)         
    [PSF_data, MTF_data] = FORREst(OS_out, settings, OP, OPD_data);

    PSF = PSF_data(1).field(1).PSFpoly; % extract point spread function.
%    PSF_diff = PSF_data(1).field(1).PSFpolydiff; % extract diffraction limited PSF.
%    ImageX = PSF_data(1).field(1).ImageX; % The x-axis of the image, in mm
%    ImageY = PSF_data(1).field(1).ImageY; % The y-axis of the image, in mm
    Strehl = PSF_data(1).field(1).Strehl; % The Strehl Ratio;
    
end