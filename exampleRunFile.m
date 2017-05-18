%% Example run file for the IDIC
% VARIABLES OPTIONS
% -------------------------------------------------------------------------
%   imagesFolder: subdirectory containing the series of images on which to run
%                 FIDIC. Images are read in alphanumeric order, which is
%                 assumed to be consistent with time steps.
%   Ext: the file extension of the input images.  All images must be of the
%        same type.
%   numImage: optional parameter to pass to img2mat limiting the number of
%             images processed from the subdirectory identified in imagesFolder
%   sSize: interrogation window (subset) size for the first iterations.
%          Must be 32,64,96, or 128 pixels and a two column
%          array (one for each dimenision) or scalar (equal for all
%          dimensions).
%   incORcum: string that defines the method of running IDIC. Options:
%             cumulative (time0 -> time1, time0 -> time2, ...)
%             (Allowable inputs: 'c','cum','cumulative')
%             or
%             incremental (time0 -> time1, time1 -> time2, ...)
%             (Allowable inputs: 'i','inc','incremental')
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u:  displacement field vector calculated from FIDIC. Format: cell array,
%      which is a 3D vector (components in x,y)  per each time point
%      (units are in pixels)
%         u{time}{1} = displacement in x-direction at t=time of size MxNxP
%         u{time}{2} = displacement in y-direction at t=time of size MxNxP
%         u{time}{3} = displacement magnitude at t=time of size MxNxP
%   cc: peak values of the cross-correlation for each interrogation
%   dm: final subset spacing in px
%
% NOTES
% -------------------------------------------------------------------------
%% Set up workspace and images

clear; close all; clc;
%dbstop if error

sSize = [64 64];
incORcum = 'c'; %use 'i' for incremental mode and 'c' for cumulative
norm_xcc = 'u'; %use 'norm' for normalized cross-correlation, considerable time-cost
ext_in = 'tif'; %Input image format
folder_in = ['.',filesep,'test_images'];
max_def_idx = 'b'; %Specify where the max deformation occurs
yn = 'y';
%use 'center' or 'c' for the center image,
%'end' or 'e' for the last image,
%'beginning' or 'b' for the first,
%or specific with an integer

%Compute basic noise floor and measurement resultion metrics
[noise_percent,meas_res,CI_disp_mean,no_im] = image_eval(folder_in,ext_in);

%Write outputs to screen and confirm run
if no_im == 0
    clc
    fprintf('IMAGE QUALITY DIAGNOSTICS\n-----------------------------------------\n')
    fprintf('Noise level: %0.2g%% \n',noise_percent)
    fprintf('Measurement resolution, x: %0.2gpx\n',meas_res(1))
    fprintf('Measurement resolution, y: %0.2gpx\n',meas_res(2))
    fprintf('Displacement confidence interval (95%%): %0.2g < %0.2gpx < %0.2g\n',...
        CI_disp_mean(1),(CI_disp_mean(1)-CI_disp_mean(2))/2,CI_disp_mean(2))
    fprintf('-----------------------------------------\n')
    yn = input('Are these acceptable? (y/n)[y] \n','s');
    
    if isempty(yn);yn='y';end
    if strcmp(yn(1),'n') || strcmp(yn(1),'N')
        return;
    end
    
else
    fprintf('IMAGE QUALITY DIAGNOSTICS UNAVAILABLE\n')
    fprintf('NO STATIC IMAGES: NOISE AND RESOLUTION UNKNOWN\n')
end

%Image cropping to get a region of interest on the specimen
[crop_nw_loc,folder_out] = imageCropping(folder_in,ext_in,sSize,max_def_idx,'on');

ext_crp = 'tif'; %output image file form, defined in image_cropping.m
resultsFolder = ['.',filesep,'Results',filesep];
numImages = 3;

%Convert input images to .mat and smooth
[cellIMG,filename,filt_opt] = img2mat(folder_out,ext_crp,'on'); %All images in "imagesFolder"
% [cellIMG,filename] = img2mat(folder_out,ext_crp,numImages); %Images 1 to
%numImages only

%% RUNNING DIC

% Estimate displacements via IDIC
[u, cc, dm] = funIDIC(filename, sSize, incORcum, norm_xcc);

% Save the results
if exist(resultsFolder,'dir') ~= 7
    mkdir(resultsFolder)
end

if no_im == 0
    %Build the reporting table struct array
    prefilt_str = strcat(filt_opt{1},', ',num2str(filt_opt{2}),', ',num2str(filt_opt{3}));
    reporting_table = struct('cameraNoise',noise_percent,'prefiltering',prefilt_str,...
        'subset',sSize,'step',dm,'xcorrType',norm_xcc,'interpolent','spline',...
        'numMeasurementPts',numel(u{1}{1}),'totalImages',length(u)+1,...
        'displacementSpatialRes',mean(sSize),'displacementResX',meas_res(1),...
        'displacementResY',meas_res(2));    
    %Save relavent workspace variables
    save(strcat(resultsFolder,'resultsFIDIC.mat'),'u','cc','cellIMG','dm','reporting_table');
else
    %Save relavent workspace variables
    save(strcat(resultsFolder,'resultsFIDIC.mat'),'u','cc','dm');
end

%% PLOTTING
close all;

def_udef = 'udef'; % set to 'udef' to plot on undeformed images,
% or 'def' to plot on deformed images

% Run the plotting routine function
FIDIC_plot(u,dm,def_udef,crop_nw_loc,folder_in,ext_in)

%save the last figure
saveas(gcf,strcat(resultsFolder,'FIDIC_plots.fig'))

if no_im == 0
%% Reporting Table
fprintf('\n-----------------------------------------\n');
fprintf('Run Parameters and Measurement Specifications\n')
fprintf('-----------------------------------------\n');
fprintf('Camera Noise \t\t %0.2g%%\n',noise_percent);
fprintf('Prefiltering \t\t %s, %0.2gx%0.2g, %0.2g\n',filt_opt{1},...
    filt_opt{2}(1),filt_opt{2}(2),filt_opt{3});
fprintf('Subset \t\t         %0.2g by %0.2gpx\n',sSize(1),sSize(2));
fprintf('Step         \t\t %0.2gpx\n',dm);
fprintf('Correlation type         %s\n',norm_xcc);
fprintf('Interpolation \t\t Spline\n');
fprintf('Measurement points \t %0.2g\n',numel(u{1}{1}));
fprintf('Total images \t\t %0.2g\n',length(u)+1);
fprintf('Displacement\n   Spatial resolution \t %0.2gpx \n   ',mean(sSize));
fprintf('Measurement res, x    %0.2g\n   ',meas_res(1));
fprintf('Measurement res, y    %0.2g\n',meas_res(2));
fprintf('-----------------------------------------\n');
else
end
%% CLEAN UP
%Clean up the current set of images from the cd
delete *IDIC_image*.mat
delete(strcat(folder_out,'*.',ext_crp));
