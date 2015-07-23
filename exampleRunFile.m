%% Example run file for the IDIC - LDTFM
% VARIABLES OPTIONS
% -------------------------------------------------------------------------
%   imagesFolder: subdirectory contain the series of images on which to run
%                 FIDIC. Images are read in alphanumeric order, which is 
%                 assumed to be consistent with time steps.
%   Ext: the file extension of the input images.  All images must be of the
%        same type.       
%   numImage: optional parameter to pass to tiff2mat limiting the number of
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

clear; close all; clc;

sSize = [64 64];
incORcum = 'c';
imagesFolder = 'test_images';
resultsFolder = '.\Results\';
Ext = 'tif';
numImages = 4;

[cellIMG,filename] = tiff2mat(imagesFolder,Ext); 
% [cellIMG,filename] = tiff2mat(imagesFolder,Ext,numImages); 

% Estimate displacements via IDIC
[u, cc, dm] = funIDIC(filename, sSize, incORcum);

if exist(resultsFolder,'dir') ~= 7
    mkdir(resultsFolder)
end
save(strcat(resultsFolder,'resultsFIDIC.mat'),'u','cc','cellIMG','dm');

%Clean up the current set of images from the cd
delete *IDIC_image*.mat


%% PLOTTING 
numInc = max(size(u));
scrsz = get(0,'ScreenSize');
sizeI = size(cellIMG{1});
plotIdx = cell(1,2);
for ii = 1:2, plotIdx{ii} = 1:dm:sizeI(ii)+1; end

close all;
% plot IDIC results

for jj = 1:numInc
    figure
    for ii = 1:3 
        subplot(2,2,ii);
        
        u_ = u{jj}{ii};
        u_field = u_(dm/2:end-dm/2,dm/2:end-dm/2);
        u_field = padarray(u_field,dm/2,nan,'both'); %Crop (some) edge effects
        [~,h] = contourf(u_field,25); colorbar; colormap jet;
        set(h,'linestyle','none'); axis image
      
        if ii == 3
          title('Displacement magnitude')  
        else
            title(['Displacement component u_',num2str(ii)])
        end
        xlabel('X_1'); ylabel('X_2');
    end
end 

