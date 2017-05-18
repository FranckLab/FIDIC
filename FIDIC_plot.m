function [] =  FIDIC_plot(u,dm,def_udef,crop_nw_loc,imagesFolder,ext_in)
%This function constructs overlay plots showing the displacement contours
%on the undeformed or deformed images
%
% INPUTS
% -------------------------------------------------------------------------
%   u: cell array of the displacement fields
%   dm: final subset spacing from FIDIC
%   def_udef: flag to specify plotting on deformed or undeformed images
%   crop_nw_loc: northwest corner location of the crop
%   imagesFolder: Directory where the complete images reside

%% Set up local vars
numInc = max(size(u));
scrsz = get(0,'ScreenSize');
sizeI = size(u{1}{1});

if strcmp(def_udef,'def')
    im_idx_offset = 1;
else
    im_idx_offset = 0;
end

%% Read in image filenames
files = dir(strcat(imagesFolder,filesep,'*',ext_in));
l = length(files);

%% plot IDIC results

for jj = 1:numInc
    
    %read in the associated image
    cur_img = imread(strcat(imagesFolder,filesep,files(jj+im_idx_offset).name));
    try
        cur_img = rgb2gray(cur_img);
    catch
    end
    img_size = size(cur_img);
    disp_img = nan*zeros(img_size);
    %     disp_mag_ave = mean(u{jj}{3}(:));
    %     disp_mag_std = std(u{jj}{3}(:));
    %     disp_mag_range = max(max(medfilt2(u{jj}{3}))) - min(min(medfilt2(u{jj}{3})));
    
    figure
    set(gcf,'position',[150,150,scrsz(3)*(7/8),scrsz(4)*3/4])
    for ii = 3:-1:1
        try
            %Set up axes for the image and for the contour plot
            a1 = axes;
            set(a1,'Units','Pixels')
            
            a2 = axes;
            set(a2,'Units','Pixels')
            
%             if ii < 3
            subplot(1,3,ii,a1);
            subplot(1,3,ii,a2);
%             end
            
            format short
            
            %get the displacement data for the current case, filter and mask it
            u_ = medfilt2(u{jj}{ii});
            u_([1:dm/2,(end-dm/2+1:end)],:) = nan;
            u_(:,[1:dm/2,(end-dm/2+1:end)]) = nan; %Crop (some) edge effects
            
            
            u_upscale = upsampleImage(u_,dm);
            %set up the image of the displacements, padded to the size of the
            %orginal images
            disp_img(crop_nw_loc(2):(crop_nw_loc(2)+size(u_upscale,1)-1),...
                crop_nw_loc(1):(crop_nw_loc(1)+size(u_upscale,2)-1)) = u_upscale;
            
            %Set up the alpha-belnding mask so that the contour underneath the
            %image can be seen
            alpha_data = ones(img_size);
            %alpha_data((crop_nw_loc(2)+dm*dm/2):(crop_nw_loc(2)+size(u_upscale,1)-dm*dm/2-1),...
            %   (crop_nw_loc(1)+dm*dm/2):(crop_nw_loc(1)+size(u_upscale,2)-dm*dm/2)-1) = 0.5;
            
            alpha_data(alpha_data>0.5) = 0.65;
            
            %Find scaling parameters needed for the colorbar
            disp_mag_ave = nanmean(u_(:));
            disp_mag_std = nanstd(u_(:));
            if ii == 3
                disp_mag_range = max(u_(:)) - min(u_(:));
                tick_labels = (disp_mag_ave-disp_mag_range):...
                    disp_mag_range/2:(disp_mag_ave+disp_mag_range);
            end
            
            %Do the image plotting
            set(gcf,'currentaxes',a2)
            p = imshow(cur_img);
            set(p,'AlphaData',alpha_data); %use the predefined transparency mask
            axis image
            
            %Do the contour plotting
            set(gcf,'currentaxes',a1)
            [~,h] = contourf(disp_img); % flip the contour about x and y
            % to be in the same reference
            %  frame as the image
            set(h,'linestyle','none');
            colormap('parula')
            axis image
            
            %set appropriate colormasks
            colormap(a1,'parula')
            colormap(a2,'gray')
            
            %Save the current axis position
            set(gcf,'CurrentAxes',a1)
            pos_a1 = get(a1,'position');
            pos_a2 = get(a2,'position');
            pos_a1(3) = pos_a1(3) - 50;
            pos_a2(3) = pos_a2(3) - 50;
            
            %set up the plots with contours and colorbar tick scaled by the
            %magnitude of the range in the Disp. Mag. plot.
            lvl_step = disp_mag_range/6; %Set number of contour levels by step size
            if ii == 3
                title('Displacement magnitude')
                
                %Define and set the contour plot levels and colorbar levels
                lvls_3 = (disp_mag_ave-disp_mag_range):...
                    lvl_step:(disp_mag_ave+disp_mag_range);
                h.LevelList = lvls_3;
                colorbar('Ticks',lvls_3,'Limits',[lvls_3(1),lvls_3(end)])
                set(a1,'position',pos_a1); % restore position
                set(a2,'position',pos_a2); % restore position
                
            elseif ii == 2
                title(['Displacement component u_',num2str(ii)])
                
                %Set the levels, using the step sizes from the magitude, but
                %center from the current plot
                lvls_2 = (disp_mag_ave-disp_mag_range):...
                    lvl_step:(disp_mag_ave+disp_mag_range);
                h.LevelList = lvls_2;
                colorbar('Ticks',lvls_2)
                set(a1,'position',pos_a1); % restore position
                set(a2,'position',pos_a2); % restore position
                
            elseif ii == 1
                title(['Displacement component u_',num2str(ii)])
                lvls_1 = (disp_mag_ave-disp_mag_range):...
                    lvl_step:(disp_mag_ave+disp_mag_range);
                h.LevelList = lvls_1;
                colorbar('Ticks',lvls_1);
                set(a1,'position',pos_a1); % restore position
                set(a2,'position',pos_a2); % restore position
            end
            
            %label the axes
            xlabel('X_1 [px]'); ylabel('X_2 [px]');
            
            %         set(gcf,'CurrentAxes',a1)
            %         p=getframe; %get a frame from the current axes.
            %         image(p.cdata); % display the data as an image
            %         axis image;
            %         impixelinfo; % turn on pixel information
            %         set(a1,'position',pos_a1); % restore position
            %         set(a2,'position',pos_a2); % restore position
        catch
        end
    end
end

function [I] = upsampleImage(comp_image, samp_scale)
%This function takes in an image and upsamples it by building blocks of
%pixels.  It takes in a compressed image and scaling factor to generate a
%final image as output


if nargin < 2
    samp_scale = 8;
end

l = size(comp_image,1);
h = size(comp_image,2);
L = l*samp_scale;
H = h*samp_scale;
I = zeros(L,H);

for ii = 1:l
    for jj = 1:h
        
        I(((ii-1)*samp_scale+1):(ii*samp_scale),...
            ((jj-1)*samp_scale+1):(jj*samp_scale)) = comp_image(ii,jj);
        
    end
end
