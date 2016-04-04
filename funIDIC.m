function [u, cc, dm] = funIDIC(varargin)
% u = funIDIC(filename, sSize, incORcum) is the main function that performs
% IDIC on a time increment of images.  
% 
% INPUTS
% -------------------------------------------------------------------------
%   filename: string for the filename prefix for the images in 
%             the current directory.  
%             Input options:
%             --- If image is not within a cell) ---
%             1) 'filename*.mat' or 'filename*' 
% 
%             --- If image is within a cell that contains multichannels ---
%             2) filename{1} = 'filename*.mat' or 'filename*' and
%                filename{2} = channel number containing images you want to
%                              run IDIC on.
%                (if the channel is not provided, i.e. length(filename) = 1
%                , then channel = 1
%                 
%   sSize: interrogation window (subset) size for the first iterations.  
%          Must be, 32,64,96, or 128 voxels and a three column
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
%   u: displacement field vector defined at every meshgrid point with 
%      spacing dm. Format: cell array, each containing a 3D matrix for each
%      time point
%         (components in x,y)
%         u{time}{1} = displacement in x-direction
%         u{time}{2} = displacement in y-direction
%         u{time}{3} = magnitude
%   cc: peak values of the cross-correlation for each interrogation point
%   dm: meshgrid spacing (8 by default)
% 
% NOTES
% -------------------------------------------------------------------------
% none
% 
% For more information please see 
%

%% ---- Opening & Reading the First Image into CPU Memory ----
[fileInfo, sSize0, incORcum, u_, u_c] = parseInputs(varargin{:});

I{1} = loadFile(fileInfo,1);

%% ---- Opening and Reading Subsequent Images ---
numImages = length(fileInfo.filename);
u = cell(numImages-1,1); uc = cell(numImages-1,1); cc = cell(numImages-1,1);
for i = 2:numImages % Reads images starting on the second image
    tStart = tic;
    I{2} = loadFile(fileInfo,i);

    %Start DIC
    disp(['Current file: ' fileInfo.filename{i}])
    [u_, u_c, cc{i-1}, dm] = IDIC(I,sSize0,u_c);
    
    % Saving iterations of the DIC
    u{i-1}{1} = u_{1};  u{i-1}{2} = u_{2};  u{i-1}{3} = u_{3};
    uc{i-1}{1} = u_c{1};  uc{i-1}{2} = u_c{2};

    u_ = num2cell(zeros(1,2));
    if strcmpi(incORcum(1),'i') == 1
        I{1} = I{2};
    else
        u_c = uc{i-1}(1:2);
    end
    
    disp(['Elapsed Time for all iterations: ',num2str(toc(tStart))]);
    
end

end

%===================================================================
function I = loadFile(fileInfo,idx)

I = load(fileInfo.filename{idx});
fieldName = fieldnames(I);
I = getfield(I,fieldName{1});
if iscell(I),
    if numel(I), I = I{1};
    else
        I = I{fileInfo.dataChannel};
    end
end
end

%================================================================
function varargout = parseInputs(varargin)
%  = parseInputs(filename, sSize, incORcum)


% Parse filenames
filename = varargin{1};
if iscell(filename)
    if length(filename) == 1, fileInfo.datachannel = 1;
    else fileInfo.datachannel = filename{2};
    end
    filename = filename{1};
end


[~,filename,~] = fileparts(filename);
filename = dir([filename,'.mat']);
fileInfo.filename = {filename.name};

if isempty(fileInfo), error('File name doesn''t exist'); end

% Ensure dimensionality of the subset size
sSize = varargin{2};
if numel(sSize) == 1,
    sSize = sSize*[1 1];
elseif numel(sSize) ~=2,
    error('Subset size must be a scalar or a two column array');
end

% Ensure range of subset size
if min(sSize) < 32 || max(sSize > 128)
   error('Subset size must be within 32 and 128 pixels'); 
end

% Ensure even subset size
% if sum(mod(sSize,4)) > 0
%     error('Subset size must be even');
% end

if sum(mod(sSize,32)) ~= 0
    error('Subset size must be 32, 64, 96, or 128 pixels in each dimension');
end

% Check run method input
incORcum  = varargin{3};

switch lower(incORcum)
    case 'cum', incORcum = 'cumulative';
    case 'inc', incORcum = 'incremental';
    case 'c', incORcum = 'cumulative';
    case 'i', incORcum = 'incremental';
    case 'incremental', incORcum = 'incremental';
    case 'cumulative', incORcum = 'incremental';
    otherwise, error('Run method must be incremental or cumulative');
end

% Initial guess of displacement field = [0 0];
u0 = num2cell(zeros(1,2));
uc0 = num2cell(zeros(1,2));
% Outputs
varargout{      1} = fileInfo;
varargout{end + 1} = sSize;
varargout{end + 1} = incORcum;
varargout{end+1} = u0;
varargout{end+1} = uc0;

end