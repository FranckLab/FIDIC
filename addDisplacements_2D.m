function [u, du, cc] = addDisplacements_2D(u0,du0,cc0,m0,dm)
% u = addDisplacements(u0,du0,cc0,m0,dm) adds displacements from previous
% iteratations to the current iterate's displacement field (du).
%
% INPUTS
% -------------------------------------------------------------------------
%   u0: past displacement field vector defined at every meshgrid point with
%      spacing dm. Format: cell array, each containing a 2D matrix
%         (components in x,y)
%         u0{1} = displacement in x-direction
%         u0{2} = displacement in y-direction
%         u0{3} = magnitude
%   du0: current displacements as
%         du0{1} = displacement in x-direction
%         du0{2} = displacement in y-direction
%   cc0: cross correlation matrix used to define interpolant locations
%   m0: mesh grid parameter
%   dm: subset spacing paramter used to set up grids
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the added displacement fields
%   du: cell containing interpolated incremental displacements
%   cc: cross-correlation matrix interpolated to fit with u and du
%
% NOTES
% -------------------------------------------------------------------------
%
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast
% iterative digital volume correlation algorithm for large deformations.
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

for i = 1:2, du0{i} = inpaint_nans(du0{i}); end % remove NaNs if present

idx = cell(1,2);
for i = 1:2, idx{i} = m0{i}(1):dm:m0{i}(end); end % construct new meshgrid

[m0_{1}, m0_{2}] = ndgrid(m0{1},m0{2});
[m{1}, m{2}] = ndgrid(idx{1},idx{2});
% sample to desired mesh spacing

try
%Try to interpolate the displacement field with griddedInterpolant.
%This function does not exist on older versions of MATLAB (2011b and
%earlier), so fail over to interpn, with is slower and prone to other
%failures.
    du = cell(1,2);
    for i = 1:2
        F = griddedInterpolant(m0_{1}, m0_{2}, du0{i}, 'linear');
        du{i} = F(m{1},m{2});
    end
    
    F = griddedInterpolant(m0_{1}, m0_{2}, cc0, 'linear');
    cc = F(m{1},m{2});
    
catch
%Attempt to use interpn instead
    disp('Attempting interpn')
    du = cell(1,2);
    for i = 1:2
        du{i} = interpn(m0_{1}, m0_{2}, du0{i}, m{1}, m{2}, 'linear');
        %du{i} = V(m{1},m{2});
    end
    
    cc = interpn(m0_{1}, m0_{2},cc0,m{1},m{2}, 'linear');
    
end

if  sum(cellfun(@numel,u0)) == 2, u = du; % on first iteration u = du
else
    u = cellfun(@plus,u0,du,'UniformOutput',0); % else u^(k) = sum(u^(k-1)) + du (see eq. 7)
    
end

end
