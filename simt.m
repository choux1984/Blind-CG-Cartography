function simt(onlyplot,fnum_in,niter)
% This is the master file.
%
%    Input:
%        ONLYPLOT   0: Simulate
%                   1: Display results (previously computed)
%        FNUM_IN    index of the figure within the function file pointed by FS below        
%

% Initializations
addpath('./lib/');
addpath('./cartographySimulator/');
initializeCartographySimulator;
%initCvx;
assert(nargin>=1,'Not enough arguments');

% Figure options (EDIT)
global plot_f_settings
plot_f_settings.docked_on=0; 
plot_f_settings.title_to_caption=1;
plot_f_settings.saveplots=0;  % write figures to files
plot_f_settings.pdf_flag=1;   % write figures to pdf
plot_f_settings.bw_maxpoints=20; 
plot_f_settings.figbasename='';  % if ~='' then figures are stored in the specified location,
                                 % else, the default folder is used. It must end with '/'.
plot_f_settings.pos=[200.0000  202.0000 .85*[592.3077  363.4615]]; % good for papers
%plot_f_settings.pos=[];
global chars_per_line
chars_per_line = 40;


% Execution options (EDIT)
%defaultFigureIndex = [3402];
defaultFigureIndex = 5201; %[1006, 3100, 3402, 3232];
if nargin<3
	niter = 500;
end
% name of the file containing the simulations
fs = CGCartographySimulations;  


% EXECUTION
if nargin < 2
	simt_compute(fs,defaultFigureIndex,niter,onlyplot)
else
	simt_compute(fs,fnum_in,niter,onlyplot)
end


end



function initCvx

cf = '../lib';
addpath([cf '/cvx'])
addpath([cf '/cvx/structures'])
addpath([cf '/cvx/lib'])
addpath([cf '/cvx/functions'])
addpath([cf '/cvx/commands'])
addpath([cf '/cvx/builtins'])


end

