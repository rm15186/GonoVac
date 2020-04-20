function [] = savepdf(scriptname,dirname)
%SAVEPDF Create a nice pdf figure.
%
%   SAVEPDF('script') will run the MATLAB script named script.m - which
%   should generate a single figure - and generate a reasonable-looking PDF
%   file from it, called script.pdf
%
%   SAVEPDF('script','dir') will save the output in a directory 'dir'
%   (default is the current directory)

if ~exist('dirname','var')
    % no output directory defined, output here
    dirname = './';
end

fontsize = 16;          % map all fonts to this size
outputsize = [18 12];   % output figure size [width height] in cm

fig = figure('PaperPosition',[0 0 outputsize],...
    'PaperSize',outputsize,...
    'PaperUnits','centimeters',...
    'Visible', 'Off');
disp(['Running ' scriptname]);
func = str2func(scriptname);
func();

axesHandle = gca;
set(axesHandle,'FontSize',fontsize)
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'FontSize',fontsize)

fname = [dirname scriptname '.pdf'];
disp(['Saving ' fname]);
print(fig,fname,'-dpdf');

disp(['Saved ' fname]);
close(fig);

end