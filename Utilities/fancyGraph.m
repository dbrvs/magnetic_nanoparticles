function fancyGraph(fig, varargin)
% Improves the visibility of a given figure by changing line-width and
% text-size of all nested objects.
% Additionaly, it may save the figure in epsc format, with the name 'filename',
% as in command 'saveas'
%
% Usage:
%       fancyGraph(fig, filename, fileformat, textsize, linewidth)
%       fancyGraph(fig, filename, textsize, linewidth)
%       fancyGraph(fig)
%       fancyGraph(fig, textsize, linewidth)
%
%
%       Without specifying textsize and linewidth, these default to 15 and
%       2, respectively.
%       File format defaults to epsc.
%
% Examples:
%       fancyGraph(gcf);  % Use current figure
%
% 2011 - Manuel Duarte

if nargin == 0
    error('Specify figure handle')
elseif nargin == 1
    textsize = 12;
    linewidth = 1;
elseif nargin == 3
    textsize = varargin{1};
    linewidth = varargin{2};
    filename = '';
    fileformat = '';
elseif nargin == 4
    filename = varargin{1};
    textsize = varargin{2};
    linewidth = varargin{3};
    fileformat = 'epsc';
elseif nargin == 5
    filename = varargin{1};
    fileformat = varargin{2};
    textsize = varargin{3};
    linewidth = varargin{4};
elseif (nargin == 2) || (nargin > 5)
    error('Invalid arguments. Refer to help fancyGraph')
end

children = findall(fig);

for k = 1:numel(children)
    try
        set(children(k), 'LineWidth',linewidth);
    end
    try
        set(children(k), 'FontSize',textsize)
    end
end

if nargin >= 4
    saveas(fig, filename, fileformat);
end
