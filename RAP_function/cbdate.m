function [h] = cbdate(varargin) 
% cbdate formats colorbar ticks as dates strings. 
% 
%% Syntax
% 
%  cbdate
%  cbdate(h)
%  cbdate(datenums)
%  cbdate(datenums,dateformat) 
%  cbdate(datestrs)
%  cbdate('horiz') 
%  h = cbdate(...) 
% 
%% Description
% 
% cbdate converts tick marks of current colorbar to date strings. 
%
% cbdate(h) specifies a colorbar handle h on which to perform cbdate. If
% no handle is specified, cbdate will try to find a current colorbar. 
% 
% cbdate(datenums) specifies datenums to set as tick marks. 
%
% cbdate(datenums,dateformat) specifies date formatting as described in the 
% documentation for datestring. For example, 'mmm yyyy' will print date
% strings as Nov 2014. 
%
% cbdate(datestrs) specifies date strings to use as colorbar ticks. 
% 
% cbdate('horiz') must be declared if the colorbar is horizontal
% (pre-R2014b). 
% 
% h = cbdate(...) returns a colorbar handle h. 
% 
%% Examples 
% % You've plotted some data with its color scaled as a function of date. The
% % plot looks like this: 
% 
% t = datenum(2014,1:24,15); 
% y = sin(t/200); 
% 
% figure('pos',[100 100 700 400])
% scatter(t,y,60,t,'filled'); 
% datetick('x','yyyy')
% 
% colorbar
% 
%% 7.361 x 10^5 does not mean much to most folks. Let's change those datenum
% % ticks on the colorbar to date strings: 
% 
% cbdate
%% Your original data were monthly data, centered on the 15th of each month,
% % so it does not make much sense to have colorbar tick labels label days.
% % Monthly resolution is about the best we can describe. So let's format
% % those tick labels to show only the month: 
% 
% cbdate('mmmm')
% 
%% Or perhaps you'd like your tick labels referenced to specific input data points. 
% % The variable |t| is in |datenum| format, so it can be used as an input directly to |cbdate|. 
% % Here we label every third month:
% 
% cbdate(t(1:3:end))
% 
%% Or perhaps you want to label two specific dates:
% 
% cbdate({'apr 1, 2014';'jan 12, 2015'})
% 
%% Perhaps your colorbar is horizontal: 
% 
% colorbar('location','north');
% cbdate('horiz','mmm-yy')
% 
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas
% at Austin's Institute for Geophysics (UTIG) on datestr(735913.51). 
% Come say hello at http://www.chadagreene.com. 
% 
% See also: datenum, datestr, colorbar. 

%% Input checks: 

assert(nargin<5,'cbdate cannot have more than four inputs.') 

%% Set defaults: 

DateFormatDeclared = false; 

% Default colorbar orientation is vertical: 
VerticalColorbar = true; 
if any(strncmpi(varargin,'hor',3))
    VerticalColorbar = false; 
end

% Default colorbar handle is the current colorbar handle: 
h = findobj('tag','Colorbar'); % <-- Only for pre 2014b???
if isempty(h)
    h = colorbar; 
else
    h = h(1); 
end


%% Parse inputs: 

for k = 1:length(varargin) 
    
    % If varargin{k} is a handle, assume it's a colorbar handle: 
    if ishandle(varargin{k})
        h = varargin{k}; 
    end
    
    % If varargin{k} is numeric and it's not a handle, assume it's datenums: 
    if isnumeric(varargin{k}) & ~ishandle(varargin{k})
        datenums = varargin{k}; 
    end
    
    % If varargin{k} is not numeric, nor the 'horizontal' tag, nor a date
    % formatting preference, assume it's date strings.  
    if ~isnumeric(varargin{k}) & size(varargin{k},1)>1 & ~strncmpi(varargin{k},'hor',3)
        datestrs = varargin{k}; 
        datenums = datenum(datestrs); 
    end
    
    % If varargin{k} is not numeric, nor the 'horizontal' tag, nor date
    % strings, assume it's date formatting: 
    if ~isnumeric(varargin{k}) & size(varargin{k},1)==1 & ~strncmpi(varargin{k},'hor',3)
        DateFormat = varargin{k}; 
        DateFormatDeclared = true; 
    end
end

% If no datenums have been declared by user, use current colorbar ticks:  
if exist('datenums','var')==0
    if verLessThan('matlab','8.4.0')
        if VerticalColorbar
            datenums = get(h,'YTick'); 
        else
            datenums = get(h,'XTick');
        end
    else
        datenums = h.Ticks; 
    end
end
    
% If no date strings have been declared by user, convert datenums to date strings:  
if ~exist('datestrs','var')
    if DateFormatDeclared
        datestrs = datestr(datenums,DateFormat);
    else
        datestrs = datestr(datenums); 
    end
end

%% Set ticks and labels: 


if verLessThan('matlab','8.4.0')
    if VerticalColorbar
        set(h,'Ytick',datenums,'YTickLabel',datestrs);
    else
        set(h,'Xtick',datenums,'XTickLabel',datestrs);
    end
else
    h.Ticks = datenums;
    h.TickLabels = datestrs;
end

%% Clean up: 

if nargout==0 
    clear h
end


end
