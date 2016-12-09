function disperror1()
%disperror1
%   Display error message if can't plot phasors

    errmsg = uicontrol('Style','text','Position',[130 400 400 20],...
    'string','ERROR! Try decreasing Load Magnitude or Line Length.',...
    'backgroundcolor',[1 1 1],'fontsize',12);

end

