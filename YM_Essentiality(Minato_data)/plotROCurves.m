

%% FIGURE SETTINGS

% Defaults for this blog post
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1;      % LineWidth
msz = 8;       % MarkerSize


% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);


pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

figure(1)
plot(FPR,SENSITIVITY,'-r*')
               % 'MarkerEdgeColor','r',...
               %  'MarkerFaceColor',[.49 1 .63],...
               % 'MarkerSize',6)
hold on
plot(x,y,'--k')
        
legend('zbar > 0.9925 (ES)... 0 <= zbar < 0.0493','Baseline', 'Location', 'SouthEast');
xlabel('1 - Specificity');
ylabel('Sensitivity');
title(METABOLIC_MODEL);

print(METABOLIC_MODEL, '-dpng', '-r300');



