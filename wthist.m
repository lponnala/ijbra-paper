function wthist(waiting_time)
%WTHIST    Create plot of datasets and fits
%   WTHIST(WAITING_TIME)
%   Creates a plot, similar to the plot in the main distribution fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1

% This function was automatically generated on 01-Sep-2009 22:09:17
 
% Data from dataset "waiting_time data":
%    Y = waiting_time
 
% Force all inputs to be column vectors
waiting_time = waiting_time(:);

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[441 243 680 475.45]);
legh_ = []; legt_ = {};   % handles and text for legend
ax_ = newplot;
set(ax_,'Box','on');
hold on;

% --- Plot data originally in dataset "waiting_time data"
t_ = ~isnan(waiting_time);
Data_ = waiting_time(t_);
[F_,X_] = ecdf(Data_,'Function','cdf'...
              );  % compute empirical cdf
Bin_.rule = 1;
[C_,E_] = dfswitchyard('dfhistbins',Data_,[],[],Bin_,F_,X_);
[N_,C_] = ecdfhist(F_,X_,'edges',E_); % empirical pdf from cdf
h_ = bar(C_,N_,'hist');
set(h_,'FaceColor','none','EdgeColor',[0 0 0],...
       'LineStyle','-', 'LineWidth',1);
xlabel('Data');
ylabel('Density')
legh_(end+1) = h_;
legt_{end+1} = 'waiting_time data';

% Nudge axis limits beyond data limits
xlim_ = get(ax_,'XLim');
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
end

x_ = linspace(xlim_(1),xlim_(2),100);

% --- Create fit "exponential distribution fit"

% Fit this distribution to get parameter values
t_ = ~isnan(waiting_time);
Data_ = waiting_time(t_);
% To use parameter estimates from the original fit:
%     p_ = [ 0.463540789907];
p_ = expfit(Data_, 0.05);
y_ = exppdf(x_,p_(1));
h_ = plot(x_,y_,'Color',[0 0 0],...
          'LineStyle','--', 'LineWidth',2,...
          'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_;
legt_{end+1} = 'exponential distribution fit';

hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'}; 
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
