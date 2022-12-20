% Q-Q plot for reaction time and choice
function QQ_plot(dataBhvr , monkey_name)
    lwd = 1.0;
    mksz = 3;
    fontsize = 11;
    
    % speed
    speed_x = dataBhvr.speed_proportionmat;
    speed_y = dataBhvr.speed_qmat;
    speed_h = figure; hold on;
    
    % filename = sprintf('Q-QPlot_%s',name);
    for vi = 1:length(speed_x)
        x_correct_coords = speed_x(vi) * ones(size(speed_y(:,1,vi)));
        x_wrong_coords = (1-speed_x(vi)) * ones(size(speed_y(:,2,vi)));
        plot(x_correct_coords, speed_y(:,1,vi), 'gx', 'MarkerSize', mksz+1, 'LineWidth', lwd);
        plot(x_wrong_coords, speed_y(:,2,vi), 'rx', 'MarkerSize', mksz+1, 'LineWidth', lwd);
    
      % fitted value
    %   En(vi) = numel(rtmat(:,vi));
    %   RT_corr = rtmat(choicemat(:,vi) == 1,vi);
    %   RT_wro = rtmat(choicemat(:,vi) == 2,vi);
    %   xr = numel(RT_corr)/(numel(RT_corr) + numel(RT_wro));
    %   q(:,1,vi) = quantile(RT_corr,qntls); % RT value on quantiles, correct trial
    %   q(:,2,vi) = quantile(RT_wro,qntls); % RT value on quantiles, error trial
    end
    
    % for qi = 1:size(q,1)
    %   xq = [flip(1-x), x]';
    %   plot(xq,[squeeze(flip(q(qi,2,:)));squeeze(q(qi,1,:))],'k-o','MarkerSize',mksz,'LineWidth',lwd/2);
    % end
    
    xlim([-.05 1.05]);
    % ylim([0.2, 1.4]);
    % yticks([.2:.4:1.4]);
    title('Speed')
    xlabel('Proportion');
    ylabel('RT (s)');
    savefigs(speed_h, monkey_name+"_Speed.eps", "C:\Users\weiyi\OneDrive\文档\LDDM Speed Accuracy Tradeoff\Data Fitting\Hanks", fontsize, [2.5 2.5]);
    
    % accuracy
    accuracy_x = dataBhvr.accuracy_proportionmat;
    accuracy_y = dataBhvr.accuracy_qmat;
    accuracy_h = figure; hold on;
    
    % filename = sprintf('Q-QPlot_%s',name);
    for vi = 1:length(accuracy_x)
        x_correct_coords = accuracy_x(vi) * ones(size(accuracy_y(:,1,vi)));
        x_wrong_coords = (1-accuracy_x(vi)) * ones(size(accuracy_y(:,2,vi)));
        plot(x_correct_coords, accuracy_y(:,1,vi), 'gx', 'MarkerSize', mksz+1, 'LineWidth', lwd);
        plot(x_wrong_coords, accuracy_y(:,2,vi), 'rx', 'MarkerSize', mksz+1, 'LineWidth', lwd);
    
      % fitted value
    %   En(vi) = numel(rtmat(:,vi));
    %   RT_corr = rtmat(choicemat(:,vi) == 1,vi);
    %   RT_wro = rtmat(choicemat(:,vi) == 2,vi);
    %   xr = numel(RT_corr)/(numel(RT_corr) + numel(RT_wro));
    %   q(:,1,vi) = quantile(RT_corr,qntls); % RT value on quantiles, correct trial
    %   q(:,2,vi) = quantile(RT_wro,qntls); % RT value on quantiles, error trial
    end
    
    % for qi = 1:size(q,1)
    %   xq = [flip(1-x), x]';
    %   plot(xq,[squeeze(flip(q(qi,2,:)));squeeze(q(qi,1,:))],'k-o','MarkerSize',mksz,'LineWidth',lwd/2);
    % end
    
    xlim([-.05 1.05]);
    % ylim([0.2, 1.4]);
    % yticks([.2:.4:1.4]);
    title('Accuracy')
    xlabel('Proportion');
    ylabel('RT (s)');
    savefigs(accuracy_h, monkey_name+"_Accuracy.eps", "C:\Users\weiyi\OneDrive\文档\LDDM Speed Accuracy Tradeoff\Data Fitting\Hanks", fontsize, [2.5 2.5]);
end