%% for reaction time by histogram
%function [speed_LL, accuracy_LL] = QMLE(data)
function LL = QMLE(simulate_data, empirical_data)
% QMLE, quantile maximum likelihood estimation
% reference: Heathcote & Australia, and Mewhort, 2002.

% speed
speed_h = figure;
%simulate_speed_rt = simulate_data.speed_rtmat;
simulate_speed_rt = simulate_data.rtmat;
%simulate_speed_choice = simulate_data.speed_choice;
simulate_speed_choice = simulate_data.choicemat;

speed_coh = empirical_data.speed_coh;
empirical_speed_qmat = empirical_data.speed_qmat;
speed_empirical_num = empirical_data.speed_observe_num;

speed_simulate_total = [];
speed_estimate_qnum = [];
speed_f = [];

for vi = 1: length(speed_coh)
    correct_RT = simulate_speed_rt(simulate_speed_choice(:, vi)==1, vi);
    %wrong_RT = simulate_speed_rt(simulate_speed_choice(:, vi)==0, vi);
    wrong_RT = simulate_speed_rt(simulate_speed_choice(:, vi)==2, vi);
    speed_simulate_total(vi) = length(simulate_speed_rt(:, vi));

    if ~any(isnan(empirical_speed_qmat(:, 1, vi)))
        temp = histogram(correct_RT, [0; empirical_speed_qmat(:, 1, vi); Inf]);
        speed_estimate_qnum(:, 1, vi) = temp.Values;
    else
        speed_estimate_qnum(:, 1, vi) = NaN(length(empirical_speed_qmat(:, 1, vi))+1, 1);
    end

    if ~any(isnan(empirical_speed_qmat(:, 2, vi)))
        temp = histogram(wrong_RT, [0; empirical_speed_qmat(:, 2, vi); Inf]);
        speed_estimate_qnum(:, 2, vi) = temp.Values;
    else
        speed_estimate_qnum(:, 2, vi) = NaN(length(empirical_speed_qmat(:, 2, vi))+1, 1);
    end

    speed_f(:, :, vi) = log(speed_estimate_qnum(:, :, vi) / speed_simulate_total(vi));
    speed_f(speed_f(:, 1, vi)==-Inf, 1, vi) = log(.1/speed_simulate_total(vi));
    speed_f(speed_f(:, 2, vi)==-Inf, 2, vi) = log(.1/speed_simulate_total(vi));
end

speed_LL = sum(speed_empirical_num(:) .* speed_f(:), 'omitnan');
LL.speed_LL = speed_LL;
close(speed_h);


% accuracy
accuracy_h = figure;
simulate_accuracy_rt = simulate_data.accuracy_rtmat;
simulate_accuracy_choice = simulate_data.accuracy_choice;

empirical_accuracy_qmat = empirical_data.accuracy_qmat;
accuracy_coh = empirical_data.accuracy_coh;
accuracy_empirical = empirical_data.accuracy_observe_num;

accuracy_total = [];
accuracy_estimate_qnum = [];
accuracy_f = [];

for vi = 1: length(accuracy_coh)
    correct_RT = simulate_speed_rt(simulate_accuracy_choice(:, vi)==1, vi);
    wrong_RT = simulate_speed_rt(simulate_accuracy_choice(:, vi)==0, vi);
    accuracy_total(vi) = length(simulate_accuracy_rt(:, vi));

    if ~any(isnan(empirical_accuracy_qmat(:, 1, vi)))
        temp = histogram(correct_RT, [0; empirical_accuracy_qmat(:, 1, vi); Inf]);
        accuracy_estimate_qnum(:, 1, vi) = temp.Values;
    else
        accuracy_estimate_qnum(:, 1, vi) = NaN(length(empirical_accuracy_qmat(:, 1, vi))+1, 1);
    end

    if ~any(isnan(empirical_accuracy_qmat(:, 2, vi)))
        temp = histogram(wrong_RT, [0; empirical_accuracy_qmat(:, 2, vi); Inf]);
        accuracy_estimate_qnum(:, 2, vi) = temp.Values;
    else
        accuracy_estimate_qnum(:, 2, vi) = NaN(length(empirical_accuracy_qmat(:, 2, vi))+1, 1);
    end

    accuracy_f(:, :, vi) = log(accuracy_estimate_qnum(:, :, vi) / accuracy_total(vi));
    accuracy_f(accuracy_f(:, 1, vi)==-Inf, 1, vi) = log(.1/accuracy_total(vi));
    accuracy_f(accuracy_f(:, 2, vi)==-Inf, 2, vi) = log(.1/accuracy_total(vi));
end

accuracy_LL = sum(accuracy_empirical(:) .* accuracy_f(:), 'omitnan');
LL.accuracy_LL = accuracy_LL;
close(accuracy_h);