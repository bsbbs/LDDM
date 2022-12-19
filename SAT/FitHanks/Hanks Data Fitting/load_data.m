function data = load_data(file_name)
% function data = LoadRoitmanData(datadir)
% addpath(datadir);
load(file_name);

% define parameters for method QMLE, 
% reference: Heathcote & Australia, and Mewhort, 2002.
qntls = .1:.1:.9; % quantile probability
expd_qntls = [0, qntls, 1];
P = expd_qntls(2:end) - expd_qntls(1:end-1);


% speed
speed = task(1).data;
speed_coh = unique(speed.coh);
speed_qmat = [];
speed_proportionmat = [];
speed_observe_num = [];
speed_observe_total = [];

for vi = 1:length(speed_coh)
    cur_RT = speed.rt(speed.coh == speed_coh(vi)) / 1000;
    speed_observe_total(vi) = length(cur_RT);

    if 0 %vi == 1
        speed_proportionmat(vi) = .5;
        speed_qmat(:, 1, vi) = quantile(cur_RT, qntls);
        speed_qmat(:, 2, vi) = quantile(cur_RT, qntls);
        speed_observe_num(:, 1, vi) = length(cur_RT, qntls) * P / 2;
        speed_observe_num(:, 2, vi) = length(cur_RT, qntls) * P / 2;

    else
        correct_RT = speed.rt(speed.coh == speed_coh(vi) & speed.cor == 1) / 1000;
        wrong_RT = speed.rt(speed.coh == speed_coh(vi) & speed.cor == 0) / 1000;
        speed_proportionmat(vi) = length(correct_RT) / length(cur_RT);
        speed_qmat(:, 1, vi) = quantile(correct_RT, qntls);
        speed_qmat(:, 2, vi) = quantile(wrong_RT, qntls); 
        speed_observe_num(:, 1, vi) = length(correct_RT) * P;
        speed_observe_num(:, 2, vi) = length(wrong_RT) * P;
    end
end

data.speed = speed;
data.speed_coh = speed_coh;
data.speed_qmat = speed_qmat;
data.speed_proportionmat = speed_proportionmat;
data.speed_observe_num = speed_observe_num;
data.speed_observe_total = speed_observe_total;


% accuracy
accuracy = task(2).data;
accuracy_coh = unique(accuracy.coh);
accuracy_qmat = [];
accuracy_proportionmat = [];
accuracy_observe_num = [];
accuracy_observe_total = [];

for vi = 1:length(accuracy_coh)
    cur_RT = accuracy.rt(accuracy.coh == accuracy_coh(vi)) / 1000;
    accuracy_observe_total(vi) = length(cur_RT);

    if vi == 1
        accuracy_proportionmat(vi) = .5;
        accuracy_qmat(:, 1, vi) = quantile(cur_RT, qntls);
        accuracy_qmat(:, 2, vi) = quantile(cur_RT, qntls);

    else
        correct_RT = accuracy.rt(accuracy.coh == accuracy_coh(vi) & accuracy.cor == 1) / 1000;
        wrong_RT = accuracy.rt(accuracy.coh == accuracy_coh(vi) & accuracy.cor == 0) / 1000;
        accuracy_proportionmat(vi) = length(correct_RT) / length(cur_RT);
        accuracy_qmat(:, 1, vi) = quantile(correct_RT, qntls);
        accuracy_qmat(:, 2, vi) = quantile(wrong_RT, qntls); 
        accuracy_observe_num(:, 1, vi) = length(correct_RT) * P;
        accuracy_observe_num(:, 2, vi) = length(wrong_RT) * P;
    end
end

data.accuracy = accuracy;
data.accuracy_coh = accuracy_coh;
data.accuracy_qmat = accuracy_qmat;
data.accuracy_proportionmat = accuracy_proportionmat;
data.accuracy_observe_num = accuracy_observe_num;
data.accuracy_observe_total = accuracy_observe_total;

end