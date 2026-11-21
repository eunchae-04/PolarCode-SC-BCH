% plot_dmin.m
% polar_dmin 이 만든 dmin_N*_R*.csv 를 모두 읽어서
% x축: information bit index, y축: # of minimum Hamming weight codewords 막대그래프로 그린다.
clear; close all;

files = dir('dmin_N*_R*.csv');
for f = 1:numel(files)
    T   = readtable(files(f).name);
    idx = T.info_bit_index;
    cnt = T.num_min_weight_codewords;
    dmin = T.d_min(1);

    tok = regexp(files(f).name, 'N(\d+)_R(\d+)_(\d+)', 'tokens', 'once');
    N = str2double(tok{1}); Rn = tok{2}; Rd = tok{3};

    figure('Name', files(f).name);
    bar(1:numel(idx), cnt, 0.6);
    set(gca, 'XTick', 1:numel(idx), 'XTickLabel', string(idx));
    xlabel('Information bit index');
    ylabel('Number of minimum Hamming weight codewords');
    title(sprintf('N = %d, R = %s/%s,  d_{min} = %d,  A_{dmin} = %d', ...
                  N, Rn, Rd, dmin, sum(cnt)));
    grid on;
    saveas(gcf, strrep(files(f).name, '.csv', '.png'));
end
