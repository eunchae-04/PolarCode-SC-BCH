% plot_perf.m
% polar_flip_sim 이 만든 perf_N*_R*.csv 를 읽어서 FER / BER 그래프를 그린다.
%   SC, SC + 뒤집기(제안), SCL 세 곡선을 비교한다.
% 폴더 구조: <root>/matlab (이 파일), <root>/results (csv 읽기 + 그래프 저장)
clear; close all;

thisDir = fileparts(mfilename('fullpath'));
resDir  = fullfile(thisDir, '..', 'results');
files = dir(fullfile(resDir, 'perf_N*_R*.csv'));
if isempty(files)
    error('perf csv를 찾을 수 없습니다: %s (먼저 src_c에서 ./polar_flip_sim 을 실행하세요)', resDir);
end

for f = 1:numel(files)
    T = readtable(fullfile(files(f).folder, files(f).name));
    tok = regexp(files(f).name, 'N(\d+)_R(\d+)_(\d+)', 'tokens', 'once');
    N = str2double(tok{1}); Rn = tok{2}; Rd = tok{3};

    % FER = 0 인 점은 로그 축에 그릴 수 없으므로 NaN 처리 (그 SNR에서는 오류가 관측되지 않음)
    z = @(v) replaceZero(v);

    figure('Name', files(f).name, 'Position', [100 100 1000 420]);

    subplot(1, 2, 1);
    semilogy(T.ebn0_db, z(T.sc_fer),   '-o', 'LineWidth', 1.5); hold on;
    semilogy(T.ebn0_db, z(T.flip_fer), '-s', 'LineWidth', 1.5);
    semilogy(T.ebn0_db, z(T.scl_fer),  '-^', 'LineWidth', 1.5);
    grid on; xlabel('E_b/N_0 (dB)'); ylabel('FER');
    title(sprintf('FER, N = %d, R = %s/%s, CRC-24', N, Rn, Rd));
    legend('SC', 'SC + flip (proposed)', 'SCL', 'Location', 'southwest');

    subplot(1, 2, 2);
    semilogy(T.ebn0_db, z(T.sc_ber),   '-o', 'LineWidth', 1.5); hold on;
    semilogy(T.ebn0_db, z(T.flip_ber), '-s', 'LineWidth', 1.5);
    semilogy(T.ebn0_db, z(T.scl_ber),  '-^', 'LineWidth', 1.5);
    grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
    title(sprintf('BER, N = %d, R = %s/%s, CRC-24', N, Rn, Rd));
    legend('SC', 'SC + flip (proposed)', 'SCL', 'Location', 'southwest');

    saveas(gcf, fullfile(resDir, strrep(files(f).name, '.csv', '.png')));
end

function v = replaceZero(v)
    v(v == 0) = NaN;
end
