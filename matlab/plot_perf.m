% plot_perf.m
% polar_flip_sim 이 만든 perf_N*_R*.csv 를 읽어서 FER / BER 그래프를 그린다.
%   SC, SC + 뒤집기(제안), SCL 세 곡선을 비교한다.
% 폴더 구조: <root>/matlab (이 파일), <root>/results (csv 읽기 + 그래프 저장)
%
% 오류가 0개 관측된 점은 FER = 0 이 아니라 "관측한 프레임 수로는 못 본 것"이므로,
% 로그 축에서 사라지지 않도록 상한(오류 0개일 때의 95% 상한 = 3 / 프레임 수)을
% 속이 빈 아래쪽 삼각형(v)으로 표시한다. 실제 FER은 그 값보다 작다.
clear; close all;

thisDir = fileparts(mfilename('fullpath'));
resDir  = fullfile(thisDir, '..', 'results');
files = dir(fullfile(resDir, 'perf_N*_R*.csv'));
if isempty(files)
    error('perf csv를 찾을 수 없습니다: %s (먼저 src_c에서 ./polar_flip_sim 을 실행하세요)', resDir);
end

col = [0.0000 0.4470 0.7410;     % SC
    0.8500 0.3250 0.0980;     % SC + flip
    0.9290 0.6940 0.1250];    % SCL
mk  = {'o', 's', '^'};
nm  = {'SC', 'SC + flip (proposed)', 'SCL'};

for f = 1:numel(files)
    T = readtable(fullfile(files(f).folder, files(f).name));
    tok = regexp(files(f).name, 'N(\d+)_R(\d+)_(\d+)', 'tokens', 'once');
    N = str2double(tok{1}); Rn = tok{2}; Rd = tok{3};
    Kmsg = round(N * str2double(Rn) / str2double(Rd));   % 메시지 비트 수 (BER 상한 계산용)
    x = T.ebn0_db; fr = T.frames;

    figure('Name', files(f).name, 'Position', [100 100 1000 420]);

    subplot(1, 2, 1);
    Y = {T.sc_fer, T.flip_fer, T.scl_fer};
    for c = 1:3, plotCurve(x, Y{c}, fr, 1, mk{c}, col(c, :), nm{c}); end
    grid on; xlabel('E_b/N_0 (dB)'); ylabel('FER');
    title(sprintf('FER, N = %d, R = %s/%s, CRC-24', N, Rn, Rd));
    addBoundLegend();

    subplot(1, 2, 2);
    Y = {T.sc_ber, T.flip_ber, T.scl_ber};
    for c = 1:3, plotCurve(x, Y{c}, fr, Kmsg, mk{c}, col(c, :), nm{c}); end
    grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
    title(sprintf('BER, N = %d, R = %s/%s, CRC-24', N, Rn, Rd));
    addBoundLegend();

    saveas(gcf, fullfile(resDir, strrep(files(f).name, '.csv', '.png')));
end

% 오류가 관측된 점은 실선으로 잇고, 0개인 점은 상한(3/전송 비트 수)을 v 마커로 표시
function plotCurve(x, y, nfr, nbits, marker, color, name)
zero = (y == 0);
semilogy(x(~zero), y(~zero), ['-' marker], 'Color', color, 'LineWidth', 1.5, ...
    'DisplayName', name); hold on;
if any(zero)
    semilogy(x(zero), 3 ./ (nfr(zero) * nbits), 'v', 'Color', color, ...
        'MarkerSize', 8, 'LineWidth', 1.3, 'HandleVisibility', 'off');
end
end

function addBoundLegend()
semilogy(NaN, NaN, 'vk', 'MarkerSize', 8, 'LineWidth', 1.3, ...
    'DisplayName', 'no error observed (95% upper bound)');
legend('Location', 'southwest');
end