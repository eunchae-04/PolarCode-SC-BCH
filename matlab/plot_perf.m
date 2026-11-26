% plot_perf.m
% polar_flip_sim 이 만든 perf_N*_R*.csv 를 읽어서 그래프를 그린다.
%   그림 1: FER / BER  (SC, SC+BCH, SC+패턴뒤집기, SC-Flip(시작 인덱스), SC-Flip(|LLR| 하위 T), SCL)
%   그림 2: SC의 "첫 번째 오류 비트"가 시작 인덱스 집합 / |LLR| 하위 T 안에 있는 비율,
%           SC가 틀린 프레임에서 "틀린 정보 비트 수" 분포
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
       0.4940 0.1840 0.5560;     % SC + BCH
       0.8500 0.3250 0.0980;     % SC + 패턴 뒤집기
       0.4660 0.6740 0.1880;     % SC-Flip (시작 인덱스)
       0.3010 0.7450 0.9330;     % SC-Flip (|LLR| 하위 T)
       0.6350 0.0780 0.1840];    % SCL
mk  = {'o', 'd', 's', 'v', 'p', '^'};
nm  = {'SC', 'SC + BCH', 'SC + pattern flip', 'SC-Flip (start idx)', 'SC-Flip (T lowest |LLR|)', 'SCL'};

for f = 1:numel(files)
    T = readtable(fullfile(files(f).folder, files(f).name));
    tok = regexp(files(f).name, 'N(\d+)_R(\d+)_(\d+)', 'tokens', 'once');
    N = str2double(tok{1}); Rn = tok{2}; Rd = tok{3};
    Kmsg = round(N * str2double(Rn) / str2double(Rd));   % 메시지 비트 수 (BER 상한 계산용)
    x = T.ebn0_db; fr = T.frames;

    % ---- 그림 1: FER / BER ----
    figure('Name', files(f).name, 'Position', [100 100 1000 420]);

    subplot(1, 2, 1);
    Y = {T.sc_fer, T.bch_fer, T.flip_fer, T.sfs_fer, T.sfl_fer, T.scl_fer};
    for c = 1:6, plotCurve(x, Y{c}, fr, 1, mk{c}, col(c, :), nm{c}); end
    grid on; xlabel('E_b/N_0 (dB)'); ylabel('FER');
    title(sprintf('FER, N = %d, R = %s/%s, CRC-24', N, Rn, Rd));
    addBoundLegend();

    subplot(1, 2, 2);
    Y = {T.sc_ber, T.bch_ber, T.flip_ber, T.sfs_ber, T.sfl_ber, T.scl_ber};
    for c = 1:6, plotCurve(x, Y{c}, fr, Kmsg, mk{c}, col(c, :), nm{c}); end
    grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
    title(sprintf('BER, N = %d, R = %s/%s, CRC-24', N, Rn, Rd));
    addBoundLegend();

    saveas(gcf, fullfile(resDir, strrep(files(f).name, '.csv', '.png')));

    % ---- 그림 2: SC 오류의 모양, BCH 정정 효과 ----
    figure('Name', [files(f).name ' (diagnostic)'], 'Position', [120 120 1000 420]);

    subplot(1, 2, 1);
    H = [T.wrong1 T.wrong2 T.wrong3_4 T.wrong5_8 T.wrong9p];
    Pr = H ./ max(sum(H, 2), 1) * 100;
    bar(x, Pr, 'stacked'); grid on;
    xlabel('E_b/N_0 (dB)'); ylabel('% of SC error frames');
    title('Wrong info bits per SC error frame');
    legend({'1 bit', '2 bits', '3-4 bits', '5-8 bits', '9+ bits'}, 'Location', 'eastoutside');
    ylim([0 100]);

    subplot(1, 2, 2);
    Ft = max(T.first_total, 1);
    bar(x, [T.first_in_start ./ Ft, T.first_in_llr ./ Ft] * 100); grid on;
    xlabel('E_b/N_0 (dB)'); ylabel('% of SC error frames'); ylim([0 100]);
    title('First wrong bit lies in candidate set');
    legend({'start-index set (d_{min} graph)', sprintf('T=%d lowest |LLR|', T.T(1))}, 'Location', 'northwest');

    saveas(gcf, fullfile(resDir, strrep(files(f).name, '.csv', '_diag.png')));
end

% 오류가 관측된 점은 실선으로 잇고, 0개인 점은 상한(3/전송 비트 수)을 v 마커로 표시
% y < 0 (해당 디코더를 사용하지 않음, 예: BCH 생략)이면 그 곡선은 그리지 않는다
function plotCurve(x, y, nfr, nbits, marker, color, name)
    if all(y < 0), return; end
    ok = (y >= 0);
    zero = ok & (y == 0);
    meas = ok & (y > 0);
    semilogy(x(meas), y(meas), ['-' marker], 'Color', color, 'LineWidth', 1.5, ...
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
