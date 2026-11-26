function plot_compare_results(resultDir)
% Plot BER/FER files written by polar_code_sc_scl_cascl_compare.c
% Usage: plot_compare_results('result_N256_R033')

if nargin < 1 || isempty(resultDir)
    resultDir = 'result';
end

if ~isfolder(resultDir)
    error('Result directory not found: %s', resultDir);
end

files = dir(fullfile(resultDir, 'result_*.txt'));
if isempty(files)
    error('No result_*.txt files found in %s', resultDir);
end

figure('Name', ['FER - ' resultDir]);
hold on; grid on; box on;
ferPlotted = false;

figure('Name', ['BER - ' resultDir]);
hold on; grid on; box on;
berPlotted = false;

for k = 1:numel(files)
    filePath = fullfile(files(k).folder, files(k).name);
    M = readmatrix(filePath, 'CommentStyle', '#');
    if isempty(M) || size(M,2) < 4
        warning('Skipping empty or malformed file: %s', filePath);
        continue;
    end

    ebno = M(:,1);
    ber = M(:,2);
    fer = M(:,3);
    frames = M(:,4);
    [ebno, order] = sort(ebno);
    ber = ber(order); fer = fer(order); frames = frames(order);
    scheme = erase(erase(files(k).name, 'result_'), '.txt');

    figure(findobj('Type','figure','Name',['FER - ' resultDir]));
    positive = fer > 0;
    if any(positive)
        semilogy(ebno(positive), fer(positive), '-o', ...
            'LineWidth', 1.3, 'DisplayName', scheme);
        ferPlotted = true;
    end
    zero = ~positive;
    if any(zero)
        semilogy(ebno(zero), 1 ./ frames(zero), 'v', ...
            'HandleVisibility', 'off');
        ferPlotted = true;
    end

    figure(findobj('Type','figure','Name',['BER - ' resultDir]));
    positive = ber > 0;
    if any(positive)
        semilogy(ebno(positive), ber(positive), '-s', ...
            'LineWidth', 1.3, 'DisplayName', scheme);
        berPlotted = true;
    end
    zero = ~positive;
    if any(zero)
        % A zero-BER point is shown at the approximate one-error resolution
        % only when total information-bit count is not present in the file.
        % The file stores frames, not K, so use 1/frames as a visual marker.
        semilogy(ebno(zero), 1 ./ frames(zero), 'v', ...
            'HandleVisibility', 'off');
        berPlotted = true;
    end
end

figure(findobj('Type','figure','Name',['FER - ' resultDir]));
xlabel('E_b/N_0 (dB)'); ylabel('FER');
title(['FER comparison: ' resultDir]);
if ferPlotted, legend('Location','southwest'); end

figure(findobj('Type','figure','Name',['BER - ' resultDir]));
xlabel('E_b/N_0 (dB)'); ylabel('BER');
title(['BER comparison: ' resultDir]);
if berPlotted, legend('Location','southwest'); end

fprintf('Finished plotting result files in %s.\n', resultDir);
end
