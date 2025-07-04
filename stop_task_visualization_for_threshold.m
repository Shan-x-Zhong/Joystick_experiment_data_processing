clear;
clc;
close all;

% Interactive visualization tool for Stop Task trial classification
% Shows X-Y trajectory, velocity, and reward analysis for each stop trial
% Press SPACE to advance through trials

baseDir = '/Users/shanzhong/Desktop/Mooncell/Angelica Cage/Postdoc/Action Regulation OR';
BASELINE_WINDOW_MS = 200;

patientFolders = dir(fullfile(baseDir, 'P*'));

for p = 1:length(patientFolders)
    patientDir = fullfile(baseDir, patientFolders(p).name);
    if ~patientFolders(p).isdir || startsWith(patientFolders(p).name, '.'), continue; end
    patientID = patientFolders(p).name;
    fprintf('\n===== Processing STOP TASK for patient: %s =====\n', patientID);

    % Calculate patient-specific thresholds
    fprintf('  Calculating thresholds for patient %s (Stop Task Only)...\n', patientID);
    try
        behavioralStopFiles = dir(fullfile(patientDir, '*stop1_behav_exp*.mat'));
        if ~isempty(behavioralStopFiles), isAppleDouble = startsWith({behavioralStopFiles.name}, '._'); behavioralStopFiles = behavioralStopFiles(~isAppleDouble); end
        if isempty(behavioralStopFiles), throw(MException('','Behavioral stop file not found.')); end
        load(fullfile(behavioralStopFiles(1).folder, behavioralStopFiles(1).name), 'stop_data');

        smoothedFiles = dir(fullfile(patientDir, 'processed', [patientID '_stop_sm*.mat']));
        if isempty(smoothedFiles), throw(MException('','Smoothed stop file not found.')); end
        load(fullfile(smoothedFiles(1).folder, smoothedFiles(1).name), 'joystickdataStop_sm_velocities');
        
        stop_baseline_noise_std = [];
        stop_trial_displacements = [];
        
        num_trials_to_process = min(size(stop_data,1), length(joystickdataStop_sm_velocities));
        for i = 1:num_trials_to_process
            if isempty(joystickdataStop_sm_velocities{1,i})||isempty(joystickdataStop_sm_velocities{1,i}{1,2}), continue; end
            x_pos_thresh = joystickdataStop_sm_velocities{1,i}{1,2};
            time_thresh = joystickdataStop_sm_velocities{1,i}{1,1};
            baseline_indices = time_thresh <= (time_thresh(1) + BASELINE_WINDOW_MS/1000);
            stop_baseline_noise_std(end+1) = std(x_pos_thresh(baseline_indices));
            
            if eq(nominal(stop_data.trial_type(i)), 'stop')
                stop_trial_displacements(end+1) = max(abs(x_pos_thresh - x_pos_thresh(1)));
            end
        end
        
        % Proactive threshold: 3x baseline noise
        proactive_threshold = 3 * mean(stop_baseline_noise_std, 'omitnan');
        
        % Displacement threshold: Otsu's method on stop trials
        data_for_otsu = stop_trial_displacements(stop_trial_displacements > 0);
        if length(unique(data_for_otsu)) > 2
            displacement_threshold = multithresh(data_for_otsu, 1);
        else
            displacement_threshold = mean(data_for_otsu);
        end
        fprintf('    -> Stop Task Proactive Threshold: %.0f | Stop Displacement Threshold: %.0f\n', proactive_threshold, displacement_threshold);
        
    catch ME
        warning('Could not calculate thresholds for patient %s. Skipping. Error: %s', patientID, ME.message);
        continue;
    end
    
    % Iterate through trials for visualization
    for trial_idx = 1:min(size(stop_data, 1), length(joystickdataStop_sm_velocities))
        
        if ~eq(nominal(stop_data(trial_idx, 5).trial_type), 'stop'), continue; end
        
        fprintf('\n--- Patient %s - Displaying Stop Trial %d ---\n', patientID, trial_idx);
        
        % Extract trial data
        if isempty(joystickdataStop_sm_velocities{1, trial_idx}), continue; end
        time = joystickdataStop_sm_velocities{1, trial_idx}{1, 1};
        x_position = joystickdataStop_sm_velocities{1, trial_idx}{1, 2};
        y_position = joystickdataStop_sm_velocities{1, trial_idx}{1, 3};
        stop_signal_time = double(stop_data.ssd(trial_idx));
        if isempty(time) || isempty(x_position), continue; end
        
        % Classify trial and calculate SSRT
        peak_disp = max(abs(x_position - x_position(1)));
        initial_classification = classifyStopTrial(peak_disp, proactive_threshold, displacement_threshold);
        final_classification = initial_classification;
        
        stop_onset_time = NaN; SSRT = NaN; reward = []; cumulativeReward = [];
        
        % Calculate SSRT for reactive and failed stops
        if strcmp(initial_classification, 'Correct (Reactive)') || strcmp(initial_classification, 'Failed Stop')
            [stop_onset_time, SSRT, reward, cumulativeReward] = findStopPoint_enhanced(time, x_position, stop_signal_time, initial_classification);
        end
        
        % Re-assign reactive trials without detectable SSRT to proactive
        if strcmp(initial_classification, 'Correct (Reactive)')
            if isnan(SSRT)
                final_classification = 'Correct (Proactive)';
                fprintf('    -> Re-assigned to Proactive: No valid SSRT detected\n');
            end
        end

        % Display classification results
        if strcmp(final_classification, 'Correct (Reactive)')
            if ~isnan(SSRT)
                fprintf('  Classification: Reactive Stop. SSRT: %.3f s\n', SSRT);
            else
                fprintf('  Classification: Reactive Stop. SSRT calculation failed.\n');
            end
        elseif strcmp(final_classification, 'Failed Stop')
            if ~isnan(SSRT)
                fprintf('  Classification: Failed Stop. SSRT: %.3f s\n', SSRT);
            else
                fprintf('  Classification: Failed Stop. SSRT calculation failed.\n');
            end
        elseif strcmp(final_classification, 'Correct (Proactive)')
            fprintf('  Classification: Proactive Stop.\n');
        end
        
        if strcmp(final_classification, 'Correct (Proactive)'), stop_onset_time = NaN; end

        % Create visualization figure
        fig = figure('Position', [100, 100, 1800, 500]);
        
        % Subplot 1: X-Y trajectory with thresholds
        ax1 = subplot(1, 3, 1);
        hold on;
        start_x = x_position(1);
        y_lims = [min(y_position) - 1000, max(y_position) + 1000];
        zone_x = [start_x - proactive_threshold, start_x + proactive_threshold, start_x + proactive_threshold, start_x - proactive_threshold];
        patch(zone_x, [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Zone of Inactivity (Proactive)');
        xline(start_x - displacement_threshold, '--', 'Color', [0 0.8 0.8], 'LineWidth', 1.5, 'DisplayName', 'Displacement Threshold');
        xline(start_x + displacement_threshold, '--', 'Color', [0 0.8 0.8], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(x_position, y_position, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Trajectory');
        plot(x_position(1), y_position(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Start');
        plot(x_position(end), y_position(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'End');
        x_at_stop_signal = interp1(time, x_position, stop_signal_time, 'linear', 'extrap');
        y_at_stop_signal = interp1(time, y_position, stop_signal_time, 'linear', 'extrap');
        plot(x_at_stop_signal, y_at_stop_signal, 'ms', 'MarkerFaceColor', 'm', 'MarkerSize', 10, 'DisplayName', 'Stop Signal (SSD)');
        
        if ~isnan(stop_onset_time)
            x_at_onset = interp1(time, x_position, stop_onset_time, 'linear', 'extrap');
            y_at_onset = interp1(time, y_position, stop_onset_time, 'linear', 'extrap');
            
            if strcmp(final_classification, 'Correct (Reactive)')
                plot(x_at_onset, y_at_onset, 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Detected Stop Onset (Success)');
            else  % Failed Stop
                plot(x_at_onset, y_at_onset, 'p', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.4 0.4], 'MarkerEdgeColor', 'k', 'DisplayName', 'Detected Stop Onset (Failed)');
            end
        end
        
        hold off; title('X-Y Trajectory with Thresholds'); xlabel('X Position'); ylabel('Y Position');
        legend('show', 'Location', 'best'); axis equal; grid on;

        % Subplot 2: X velocity
        x_velocity = gradient(x_position(:)) ./ gradient(time(:));
        subplot(1, 3, 2);
        plot(time, x_velocity, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'DisplayName', 'X Velocity'); hold on;
        xline(stop_signal_time, 'g--', 'LineWidth', 2, 'DisplayName', 'Stop Signal (SSD)');
        
        if ~isnan(stop_onset_time)
            if strcmp(final_classification, 'Correct (Reactive)')
                xline(stop_onset_time, 'r--', 'LineWidth', 2, 'DisplayName','Detected Stop Onset (Success)');
            else  % Failed Stop
                xline(stop_onset_time, '--', 'Color', [1 0.4 0.4], 'LineWidth', 2, 'DisplayName','Detected Stop Onset (Failed)');
            end
        end
        
        hold off; title('X Velocity vs. Time'); xlabel('Time (s)'); ylabel('X Velocity');
        legend('show', 'Location', 'best'); grid on;
        
        % Subplot 3: Reward analysis
        ax_reward = subplot(1, 3, 3);
        if ~isempty(reward) && (strcmp(final_classification, 'Correct (Reactive)') || strcmp(final_classification, 'Failed Stop'))
            hold on;
            time_for_reward = time(1:length(reward));
            yyaxis(ax_reward, 'left'); 
            plot(ax_reward, time_for_reward, reward, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Instant Reward');
            ylabel(ax_reward, 'Instantaneous Reward'); ax_reward.YAxis(1).Color = 'k'; yline(0, 'k--', 'HandleVisibility', 'off');
            
            yyaxis(ax_reward, 'right'); 
            plot(ax_reward, time_for_reward, cumulativeReward, 'LineWidth', 2.5, 'Color', [0.9290 0.6940 0.1250], 'DisplayName', 'Cumulative Reward');
            ylabel(ax_reward, 'Cumulative Reward'); ax_reward.YAxis(2).Color = [0.9290 0.6940 0.1250];
            
            xline(ax_reward, stop_signal_time, 'g--', 'LineWidth', 2, 'DisplayName','Stop Signal');
            
            if ~isnan(stop_onset_time)
                if strcmp(final_classification, 'Correct (Reactive)')
                    xline(ax_reward, stop_onset_time, 'r--', 'LineWidth', 2, 'DisplayName','Detected Stop Onset (Success)');
                else  % Failed Stop
                    xline(ax_reward, stop_onset_time, '--', 'Color', [1 0.4 0.4], 'LineWidth', 2, 'DisplayName','Detected Stop Onset (Failed)');
                end
                reward_at_onset = interp1(time_for_reward, cumulativeReward, stop_onset_time, 'linear');
                
                if strcmp(final_classification, 'Correct (Reactive)')
                    plot(ax_reward, stop_onset_time, reward_at_onset, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'Onset Point (Success)');
                else  % Failed Stop
                    plot(ax_reward, stop_onset_time, reward_at_onset, 'p', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.4 0.4], 'MarkerEdgeColor', 'k', 'DisplayName', 'Onset Point (Failed)');
                end
            end
            
            grid on; legend('show', 'Location', 'best'); hold off;
        else
            text(0.5, 0.5, 'No Onset to Detect for this Trial', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'blue');
            set(gca, 'XTick', [], 'YTick', []); box on;
        end
        title('Reward Analysis'); xlabel('Time (s)');
        
        % Figure title with SSRT
        ssrt_str = 'N/A'; 
        if ~isnan(SSRT), ssrt_str = sprintf('%.3f s', SSRT); end
        sgtitle(sprintf('%s - Stop Trial %d - Classification: %s (SSRT: %s)\nPRESS SPACE TO CONTINUE', ...
            patientID, trial_idx, final_classification, ssrt_str), 'FontSize', 16, 'FontWeight', 'bold');

        % Wait for space bar
        k = 0; is_fig_open = true;
        while k ~= 32
            try 
                if ~ishandle(fig), is_fig_open = false; break; end
                waitforbuttonpress; 
                if ~ishandle(fig), is_fig_open = false; break; end
                k = get(fig, 'CurrentCharacter'); 
                if isempty(k), k = 0; end
            catch
                is_fig_open = false; break; 
            end
        end
        if ishandle(fig), close(fig); end
        if ~is_fig_open, fprintf('Figure was closed. Exiting loop for this patient.\n'); break; end
    end
end
disp('--- Stop task display review completed. ---');

%% Helper Functions

function classification_str = classifyStopTrial(peak_disp, proactive_thresh, displacement_thresh)
    if peak_disp < proactive_thresh, classification_str = 'Correct (Proactive)';
    elseif peak_disp < displacement_thresh, classification_str = 'Correct (Reactive)';
    else, classification_str = 'Failed Stop';
    end
end

function [stopTime, ssrt, reward, cumulativeReward] = findStopPoint_enhanced(time, x_pos, stop_signal_time, classification_str)
    % Enhanced stop onset detection for both successful and failed stops
    stopTime = NaN; ssrt = NaN; reward = []; cumulativeReward = [];
    
    if strcmp(classification_str, 'Failed Stop')
        HORIZONTAL_DISPLACEMENT_THRESHOLD = 25;
    else
        HORIZONTAL_DISPLACEMENT_THRESHOLD = 50;
    end
    
    time = time(:); x_pos = x_pos(:);
    x_displacement = abs(x_pos - x_pos(1));
    
    gate_open_idx = find(x_displacement > HORIZONTAL_DISPLACEMENT_THRESHOLD, 1, 'first');
    if isempty(gate_open_idx), return; end
    
    vx = gradient(x_pos) ./ gradient(time); 
    vx(isnan(vx)) = 0;
    speed = abs(vx);
    
    if strcmp(classification_str, 'Failed Stop')
        onset_threshold = 0.05 * max(speed(gate_open_idx:end));
    else
        onset_threshold = 0.10 * max(speed(gate_open_idx:end));
    end
    
    onset_idx_relative = find(speed(gate_open_idx:end) > onset_threshold, 1, 'first');
    if isempty(onset_idx_relative), onset_idx = gate_open_idx; 
    else, onset_idx = gate_open_idx + onset_idx_relative - 1; end
    
    dir_window_end = min(length(vx), onset_idx + round(0.150 / mean(diff(time))));
    initial_vx_sign = sign(mean(vx(onset_idx:dir_window_end)));
    if initial_vx_sign == 0, initial_vx_sign = 1; end
    
    reward = vx * initial_vx_sign;
    cumulativeReward = cumsum(reward);
    
    if strcmp(classification_str, 'Failed Stop')
        signal_idx = find(time >= stop_signal_time, 1, 'first');
        if isempty(signal_idx), return; end
        
        [peak_vals, peak_indices_relative] = findpeaks(cumulativeReward(signal_idx:end));
        if ~isempty(peak_vals)
            [~, max_peak_idx] = max(peak_vals);
            peakIdx = signal_idx + peak_indices_relative(max_peak_idx) - 1;
        else
            [~, max_idx_relative] = max(cumulativeReward(signal_idx:end));
            peakIdx = signal_idx + max_idx_relative - 1;
        end
    else
        [~, peakIdx] = max(cumulativeReward);
    end
    
    detected_stop_time = time(peakIdx);
    
    if detected_stop_time < stop_signal_time, return; end 
    
    calculated_ssrt = detected_stop_time - stop_signal_time;
    
    MIN_RT = 0.100; 
    if strcmp(classification_str, 'Failed Stop')
        MAX_RT = 1.500;
    else
        MAX_RT = 1.000;
    end
    
    if calculated_ssrt >= MIN_RT && calculated_ssrt <= MAX_RT
        stopTime = detected_stop_time;
        ssrt = calculated_ssrt;
    end
end