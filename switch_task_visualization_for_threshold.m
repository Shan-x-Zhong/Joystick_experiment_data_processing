clear;
clc;
close all;

% Interactive visualization tool for Switch Task trial classification
% Shows X-Y trajectory, velocities, and reward analysis for each switch trial
% Uses separate thresholds: switch task proactive, stop task displacement
% Press SPACE to advance through trials

baseDir = '/Users/shanzhong/Desktop/Mooncell/Angelica Cage/Postdoc/Action Regulation OR'; 
BASELINE_WINDOW_MS = 200;
Y_POSITION_LIMIT = 25000;

patientFolders = dir(fullfile(baseDir, 'P*'));

for p = 1:length(patientFolders)
    patientDir = fullfile(baseDir, patientFolders(p).name);
    if ~patientFolders(p).isdir || startsWith(patientFolders(p).name, '.'), continue; end
    patientID = patientFolders(p).name;
    fprintf('\n===== Processing SWITCH TASK for patient: %s =====\n', patientID);

    % Calculate patient-specific thresholds with separate approach
    fprintf('  Calculating separate thresholds for patient %s...\n', patientID);
    try
        % Load behavioral and kinematic data
        stop_files = dir(fullfile(patientDir, '*stop1_behav_exp*.mat'));
        switch_files = dir(fullfile(patientDir, '*switch_behav_exp*.mat'));
        if isempty(stop_files) || isempty(switch_files), throw(MException('','Behavioral file(s) for thresholding not found.')); end
        load(fullfile(stop_files(1).folder, stop_files(1).name), 'stop_data');
        load(fullfile(switch_files(1).folder, switch_files(1).name), 'switch_data');
        stop_sm_files = dir(fullfile(patientDir, 'processed', [patientID '_stop_sm*.mat']));
        switch_sm_files = dir(fullfile(patientDir, 'processed', [patientID '_switch_sm*.mat']));
        if isempty(stop_sm_files) || isempty(switch_sm_files), throw(MException('','Smoothed file(s) for thresholding not found.')); end
        load(fullfile(stop_sm_files(1).folder, stop_sm_files(1).name), 'joystickdataStop_sm_velocities');
        load(fullfile(switch_sm_files(1).folder, switch_sm_files(1).name), 'joystickdataSwitch_sm_velocities');
        
        % Calculate separate baseline noise from each task
        stop_baseline_noise_std = [];
        switch_baseline_noise_std = [];
        stop_trial_displacements = [];
        
        % Stop task: baseline noise + displacement data
        for i = 1:min(size(stop_data,1), length(joystickdataStop_sm_velocities))
            if isempty(joystickdataStop_sm_velocities{1,i})||isempty(joystickdataStop_sm_velocities{1,i}{1,2}), continue; end
            x_pos = joystickdataStop_sm_velocities{1,i}{1,2}; time = joystickdataStop_sm_velocities{1,i}{1,1};
            stop_baseline_noise_std(end+1) = std(x_pos(time <= (time(1) + BASELINE_WINDOW_MS/1000)));
            if eq(nominal(stop_data.trial_type(i)), 'stop'), stop_trial_displacements(end+1) = max(abs(x_pos - x_pos(1))); end
        end
        
        % Switch task: baseline noise only
        for i = 1:min(size(switch_data,1), length(joystickdataSwitch_sm_velocities))
            if isempty(joystickdataSwitch_sm_velocities{1,i})||isempty(joystickdataSwitch_sm_velocities{1,i}{1,2}), continue; end
            x_pos = joystickdataSwitch_sm_velocities{1,i}{1,2}; time = joystickdataSwitch_sm_velocities{1,i}{1,1};
            switch_baseline_noise_std(end+1) = std(x_pos(time <= (time(1) + BASELINE_WINDOW_MS/1000)));
        end
        
        if isempty(stop_trial_displacements), throw(MException('','No "stop" type trials found to calculate displacement threshold.')); end
        
        % Proactive threshold: 3x switch task baseline noise
        proactive_threshold = 3 * mean(switch_baseline_noise_std, 'omitnan');
        
        % Displacement threshold: Otsu's method on stop trial displacements
        data_for_otsu = stop_trial_displacements(stop_trial_displacements > 0);
        if length(unique(data_for_otsu)) > 2, displacement_threshold = multithresh(data_for_otsu, 1);
        else, displacement_threshold = mean(data_for_otsu); end
        
        fprintf('    -> Switch Task Proactive Threshold: %.0f | Stop Displacement Threshold: %.0f\n', proactive_threshold, displacement_threshold);
    catch ME, warning('Could not calculate thresholds for patient %s. Skipping. Error: %s', patientID, ME.message); continue; end
    
    % Iterate through switch trials for visualization
    for trial_idx = 1:min(size(switch_data, 1), length(joystickdataSwitch_sm_velocities))
        
        if ~eq(nominal(switch_data.trial_type(trial_idx)), 'switch'), continue; end
        
        % Extract trial data
        time = joystickdataSwitch_sm_velocities{1,trial_idx}{1,1}; 
        x_position = joystickdataSwitch_sm_velocities{1,trial_idx}{1,2}; 
        y_position = joystickdataSwitch_sm_velocities{1,trial_idx}{1,3};
        if isempty(time) || length(time) < 10, continue; end
        switch_signal_time = double(switch_data.swsd(trial_idx));
        trial_direction = char(switch_data.direction(trial_idx));

        % Classify trial and calculate switch onset
        peak_disp = max(abs(x_position - x_position(1)));
        initial_classification = classifySwitchTrial(peak_disp, proactive_threshold, displacement_threshold);
        final_classification = initial_classification;
        
        switch_onset_time = NaN; reward = []; cumulativeReward = [];
        
        % Calculate switch onset for reactive and failed trials
        if strcmp(initial_classification, 'Correct (Reactive)') || strcmp(initial_classification, 'Failed Switch')
            [switch_onset_time, reward, cumulativeReward] = findSwitchPoint_hybrid(time, x_position, y_position, switch_signal_time, trial_direction, Y_POSITION_LIMIT);
        end
        
        % Re-assignment logic with spatial constraints
        if strcmp(initial_classification, 'Correct (Reactive)')
            if ~isnan(switch_onset_time)
                % Check X-axis constraint (movement must exceed proactive threshold)
                x_at_onset = interp1(time, x_position, switch_onset_time, 'linear', 'extrap');
                x_constraint_failed = abs(x_at_onset - x_position(1)) < proactive_threshold;
                
                % Check Y-axis constraint (vertical movement within limits)
                y_at_onset = interp1(time, y_position, switch_onset_time, 'linear', 'extrap');
                y_constraint_failed = abs(y_at_onset) > Y_POSITION_LIMIT;
                
                % Re-assign to Proactive if constraints fail
                if x_constraint_failed || y_constraint_failed
                    final_classification = 'Correct (Proactive)';
                    if x_constraint_failed && y_constraint_failed
                        fprintf('    -> Re-assigned to Proactive: X-position too close to start AND Y-position out of bounds\n');
                    elseif x_constraint_failed
                        fprintf('    -> Re-assigned to Proactive: X-position too close to start position\n');
                    else
                        fprintf('    -> Re-assigned to Proactive: Y-position out of bounds (%.0f)\n', y_at_onset);
                    end
                end
            else
                final_classification = 'Correct (Proactive)';
                fprintf('    -> Re-assigned to Proactive: No switch onset detected\n');
            end
        end

        % Display classification results
        fprintf('\n--- Patient %s - Displaying Switch Trial %d ---\n', patientID, trial_idx);
        RT_auto = NaN;
        
        if strcmp(final_classification, 'Correct (Reactive)')
            if ~isnan(switch_onset_time)
                RT_auto = switch_onset_time - switch_signal_time; 
                fprintf('  Classification: Reactive Switch. RT: %.3f s\n', RT_auto);
            else
                fprintf('  Classification: Reactive Switch. RT calculation failed.\n');
            end
        elseif strcmp(final_classification, 'Failed Switch')
            if ~isnan(switch_onset_time)
                RT_auto = switch_onset_time - switch_signal_time; 
                fprintf('  Classification: Failed Switch. RT: %.3f s\n', RT_auto);
            else
                fprintf('  Classification: Failed Switch. RT calculation failed.\n');
            end
        elseif strcmp(final_classification, 'Correct (Proactive)')
            fprintf('  Classification: Proactive Switch.\n');
        end
        
        if strcmp(final_classification, 'Correct (Proactive)'), switch_onset_time = NaN; end
        
        % Create visualization figure
        fig = figure('Position', [100, 100, 1600, 800]);
        
        % Subplot 1: X-Y trajectory with thresholds
        subplot(2,2,1); hold on;
        start_x = x_position(1); y_lims = [min(y_position) - 1000, max(y_position) + 1000];
        zone_x = [start_x - proactive_threshold, start_x + proactive_threshold, start_x + proactive_threshold, start_x - proactive_threshold];
        patch(zone_x, [y_lims(1) y_lims(1) y_lims(2) y_lims(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', 'Zone of Inactivity');
        xline(start_x - displacement_threshold, '--', 'Color', [0 0.8 0.8], 'LineWidth', 1.5, 'DisplayName', 'Displacement Threshold');
        xline(start_x + displacement_threshold, '--', 'Color', [0 0.8 0.8], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        % Y-axis constraint visualization
        yline(Y_POSITION_LIMIT, ':', 'Color', [1 0.5 0], 'LineWidth', 2, 'DisplayName', 'Y-Position Limit');
        yline(-Y_POSITION_LIMIT, ':', 'Color', [1 0.5 0], 'LineWidth', 2, 'HandleVisibility', 'off');
        
        plot(x_position, y_position, 'k-', 'LineWidth', 1, 'DisplayName', 'Trajectory');
        plot(x_position(1), y_position(1), 'go', 'MarkerFaceColor','g', 'MarkerSize', 8, 'DisplayName','Start');
        plot(x_position(end), y_position(end), 'ro','MarkerFaceColor','r', 'MarkerSize', 8, 'DisplayName','End');
        x_at_signal = interp1(time, x_position, switch_signal_time, 'linear', 'extrap'); 
        y_at_signal = interp1(time, y_position, switch_signal_time, 'linear', 'extrap');
        plot(x_at_signal, y_at_signal, 'ms', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'Switch Signal');
        
        if ~isnan(switch_onset_time)
            x_at_onset = interp1(time, x_position, switch_onset_time, 'linear', 'extrap'); 
            y_at_onset = interp1(time, y_position, switch_onset_time, 'linear', 'extrap'); 
            
            if strcmp(final_classification, 'Correct (Reactive)')
                plot(x_at_onset, y_at_onset, 'p', 'MarkerSize', 14, 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'MarkerEdgeColor', 'k', 'DisplayName', 'Detected Switch Onset (Success)');
            else  % Failed Switch
                plot(x_at_onset, y_at_onset, 'p', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.4 0.4], 'MarkerEdgeColor', 'k', 'DisplayName', 'Detected Switch Onset (Failed)');
            end
        end
        
        title('X-Y Trajectory with Thresholds'); xlabel('X Position'); ylabel('Y Position'); 
        legend('show', 'Location', 'best'); axis equal; grid on; set(gca, 'YDir', 'reverse'); hold off;

        % Subplot 2 & 3: X and Y velocities
        x_velocity = gradient(x_position) ./ gradient(time); y_velocity = gradient(y_position) ./ gradient(time);
        subplot(2,2,3); plot(time, x_velocity, 'b-'); hold on; 
        xline(switch_signal_time, 'g--', 'LineWidth', 2, 'DisplayName', 'Switch Signal'); 
        if ~isnan(switch_onset_time)
            if strcmp(final_classification, 'Correct (Reactive)')
                xline(switch_onset_time, 'r--', 'LineWidth', 2, 'DisplayName', 'Switch Onset (Success)'); 
            else
                xline(switch_onset_time, '--', 'Color', [1 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Switch Onset (Failed)'); 
            end
        end
        title('X Velocity'); xlabel('Time (s)'); ylabel('X Velocity'); legend('show'); grid on; hold off;
        
        subplot(2,2,2); plot(time, y_velocity, 'b-'); hold on; 
        xline(switch_signal_time, 'g--', 'LineWidth', 2, 'DisplayName', 'Switch Signal'); 
        if ~isnan(switch_onset_time)
            if strcmp(final_classification, 'Correct (Reactive)')
                xline(switch_onset_time, 'r--', 'LineWidth', 2, 'DisplayName', 'Switch Onset (Success)'); 
            else
                xline(switch_onset_time, '--', 'Color', [1 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Switch Onset (Failed)'); 
            end
        end
        title('Y Velocity'); xlabel('Time (s)'); ylabel('Y Velocity'); legend('show'); grid on; hold off;
        
        % Subplot 4: Reward analysis
        ax_reward = subplot(2,2,4);
        if ~isempty(reward) && (strcmp(final_classification, 'Correct (Reactive)') || strcmp(final_classification, 'Failed Switch'))
            hold on; time_for_reward = time(1:length(reward));
            yyaxis(ax_reward, 'left'); 
            plot(ax_reward, time_for_reward, reward, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Instant Reward'); 
            ylabel('Instant Reward'); ax_reward.YAxis(1).Color = 'k';
            
            yyaxis(ax_reward, 'right'); 
            plot(ax_reward, time_for_reward, cumulativeReward, 'LineWidth', 2.5, 'DisplayName', 'Cumulative Reward'); 
            ylabel('Cumulative Reward'); ax_reward.YAxis(2).Color = [0.9290 0.6940 0.1250];
            
            xline(ax_reward, switch_signal_time, 'g--', 'LineWidth', 2, 'DisplayName', 'Switch Signal');
            
            if ~isnan(switch_onset_time)
                if strcmp(final_classification, 'Correct (Reactive)')
                    xline(ax_reward, switch_onset_time, 'r--', 'LineWidth', 2, 'DisplayName', 'Switch Onset (Success)');
                else
                    xline(ax_reward, switch_onset_time, '--', 'Color', [1 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Switch Onset (Failed)');
                end
                reward_at_onset = interp1(time_for_reward, cumulativeReward, switch_onset_time, 'linear', 'extrap');
                
                if strcmp(final_classification, 'Correct (Reactive)')
                    plot(ax_reward, switch_onset_time, reward_at_onset, 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'DisplayName', 'Onset Point (Success)');
                else
                    plot(ax_reward, switch_onset_time, reward_at_onset, 'p', 'MarkerSize', 14, 'MarkerFaceColor', [1 0.4 0.4], 'MarkerEdgeColor', 'k', 'DisplayName', 'Onset Point (Failed)');
                end
            end
            
            legend('show', 'Location', 'best'); grid on; hold off;
        else
            text(0.5, 0.5, 'Proactive Trial - No Onset Analysis', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', 'blue'); 
            set(gca, 'XTick', [], 'YTick', []); box on;
        end
        title('Reward Analysis (Weighted Logic)'); xlabel('Time (s)');
        
        % Figure title with RT information
        rt_string_title = 'N/A'; 
        if ~isnan(RT_auto), rt_string_title = sprintf('%.3fs', RT_auto); end
        direction_str_formatted = strrep(trial_direction, '_', '-');
        sgtitle(sprintf('%s - Switch Trial %d (%s) - Classification: %s (RT: %s)\nPRESS SPACE', ...
            patientID, trial_idx, direction_str_formatted, final_classification, rt_string_title), 'FontSize', 14, 'FontWeight', 'bold');
        
        % Wait for space bar
        k = 0; valid_fig = true; 
        while k~=32
            try 
                if ~ishandle(fig), valid_fig=false; break; end
                waitforbuttonpress; 
                if ~ishandle(fig), valid_fig=false; break; end
                k=get(fig,'CurrentCharacter'); 
                if isempty(k), k=0; end
            catch 
                valid_fig=false; break; 
            end
        end
        if ishandle(fig), close(fig); end
        if ~valid_fig, break; end
    end
end
disp('Switch task display review completed.');

%% Helper Functions

function [switchTime, reward, cumulativeReward] = findSwitchPoint_hybrid(time, x_pos, y_pos, switch_signal_time, trial_direction, y_limit)
    % Enhanced switch onset detection for both successful and failed trials
    % Uses weighted reward function: favors horizontal movement, penalizes vertical drift
    reward = []; cumulativeReward = []; switchTime = NaN;
    time = time(:); x_pos = x_pos(:); y_pos = y_pos(:);
    
    vx = gradient(x_pos) ./ gradient(time); 
    vy = gradient(y_pos) ./ gradient(time);
    vx(isnan(vx)) = 0; vy(isnan(vy)) = 0;
    
    if contains(trial_direction, 'left'), initial_vx_sign = -1;
    elseif contains(trial_direction, 'right'), initial_vx_sign = 1;
    else, warning('Invalid horizontal direction: %s', trial_direction); return; end
    
    % Weighted reward: 3x horizontal movement - 1x vertical movement
    weight_x = 3.0; 
    weight_y = 1.0;
    reward = (weight_x * (vx * initial_vx_sign)) - (weight_y * abs(vy));
    
    cumulativeReward = cumsum(reward);
    time_for_reward = time;
    
    signal_idx = find(time_for_reward >= switch_signal_time, 1, 'first');
    if isempty(signal_idx), return; end
    
    [peak_vals, peak_indices_relative] = findpeaks(cumulativeReward(signal_idx:end));
    if isempty(peak_vals), return; end
    
    [~, sort_order] = sort(peak_vals, 'descend');
    sorted_peak_indices = signal_idx + peak_indices_relative(sort_order) - 1;
    
    MIN_REACTION_TIME = 0.100;
    MAX_REACTION_TIME = 1.500;
    
    for i = 1:length(sorted_peak_indices)
        candidate_idx = sorted_peak_indices(i);
        candidate_time = time_for_reward(candidate_idx);
        
        reaction_time = candidate_time - switch_signal_time;
        is_temporally_valid = reaction_time >= MIN_REACTION_TIME && reaction_time <= MAX_REACTION_TIME;
        
        y_at_candidate = interp1(time, y_pos, candidate_time, 'linear', 'extrap');
        is_spatially_valid = abs(y_at_candidate) <= y_limit * 2;
        
        if is_temporally_valid && is_spatially_valid
            switchTime = candidate_time; 
            return; 
        end
    end
end

function classification_str = classifySwitchTrial(peak_disp, proactive_thresh, displacement_thresh)
    % Three-tier classification based on peak displacement
    if peak_disp < proactive_thresh
        classification_str = 'Correct (Proactive)';
    elseif peak_disp < displacement_thresh
        classification_str = 'Correct (Reactive)';
    else
        classification_str = 'Failed Switch';
    end
end