clear;
clc;
close all;

% =======================================================================
% ACTION REGULATION TASK ANALYSIS
% Batch processing script for Stop and Switch task behavioral data
% Classifies trials as Proactive, Reactive, or Failed based on kinematics
% =======================================================================

baseDir = '/Users/shanzhong/Desktop/Mooncell/Angelica Cage/Postdoc/Action Regulation OR';
GO_TRIAL_RT_THRESHOLD = 2.5; % seconds
BASELINE_WINDOW_MS = 200; % Use first 200ms of each trial to measure noise
Y_POSITION_LIMIT = 25000; % Y-axis limit for switch task

% Initialize global counters for aggregated results across all patients
fprintf('--- Initializing Global Counters for Final Report ---\n');
fields = {'StopTask_Go', 'StopTask_Stop', 'SwitchTask_Go', 'SwitchTask_Switch'};
global_counts = struct();
for i = 1:length(fields)
    global_counts.(fields{i}).Correct_Proactive = 0;
    global_counts.(fields{i}).Correct_Reactive = 0;
    global_counts.(fields{i}).Correct = 0;
    global_counts.(fields{i}).Incorrect = 0;
end

% Main processing loop - iterate through all patient folders
patientFolders = dir(fullfile(baseDir, 'P*'));
for p = 1:length(patientFolders)
    patientDir = fullfile(baseDir, patientFolders(p).name);
    patientID = patientFolders(p).name;
    fprintf('\n\n===== Processing Patient: %s =====\n', patientID);

    % Initialize patient-specific counters
    patient_counts = global_counts;

    % ===================================================================
    % STOP TASK DATA PROCESSING
    % ===================================================================
    fprintf('  1. Processing Stop Task data...\n');
    stop_data = []; 
    joystickdataStop_sm_velocities = []; 
    stop_baseline_noise_std = []; 
    stop_trial_displacements = [];
    stop_kinematic_data = struct();
    
    try
        % Load behavioral data files
        stop_files = dir(fullfile(patientDir, '*stop1_behav_exp*.mat'));
        if ~isempty(stop_files)
            isAppleDouble = startsWith({stop_files.name}, '._'); % Filter out system files
            stop_files = stop_files(~isAppleDouble);
        end
        
        % Load processed kinematic data
        stop_sm_files = dir(fullfile(patientDir, 'processed', [patientID '_stop_sm*.mat']));
        
        if ~isempty(stop_files) && ~isempty(stop_sm_files)
            load(fullfile(stop_files(1).folder, stop_files(1).name), 'stop_data');
            load(fullfile(stop_sm_files(1).folder, stop_sm_files(1).name), 'joystickdataStop_sm_velocities');
            
            % Extract kinematic data for each trial
            types = nominal(stop_data.trial_type);
            for i = 1:min(length(types), length(joystickdataStop_sm_velocities))
                if isempty(joystickdataStop_sm_velocities{1,i}) || isempty(joystickdataStop_sm_velocities{1,i}{1,2})
                    continue; 
                end
                
                % Extract time and position data
                time_data = joystickdataStop_sm_velocities{1,i}{1,1};
                x_pos_data = joystickdataStop_sm_velocities{1,i}{1,2};
                y_pos_data = joystickdataStop_sm_velocities{1,i}{1,3};
                
                % Store kinematic data structure
                stop_kinematic_data(i).time = time_data;
                stop_kinematic_data(i).x_pos = x_pos_data;
                stop_kinematic_data(i).y_pos = y_pos_data;
                stop_kinematic_data(i).trial_type = char(types(i));
                
                % Calculate peak displacement and baseline noise
                peak_disp = max(abs(x_pos_data - x_pos_data(1)));
                stop_baseline_noise_std(end+1) = std(x_pos_data(time_data <= (time_data(1) + BASELINE_WINDOW_MS/1000)));
                
                % Collect displacement data for stop trials (for threshold calculation)
                if eq(types(i), 'stop')
                    stop_trial_displacements(end+1) = peak_disp;
                end
            end
        else
            fprintf('    -> Stop task files not found for patient %s\n', patientID);
        end
    catch ME
        warning('Could not process Stop Task data for patient %s: %s', patientID, ME.message);
    end

    % ===================================================================
    % SWITCH TASK DATA PROCESSING
    % ===================================================================
    fprintf('  2. Processing Switch Task data...\n');
    switch_data = []; 
    joystickdataSwitch_sm_velocities = [];
    switch_baseline_noise_std = [];
    switch_kinematic_data = struct();
    
    try
        % Load behavioral data files
        switch_files = dir(fullfile(patientDir, '*switch_behav_exp*.mat'));
        if ~isempty(switch_files)
            isAppleDouble = startsWith({switch_files.name}, '._'); % Filter out system files
            switch_files = switch_files(~isAppleDouble);
        end
        
        % Load processed kinematic data
        switch_sm_files = dir(fullfile(patientDir, 'processed', [patientID '_switch_sm*.mat']));
        
        if ~isempty(switch_files) && ~isempty(switch_sm_files)
            load(fullfile(switch_files(1).folder, switch_files(1).name), 'switch_data');
            load(fullfile(switch_sm_files(1).folder, switch_sm_files(1).name), 'joystickdataSwitch_sm_velocities');
            
            % Extract kinematic data for each trial
            types = nominal(switch_data.trial_type);
            for i = 1:min(length(types), length(joystickdataSwitch_sm_velocities))
                if isempty(joystickdataSwitch_sm_velocities{1,i}) || isempty(joystickdataSwitch_sm_velocities{1,i}{1,2})
                    continue; 
                end
                
                % Extract time and position data
                time_data = joystickdataSwitch_sm_velocities{1,i}{1,1};
                x_pos_data = joystickdataSwitch_sm_velocities{1,i}{1,2};
                y_pos_data = joystickdataSwitch_sm_velocities{1,i}{1,3};
                
                % Store kinematic data structure
                switch_kinematic_data(i).time = time_data;
                switch_kinematic_data(i).x_pos = x_pos_data;
                switch_kinematic_data(i).y_pos = y_pos_data;
                switch_kinematic_data(i).trial_type = char(types(i));
                
                % Calculate baseline noise for threshold setting
                switch_baseline_noise_std(end+1) = std(x_pos_data(time_data <= (time_data(1) + BASELINE_WINDOW_MS/1000)));
            end
        else
            fprintf('    -> Switch task files not found for patient %s\n', patientID);
        end
    catch ME
        warning('Could not process Switch Task data for patient %s: %s', patientID, ME.message);
    end

    % ===================================================================
    % THRESHOLD CALCULATION (Patient-specific)
    % ===================================================================
    fprintf('  3. Calculating separate thresholds for this patient...\n');
    
    % Check if we have sufficient data
    if isempty(stop_baseline_noise_std) || isempty(stop_trial_displacements)
        warning('Insufficient stop task data for patient %s', patientID);
        continue;
    end
    
    if isempty(switch_baseline_noise_std)
        warning('Insufficient switch task data for patient %s', patientID);
        continue;
    end

    % Calculate proactive thresholds (3x baseline noise)
    stop_proactive_threshold = 3 * mean(stop_baseline_noise_std, 'omitnan');
    switch_proactive_threshold = 3 * mean(switch_baseline_noise_std, 'omitnan');

    % Calculate displacement thresholds using Otsu's method
    data_for_stop_thresh = stop_trial_displacements(stop_trial_displacements > 0);
    if length(unique(data_for_stop_thresh)) > 2
        stop_displacement_threshold = multithresh(data_for_stop_thresh, 1);
    else
        stop_displacement_threshold = mean(data_for_stop_thresh);
    end
    
    % Calculate switch trial displacements for threshold
    switch_trial_displacements = [];
    if exist('switch_data', 'var') && ~isempty(switch_data)
        types = nominal(switch_data.trial_type);
        for i = 1:min(length(types), length(joystickdataSwitch_sm_velocities))
            if isempty(joystickdataSwitch_sm_velocities{1,i}) || isempty(joystickdataSwitch_sm_velocities{1,i}{1,2})
                continue; 
            end
            x_pos_data = joystickdataSwitch_sm_velocities{1,i}{1,2};
            if eq(types(i), 'switch')
                switch_trial_displacements(end+1) = max(abs(x_pos_data - x_pos_data(1)));
            end
        end
    end
    
    if isempty(switch_trial_displacements)
        warning('No switch trials found for displacement threshold calculation for patient %s', patientID);
        continue;
    end
    
    % Calculate switch displacement threshold
    data_for_switch_thresh = switch_trial_displacements(switch_trial_displacements > 0);
    if length(unique(data_for_switch_thresh)) > 2
        switch_displacement_threshold = multithresh(data_for_switch_thresh, 1);
    else
        switch_displacement_threshold = mean(data_for_switch_thresh);
    end

    % Display calculated thresholds
    fprintf('     -> Stop: Proactive=%.0f, Displacement=%.0f\n', stop_proactive_threshold, stop_displacement_threshold);
    fprintf('     -> Switch: Proactive=%.0f, Displacement=%.0f\n', switch_proactive_threshold, switch_displacement_threshold);
    fprintf('     -> Displacement Difference: %.0f (Stop-Switch)\n', stop_displacement_threshold - switch_displacement_threshold);

    % ===================================================================
    % STOP TASK TRIAL CLASSIFICATION
    % ===================================================================
    fprintf('  4. Classifying Stop Task trials with onset detection...\n');
    
    if exist('stop_data', 'var') && ~isempty(stop_data)
        % Initialize classification arrays
        stop_classification = zeros(size(stop_data, 1), 1);
        stop_ssrt_values = NaN(size(stop_data, 1), 1);
        
        types = nominal(stop_data.trial_type);
        go_indices = find(eq(types, 'go'));
        stop_indices = find(eq(types, 'stop'));
        
        % Classify Go trials (using STOP task displacement threshold)
        for i = 1:length(go_indices)
            idx = go_indices(i);
            if idx > length(stop_kinematic_data) || isempty(stop_kinematic_data(idx).x_pos)
                continue; 
            end
            
            peak_disp = max(abs(stop_kinematic_data(idx).x_pos - stop_kinematic_data(idx).x_pos(1)));
            rt = stop_data.rt(idx);
            
            % Classification: Correct if displacement > threshold AND RT < 2.5s
            if (peak_disp >= stop_displacement_threshold) && (rt < GO_TRIAL_RT_THRESHOLD)
                stop_classification(idx) = 1; % Correct Go
            else
                stop_classification(idx) = 0; % Incorrect Go
            end
        end
        
        % Classify Stop trials with onset detection
        for i = 1:length(stop_indices)
            idx = stop_indices(i);
            if idx > length(stop_kinematic_data) || isempty(stop_kinematic_data(idx).x_pos)
                continue; 
            end
            
            time_data = stop_kinematic_data(idx).time;
            x_pos_data = stop_kinematic_data(idx).x_pos;
            stop_signal_time = double(stop_data.ssd(idx));
            
            peak_disp = max(abs(x_pos_data - x_pos_data(1)));
            
            % Initial classification based on peak displacement
            if peak_disp < stop_proactive_threshold
                initial_class = 'Proactive';
            elseif peak_disp < stop_displacement_threshold
                initial_class = 'Reactive';
            else
                initial_class = 'Failed';
            end
            
            final_class = initial_class;
            ssrt_value = NaN;
            
            % For Reactive/Failed trials, detect stop onset using enhanced algorithm
            if strcmp(initial_class, 'Reactive') || strcmp(initial_class, 'Failed')
                [stop_onset_time, ssrt_value, ~, ~] = findStopPoint_enhanced(time_data, x_pos_data, stop_signal_time, initial_class);
                
                % Re-classify Reactive trials without detectable SSRT as Proactive
                if strcmp(initial_class, 'Reactive') && isnan(ssrt_value)
                    final_class = 'Proactive';
                    fprintf('    -> Trial %d re-assigned to Proactive (no SSRT detected)\n', idx);
                end
            end
            
            stop_ssrt_values(idx) = ssrt_value;
            
            % Final classification assignment
            if strcmp(final_class, 'Proactive')
                stop_classification(idx) = 2; % Proactive Correct
            elseif strcmp(final_class, 'Reactive')
                stop_classification(idx) = 1; % Reactive Correct
            else % Failed
                stop_classification(idx) = 0; % Incorrect
            end
        end
        
        % Update patient counts
        patient_counts.StopTask_Go.Correct = sum(stop_classification(go_indices) == 1);
        patient_counts.StopTask_Go.Incorrect = sum(stop_classification(go_indices) == 0);
        patient_counts.StopTask_Stop.Correct_Proactive = sum(stop_classification(stop_indices) == 2);
        patient_counts.StopTask_Stop.Correct_Reactive = sum(stop_classification(stop_indices) == 1);
        patient_counts.StopTask_Stop.Incorrect = sum(stop_classification(stop_indices) == 0);
        
        % Save classification results
        save(fullfile(patientDir, 'processed', [patientID '_stop_classification_complete.mat']), ...
            'stop_classification', 'stop_ssrt_values');
    end

    % ===================================================================
    % SWITCH TASK TRIAL CLASSIFICATION
    % ===================================================================
    fprintf('  5. Classifying Switch Task trials with onset detection...\n');
    
    if exist('switch_data', 'var') && ~isempty(switch_data)
        % Initialize classification arrays
        switch_classification = zeros(size(switch_data, 1), 1);
        switch_rt_values = NaN(size(switch_data, 1), 1);
        
        types = nominal(switch_data.trial_type);
        go_indices = find(eq(types, 'go'));
        switch_indices = find(eq(types, 'switch'));
        
        % Classify Go trials
        for i = 1:length(go_indices)
            idx = go_indices(i);
            if idx > length(switch_kinematic_data) || isempty(switch_kinematic_data(idx).x_pos)
                continue; 
            end
            
            peak_disp = max(abs(switch_kinematic_data(idx).x_pos - switch_kinematic_data(idx).x_pos(1)));
            rt = switch_data.rt(idx);
            
            % Classification: Correct if displacement > threshold AND RT < 2.5s
            if (peak_disp >= switch_displacement_threshold) && (rt < GO_TRIAL_RT_THRESHOLD)
                switch_classification(idx) = 1; % Correct Go
            else
                switch_classification(idx) = 0; % Incorrect Go
            end
        end
        
        % Classify Switch trials with onset detection
        for i = 1:length(switch_indices)
            idx = switch_indices(i);
            if idx > length(switch_kinematic_data) || isempty(switch_kinematic_data(idx).x_pos)
                continue; 
            end
            
            time_data = switch_kinematic_data(idx).time;
            x_pos_data = switch_kinematic_data(idx).x_pos;
            y_pos_data = switch_kinematic_data(idx).y_pos;
            switch_signal_time = double(switch_data.swsd(idx));
            trial_direction = char(switch_data.direction(idx));
            
            peak_disp = max(abs(x_pos_data - x_pos_data(1)));
            
            % Initial classification based on peak displacement
            if peak_disp < switch_proactive_threshold
                initial_class = 'Proactive';
            elseif peak_disp < switch_displacement_threshold
                initial_class = 'Reactive';
            else
                initial_class = 'Failed';
            end
            
            final_class = initial_class;
            rt_value = NaN;
            
            % For Reactive/Failed trials, detect switch onset using hybrid algorithm
            if strcmp(initial_class, 'Reactive') || strcmp(initial_class, 'Failed')
                [switch_onset_time, ~, ~] = findSwitchPoint_hybrid(time_data, x_pos_data, y_pos_data, switch_signal_time, trial_direction, Y_POSITION_LIMIT);
                
                if ~isnan(switch_onset_time)
                    rt_value = switch_onset_time - switch_signal_time;
                end
                
                % Apply additional constraints for Reactive trials
                if strcmp(initial_class, 'Reactive')
                    if ~isnan(switch_onset_time)
                        % Check X-axis constraint (movement must exceed proactive threshold)
                        x_at_onset = interp1(time_data, x_pos_data, switch_onset_time, 'linear', 'extrap');
                        x_constraint_failed = abs(x_at_onset - x_pos_data(1)) < switch_proactive_threshold;
                        
                        % Check Y-axis constraint (vertical movement within limits)
                        y_at_onset = interp1(time_data, y_pos_data, switch_onset_time, 'linear', 'extrap');
                        y_constraint_failed = abs(y_at_onset) > Y_POSITION_LIMIT;
                        
                        if x_constraint_failed || y_constraint_failed
                            final_class = 'Proactive';
                            fprintf('    -> Trial %d re-assigned to Proactive (constraint failed)\n', idx);
                        end
                    else
                        final_class = 'Proactive';
                        fprintf('    -> Trial %d re-assigned to Proactive (no onset detected)\n', idx);
                    end
                end
            end
            
            switch_rt_values(idx) = rt_value;
            
            % Final classification assignment
            if strcmp(final_class, 'Proactive')
                switch_classification(idx) = 2; % Proactive Correct
            elseif strcmp(final_class, 'Reactive')
                switch_classification(idx) = 1; % Reactive Correct
            else % Failed
                switch_classification(idx) = 0; % Incorrect
            end
        end
        
        % Update patient counts
        patient_counts.SwitchTask_Go.Correct = sum(switch_classification(go_indices) == 1);
        patient_counts.SwitchTask_Go.Incorrect = sum(switch_classification(go_indices) == 0);
        patient_counts.SwitchTask_Switch.Correct_Proactive = sum(switch_classification(switch_indices) == 2);
        patient_counts.SwitchTask_Switch.Correct_Reactive = sum(switch_classification(switch_indices) == 1);
        patient_counts.SwitchTask_Switch.Incorrect = sum(switch_classification(switch_indices) == 0);
        
        % Save classification results
        save(fullfile(patientDir, 'processed', [patientID '_switch_classification_complete.mat']), ...
            'switch_classification', 'switch_rt_values');
    end
    
    % Display patient-specific results
    fprintf('  -> Classification Results for %s:\n', patientID);
    fprintf('     - Stop Task (Go):    Correct: %d, Incorrect: %d\n', patient_counts.StopTask_Go.Correct, patient_counts.StopTask_Go.Incorrect);
    fprintf('     - Stop Task (Stop):  Proactive: %d, Reactive: %d, Failed: %d\n', patient_counts.StopTask_Stop.Correct_Proactive, patient_counts.StopTask_Stop.Correct_Reactive, patient_counts.StopTask_Stop.Incorrect);
    fprintf('     - Switch Task (Go):  Correct: %d, Incorrect: %d\n', patient_counts.SwitchTask_Go.Correct, patient_counts.SwitchTask_Go.Incorrect);
    fprintf('     - Switch Task (Switch): Proactive: %d, Reactive: %d, Failed: %d\n', patient_counts.SwitchTask_Switch.Correct_Proactive, patient_counts.SwitchTask_Switch.Correct_Reactive, patient_counts.SwitchTask_Switch.Incorrect);

    % Aggregate counts to global totals
    for i = 1:length(fields)
        f_names = fieldnames(patient_counts.(fields{i}));
        for j = 1:length(f_names)
            if isfield(global_counts.(fields{i}), f_names{j})
                global_counts.(fields{i}).(f_names{j}) = global_counts.(fields{i}).(f_names{j}) + patient_counts.(fields{i}).(f_names{j});
            end
        end
    end
    
    % Generate individual patient summary plot
    fprintf('  6. Generating summary plot for this patient...\n');
    generateSummaryPlot(patient_counts, sprintf('Complete Classification for Patient %s', patientID));
end

% ===================================================================
% FINAL AGGREGATE REPORTING
% ===================================================================
fprintf('\n\n===== FINAL AGGREGATE REPORTS (ALL PATIENTS) =====\n');
fprintf('\n--- Total Classification Counts (Complete Method) ---\n');
fprintf('Stop Task (Go):       Correct: %d, Incorrect: %d\n', global_counts.StopTask_Go.Correct, global_counts.StopTask_Go.Incorrect);
fprintf('Stop Task (Stop):     Proactive: %d, Reactive: %d, Failed: %d\n', global_counts.StopTask_Stop.Correct_Proactive, global_counts.StopTask_Stop.Correct_Reactive, global_counts.StopTask_Stop.Incorrect);
fprintf('Switch Task (Go):     Correct: %d, Incorrect: %d\n', global_counts.SwitchTask_Go.Correct, global_counts.SwitchTask_Go.Incorrect);
fprintf('Switch Task (Switch): Proactive: %d, Reactive: %d, Failed: %d\n', global_counts.SwitchTask_Switch.Correct_Proactive, global_counts.SwitchTask_Switch.Correct_Reactive, global_counts.SwitchTask_Switch.Incorrect);

% Generate aggregate summary plot
fprintf('\n--- Aggregate Classification Summary Figure ---\n');
generateSummaryPlot(global_counts, 'Complete Classification Summary (All Patients)');
fprintf('\n\n===== Complete batch processing finished. =====\n');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [stopTime, ssrt, reward, cumulativeReward] = findStopPoint_enhanced(time, x_pos, stop_signal_time, classification_str)
    % Enhanced stop onset detection algorithm
    % Works for both successful reactive stops and failed stop attempts
    % Uses cumulative reward function to detect movement cessation
    
    stopTime = NaN; ssrt = NaN; reward = []; cumulativeReward = [];
    
    % Movement detection threshold - more sensitive for failed stops
    if strcmp(classification_str, 'Failed')
        HORIZONTAL_DISPLACEMENT_THRESHOLD = 25;  % Lower threshold for failed stops
    else
        HORIZONTAL_DISPLACEMENT_THRESHOLD = 50;  % Original threshold for successful stops
    end
    
    time = time(:); x_pos = x_pos(:);
    x_displacement = abs(x_pos - x_pos(1));
    
    % Find when significant movement begins
    gate_open_idx = find(x_displacement > HORIZONTAL_DISPLACEMENT_THRESHOLD, 1, 'first');
    if isempty(gate_open_idx), return; end
    
    % Calculate velocities
    vx = gradient(x_pos) ./ gradient(time); 
    vx(isnan(vx)) = 0;
    speed = abs(vx);
    
    % Determine onset threshold - more sensitive for failed stops
    if strcmp(classification_str, 'Failed')
        onset_threshold = 0.05 * max(speed(gate_open_idx:end));  % Lower threshold
    else
        onset_threshold = 0.10 * max(speed(gate_open_idx:end));  % Original threshold
    end
    
    % Find movement onset
    onset_idx_relative = find(speed(gate_open_idx:end) > onset_threshold, 1, 'first');
    if isempty(onset_idx_relative)
        onset_idx = gate_open_idx; 
    else
        onset_idx = gate_open_idx + onset_idx_relative - 1; 
    end
    
    % Determine initial movement direction
    dir_window_end = min(length(vx), onset_idx + round(0.150 / mean(diff(time))));
    initial_vx_sign = sign(mean(vx(onset_idx:dir_window_end)));
    if initial_vx_sign == 0, initial_vx_sign = 1; end
    
    % Create reward signal (favors movement in initial direction)
    reward = vx * initial_vx_sign;
    cumulativeReward = cumsum(reward);
    
    % Find the stopping point (peak of cumulative reward)
    if strcmp(classification_str, 'Failed')
        % For failed stops, look for peaks after stop signal
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
        % For successful reactive stops, use global maximum
        [~, peakIdx] = max(cumulativeReward);
    end
    
    detected_stop_time = time(peakIdx);
    
    % Validate that stop attempt occurred after stop signal
    if detected_stop_time < stop_signal_time, return; end 
    
    calculated_ssrt = detected_stop_time - stop_signal_time;
    
    % Set constraints based on trial type
    MIN_RT = 0.100; 
    if strcmp(classification_str, 'Failed')
        MAX_RT = 1.500;  % More lenient for failed stops
    else
        MAX_RT = 1.000;  % Original constraint for successful stops
    end
    
    % Validate reaction time is within reasonable bounds
    if calculated_ssrt >= MIN_RT && calculated_ssrt <= MAX_RT
        stopTime = detected_stop_time;
        ssrt = calculated_ssrt;
    end
end

function [switchTime, reward, cumulativeReward] = findSwitchPoint_hybrid(time, x_pos, y_pos, switch_signal_time, trial_direction, y_limit)
    % Hybrid switch onset detection algorithm
    % Uses weighted reward function considering both X and Y movement
    % Finds switch onset by detecting peaks in cumulative reward after switch signal
    
    reward = []; cumulativeReward = []; switchTime = NaN;
    time = time(:); x_pos = x_pos(:); y_pos = y_pos(:);
    
    % Calculate velocities
    vx = gradient(x_pos) ./ gradient(time); 
    vy = gradient(y_pos) ./ gradient(time);
    vx(isnan(vx)) = 0; vy(isnan(vy)) = 0;
    
    % Determine expected direction
    if contains(trial_direction, 'left')
        initial_vx_sign = -1;
    elseif contains(trial_direction, 'right')
        initial_vx_sign = 1;
    else
        warning('Invalid horizontal direction: %s', trial_direction); 
        return; 
    end
    
    % Weighted reward function (favors horizontal movement, penalizes vertical)
    weight_x = 3.0; 
    weight_y = 1.0;
    reward = (weight_x * (vx * initial_vx_sign)) - (weight_y * abs(vy));
    
    cumulativeReward = cumsum(reward);
    time_for_reward = time;
    
    % Find switch signal time index
    signal_idx = find(time_for_reward >= switch_signal_time, 1, 'first');
    if isempty(signal_idx), return; end
    
    % Find peaks in cumulative reward after switch signal
    [peak_vals, peak_indices_relative] = findpeaks(cumulativeReward(signal_idx:end));
    if isempty(peak_vals), return; end
    
    % Sort peaks by magnitude (highest first)
    [~, sort_order] = sort(peak_vals, 'descend');
    sorted_peak_indices = signal_idx + peak_indices_relative(sort_order) - 1;
    
    % Apply constraints to find valid switch onset
    MIN_REACTION_TIME = 0.100;
    MAX_REACTION_TIME = 1.500;
    
    for i = 1:length(sorted_peak_indices)
        candidate_idx = sorted_peak_indices(i);
        candidate_time = time_for_reward(candidate_idx);
        
        % Check temporal validity
        reaction_time = candidate_time - switch_signal_time;
        is_temporally_valid = reaction_time >= MIN_REACTION_TIME && reaction_time <= MAX_REACTION_TIME;
        
        % Check spatial validity (relaxed for failed trials)
        y_at_candidate = interp1(time, y_pos, candidate_time, 'linear', 'extrap');
        is_spatially_valid = abs(y_at_candidate) <= y_limit * 2;
        
        if is_temporally_valid && is_spatially_valid
            switchTime = candidate_time; 
            return; 
        end
    end
end

function generateSummaryPlot(counts, figure_title)
    % Generate stacked bar chart summarizing trial classifications
    % Shows breakdown of Correct/Incorrect trials by task type and condition
    
    fields = {'StopTask_Go', 'StopTask_Stop', 'SwitchTask_Go', 'SwitchTask_Switch'};
    figure('Position', [100, 100, 1200, 700]);

    % Prepare data matrix for stacked bar chart
    % Columns: 1=Correct(Go), 2=Incorrect, 3=Correct(Reactive), 4=Correct(Proactive)
    bar_data = zeros(length(fields), 4);
    for i = 1:length(fields)
        fname = fields{i};
        if contains(fname, '_Go')
            bar_data(i, 1) = counts.(fname).Correct;
            bar_data(i, 2) = counts.(fname).Incorrect;
        else % Stop or Switch trials
            bar_data(i, 2) = counts.(fname).Incorrect;
            bar_data(i, 3) = counts.(fname).Correct_Reactive;
            bar_data(i, 4) = counts.(fname).Correct_Proactive;
        end
    end

    % Create labels and plot
    bar_labels = strrep(fields, '_', ' ');

    b = bar(bar_data, 'stacked');
    b(1).DisplayName = 'Correct (Go)'; b(1).FaceColor = '#2b83ba';
    b(2).DisplayName = 'Incorrect / Failed'; b(2).FaceColor = '#d7191c';
    b(3).DisplayName = 'Correct (Reactive)'; b(3).FaceColor = '#abdda4';
    b(4).DisplayName = 'Correct (Proactive)'; b(4).FaceColor = '#2ca25f';

    % Format plot
    set(gca, 'XTickLabel', bar_labels);
    xtickangle(15);
    ylabel('Number of Trials');
    title(figure_title, 'FontSize', 16, 'FontWeight', 'bold');
    legend('Location', 'northeast');
    grid on;
    box off;
end