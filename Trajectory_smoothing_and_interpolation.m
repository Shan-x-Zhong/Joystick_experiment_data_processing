clear;
rng default;
clc;
close all;

baseDir = '/Users/x/Desktop';

patientFolders = dir(fullfile(baseDir, 'P*'));
normtime = 100000; % Normalization factor for time

for p = 1:length(patientFolders)
    patientDir = fullfile(baseDir, patientFolders(p).name);
    fprintf('Processing patient: %s\n', patientFolders(p).name);

    % ----- Process STOP task files -----
    stopFiles = dir(fullfile(patientDir, '*stop1_behav_exp*.mat'));
    for f = 1:length(stopFiles)
        fileName = fullfile(patientDir, stopFiles(f).name);
        fprintf('  Processing STOP file: %s\n', stopFiles(f).name);

        load(fileName); 

        joystickdataStop_sm_velocities = cell(1, length(joystickdataStop));
        mov_onset_stop = zeros(length(joystickdataStop), 1);

        for i = 1:length(joystickdataStop)
            if isempty(joystickdataStop{1,i}) || ~iscell(joystickdataStop{1,i}) || length(joystickdataStop{1,i}) < 3
                fprintf('    Skipping invalid trial %d in STOP task\n', i);
                joystickdataStop_sm_velocities{1,i} = cell(1,6); 
                mov_onset_stop(i) = NaN;
                continue;
            end

            try
                dat = []; 
                dat(:,1) = joystickdataStop{1,i}{1,1}; % Original time data
                dat(:,2) = joystickdataStop{1,i}{1,2}; % Original x position data
                dat(:,3) = joystickdataStop{1,i}{1,3}; % Original y position data

                if isempty(dat) || size(dat,1) == 0 || size(dat,2) < 3
                    fprintf('    Skipping empty data in trial %d in STOP task\n', i);
                    joystickdataStop_sm_velocities{1,i} = cell(1,6);
                    mov_onset_stop(i) = NaN;
                    continue;
                end

                [time_n, x_pos_n, y_pos_n, total_vel_n, x_vel_n, y_vel_n] = smoothtrajectory3(dat, normtime);

            catch ME 
                fprintf('    Error processing trial %d in STOP task: %s\n', i, ME.message);
                joystickdataStop_sm_velocities{1,i} = cell(1,6);
                mov_onset_stop(i) = NaN;
                continue;
            end

            joystickdataStop_sm_velocities{1,i}{1,1} = time_n;       % Normalized Time
            joystickdataStop_sm_velocities{1,i}{1,2} = x_pos_n;     % Normalized X position
            joystickdataStop_sm_velocities{1,i}{1,3} = y_pos_n;     % Normalized Y position
            joystickdataStop_sm_velocities{1,i}{1,4} = total_vel_n; % Total velocity
            joystickdataStop_sm_velocities{1,i}{1,5} = x_vel_n;     % X velocity
            joystickdataStop_sm_velocities{1,i}{1,6} = y_vel_n;     % Y velocity
        end

        saveDir = fullfile(patientDir, 'processed');
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        save(fullfile(saveDir, [patientFolders(p).name '_stop_onset.mat']), 'mov_onset_stop');
        save(fullfile(saveDir, [patientFolders(p).name '_stop_sm.mat']), 'joystickdataStop_sm_velocities');
    end

    % ----- Process SWITCH task files -----
    switchFiles = dir(fullfile(patientDir, '*switch_behav_exp*.mat'));
    for f = 1:length(switchFiles)
        fileName = fullfile(patientDir, switchFiles(f).name);
        fprintf('  Processing SWITCH file: %s\n', switchFiles(f).name);

        load(fileName);

        joystickdataSwitch_sm_velocities = cell(1, length(joystickdataSwitch));
        mov_onset_switch = zeros(length(joystickdataSwitch), 1);

        for i = 1:length(joystickdataSwitch)
            if isempty(joystickdataSwitch{1,i}) || ~iscell(joystickdataSwitch{1,i}) || length(joystickdataSwitch{1,i}) < 3
                fprintf('    Skipping invalid trial %d in SWITCH task\n', i);
                joystickdataSwitch_sm_velocities{1,i} = cell(1,6);
                mov_onset_switch(i) = NaN;
                continue;
            end

            try
                dat = [];
                dat(:,1) = joystickdataSwitch{1,i}{1,1};
                dat(:,2) = joystickdataSwitch{1,i}{1,2};
                dat(:,3) = joystickdataSwitch{1,i}{1,3};

                if isempty(dat) || size(dat,1) == 0 || size(dat,2) < 3
                    fprintf('    Skipping empty data in trial %d in SWITCH task\n', i);
                    joystickdataSwitch_sm_velocities{1,i} = cell(1,6);
                    mov_onset_switch(i) = NaN;
                    continue;
                end

                [time_n, x_pos_n, y_pos_n, total_vel_n, x_vel_n, y_vel_n] = smoothtrajectory3(dat, normtime);
            catch ME
                fprintf('    Error processing trial %d in SWITCH task: %s\n', i, ME.message);
                joystickdataSwitch_sm_velocities{1,i} = cell(1,6);
                mov_onset_switch(i) = NaN;
                continue;
            end

            joystickdataSwitch_sm_velocities{1,i}{1,1} = time_n;
            joystickdataSwitch_sm_velocities{1,i}{1,2} = x_pos_n;
            joystickdataSwitch_sm_velocities{1,i}{1,3} = y_pos_n;
            joystickdataSwitch_sm_velocities{1,i}{1,4} = total_vel_n;
            joystickdataSwitch_sm_velocities{1,i}{1,5} = x_vel_n;
            joystickdataSwitch_sm_velocities{1,i}{1,6} = y_vel_n;

        end

        saveDir = fullfile(patientDir, 'processed');
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        save(fullfile(saveDir, [patientFolders(p).name '_switch_onset.mat']), 'mov_onset_switch');
        save(fullfile(saveDir, [patientFolders(p).name '_switch_sm.mat']), 'joystickdataSwitch_sm_velocities');
    end

    % ----- Process GO task files -----
    goFiles = dir(fullfile(patientDir, '*gotask_behav_exp*.mat'));
    for f = 1:length(goFiles)
        fileName = fullfile(patientDir, goFiles(f).name);
        fprintf('  Processing GO file: %s\n', goFiles(f).name);

        load(fileName);

        joystickdataGo_sm_velocities = cell(1, length(joystickdataGo));
        mov_onset_go = zeros(length(joystickdataGo), 1);

        for i = 1:length(joystickdataGo)
            if isempty(joystickdataGo{1,i}) || ~iscell(joystickdataGo{1,i}) || length(joystickdataGo{1,i}) < 3
                fprintf('    Skipping invalid trial %d in GO task\n', i);
                joystickdataGo_sm_velocities{1,i} = cell(1,6);
                mov_onset_go(i) = NaN;
                continue;
            end

            try
                dat = [];
                dat(:,1) = joystickdataGo{1,i}{1,1};
                dat(:,2) = joystickdataGo{1,i}{1,2};
                dat(:,3) = joystickdataGo{1,i}{1,3};

                if isempty(dat) || size(dat,1) == 0 || size(dat,2) < 3
                    fprintf('    Skipping empty data in trial %d in GO task\n', i);
                    joystickdataGo_sm_velocities{1,i} = cell(1,6);
                    mov_onset_go(i) = NaN;
                    continue;
                end

                [time_n, x_pos_n, y_pos_n, total_vel_n, x_vel_n, y_vel_n] = smoothtrajectory3(dat, normtime);
            catch ME
                fprintf('    Error processing trial %d in GO task: %s\n', i, ME.message);
                joystickdataGo_sm_velocities{1,i} = cell(1,6);
                mov_onset_go(i) = NaN;
                continue;
            end

            joystickdataGo_sm_velocities{1,i}{1,1} = time_n;
            joystickdataGo_sm_velocities{1,i}{1,2} = x_pos_n;
            joystickdataGo_sm_velocities{1,i}{1,3} = y_pos_n;
            joystickdataGo_sm_velocities{1,i}{1,4} = total_vel_n;
            joystickdataGo_sm_velocities{1,i}{1,5} = x_vel_n;
            joystickdataGo_sm_velocities{1,i}{1,6} = y_vel_n;

        end

        saveDir = fullfile(patientDir, 'processed');
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        save(fullfile(saveDir, [patientFolders(p).name '_go_onset.mat']), 'mov_onset_go');
        save(fullfile(saveDir, [patientFolders(p).name '_go_sm.mat']), 'joystickdataGo_sm_velocities');
    end
end


function [time_n, x_pos_n, y_pos_n, total_vel_n, x_vel_n, y_vel_n] = smoothtrajectory3(dat, normtime_points)
    t_original = dat(:,1);  % Actual time values (not length!)
    t_min = min(t_original);
    t_max = max(t_original);
    
    % Fit splines directly against actual time
    cp_x_pos_fit = csaps(t_original, dat(:,2), 0.999999); 
    cp_y_pos_fit = csaps(t_original, dat(:,3), 0.999999);
    
    % Derivatives w.r.t actual time
    dcp_x_pos = fnder(cp_x_pos_fit); 
    dcp_y_pos = fnder(cp_y_pos_fit); 
    
    % Uniform resampling in TRUE TIME domain
    time_n = linspace(t_min, t_max, normtime_points);
    
    % Evaluate at true time points
    x_pos_n = fnval(cp_x_pos_fit, time_n);
    y_pos_n = fnval(cp_y_pos_fit, time_n);
    x_vel_n = fnval(dcp_x_pos, time_n); 
    y_vel_n = fnval(dcp_y_pos, time_n); 
    
    total_vel_n = sqrt(x_vel_n.^2 + y_vel_n.^2);
end