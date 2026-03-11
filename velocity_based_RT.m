function [RT_velocity, t_init, t_peak] = velocity_based_RT(time_vector, position_vector, x_velocity, displacement_threshold, rt_velocity_fraction)
% Detects movement onset RT using velocity threshold.
%
%   t_thres : first sample where |displacement| >= displacement_threshold
%   t_peak  : sample of max absolute displacement from t_thres onward
%   t_init  : last sample before t_thres where velocity had opposite sign
%             to movement direction (falls back to trial start if none found)
%
% RT detection within [t_init, t_peak]:
%   - max_abs_vel = max of |velocity| within the window
%   - threshold   = rt_velocity_fraction * max_abs_vel
%   - Search forward from t_init, return first sample >= threshold (no interpolation)
%
% INPUTS
%   time_vector            : time in seconds — assumed to be entirely post-stimulus
%   position_vector        : baseline-corrected x position (x - x(1))
%   x_velocity             : smoothed x velocity
%   displacement_threshold : patient-specific Otsu threshold (same units as position)
%   rt_velocity_fraction   : fraction of peak velocity for onset (e.g. 0.10)
%
% OUTPUTS
%   RT_velocity : reaction time in seconds (NaN if not detected)
%   t_init      : movement initiation time (NaN if not detected)
%   t_peak      : peak displacement time   (NaN if not detected)

    RT_velocity = NaN;
    t_init      = NaN;
    t_peak      = NaN;

    % --- Basic validation ---
    if isempty(time_vector) || isempty(position_vector) || isempty(x_velocity)
        return;
    end
    if length(time_vector) ~= length(position_vector) || ...
       length(time_vector) ~= length(x_velocity)
        return;
    end

    % --- t_thres: first sample where |displacement| crosses threshold ---
    idx_thres = find(abs(position_vector) >= displacement_threshold, 1, 'first');
    if isempty(idx_thres)
        return;
    end

    % --- t_peak: max absolute displacement from t_thres onward ---
    [~, idx_peak_rel] = max(abs(position_vector(idx_thres:end)));
    idx_peak = idx_thres + idx_peak_rel - 1;
    t_peak   = time_vector(idx_peak);

    if t_peak <= time_vector(idx_thres)
        return;
    end

    % --- t_init: last opposite-sign velocity sample before t_thres ---
    movement_sign = sign(mean(x_velocity(idx_thres:idx_peak)));
    if movement_sign == 0
        movement_sign = sign(x_velocity(idx_thres));
    end
    if movement_sign == 0
        movement_sign = 1;
    end

    opposite_idx = find(sign(x_velocity(1:idx_thres)) == -movement_sign, 1, 'last');
    if isempty(opposite_idx)
        idx_init = 1;
    else
        idx_init = opposite_idx;
    end

    if time_vector(idx_init) >= time_vector(idx_thres)
        idx_init = 1;
    end

    t_init = time_vector(idx_init);

    if t_peak <= t_init
        return;
    end

    % --- Velocity RT: forward search within [t_init, t_peak] ---
    win_idx       = idx_init:idx_peak;
    vel_win       = abs(x_velocity(win_idx));
    time_win      = time_vector(win_idx);

    max_abs_vel   = max(vel_win);
    if max_abs_vel == 0
        return;
    end

    vel_threshold = rt_velocity_fraction * max_abs_vel;

    cross_idx = find(vel_win >= vel_threshold, 1, 'first');
    if isempty(cross_idx)
        return;
    end

    RT_velocity = time_win(cross_idx);

end