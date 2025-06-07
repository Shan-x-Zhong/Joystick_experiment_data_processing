function reaction_time = kinematicOnsetMJT(time_vector, position_vector, x_velocity_vector, target_pos_threshold, rt_velocity_fraction)
% Calculates Reaction Time using a Minimum Jerk Trajectory (MJT) method.
%
%
%   ALGORITHM LOGIC:
%   1.  Identify Target Crossing: Finds the first time the absolute position
%       crosses a specified target threshold. This serves as a robust anchor
%       point within a confirmed, large-scale movement.
%   2.  Determine Movement Bounds: Works backward and forward from the target
%       crossing point to find the start (ts) and end (te) of the main
%       movement by identifying the nearest zero-crossings of the
%       signed x-velocity.
%   3.  Fit MJT: Directly calculates the parameters for an idealized, smooth
%       Minimum Jerk Trajectory between the determined start and end bounds
%       (ts, te).
%   4.  Calculate RT: The final reaction time is determined as the moment
%       the fitted MJT's velocity first reaches a specified fraction (e.g., 10%)
%       of its own peak velocity.
%
%   INPUTS:
%   - time_vector:        (N x 1) vector of time points for the trial.
%                         Assumed to start at or before stimulus onset (t=0).
%   - position_vector:    (N x 1) vector of x-position data.
%   - x_velocity_vector:  (N x 1) vector of signed x-velocity data.
%   - target_pos_threshold: Scalar value for the position threshold (e.g., 25000).
%   - rt_velocity_fraction: Scalar fraction for RT threshold (e.g., 0.10 for 10%).
%
%   OUTPUT:
%   - reaction_time:      Scalar value of the calculated RT. Returns NaN if RT
%                         cannot be determined (e.g., target not crossed,
%                         movement segment invalid).

reaction_time = NaN;

if isempty(time_vector) || isempty(position_vector) || isempty(x_velocity_vector) || ...
   length(time_vector) < 100 
    fprintf('  RT_FUNC_ERROR: Input data is empty or too short.\n');
    return;
end

if ~(length(time_vector) == length(position_vector) && length(time_vector) == length(x_velocity_vector))
    fprintf('  RT_FUNC_ERROR: Input vectors must be of the same length.\n');
    return;
end

% --- 1. Find Target Crossing ---
post_stim_indices = find(time_vector >= 0);
if length(post_stim_indices) < 50
    fprintf('  RT_FUNC_NOTE: Insufficient post-stimulus data.\n');
    return;
end

time_ps = time_vector(post_stim_indices);
x_pos_ps = position_vector(post_stim_indices);
x_vel_ps = x_velocity_vector(post_stim_indices);

idx_cross_threshold = find(abs(x_pos_ps) >= target_pos_threshold, 1, 'first');

if isempty(idx_cross_threshold)
    fprintf('  RT_FUNC_NOTE: Did not cross target threshold (+/-%d).\n', target_pos_threshold);
    return;
end
T_cross_threshold = time_ps(idx_cross_threshold);
fprintf('  RT_FUNC_INFO: Crossed target threshold at t=%.3f s.\n', T_cross_threshold);

% --- 2. Determine Movement Bounds (ts, te) using Velocity Zero-Crossings ---
% Find last zero-crossing BEFORE the threshold crossing (for ts)
indices_before_crossing = 1:idx_cross_threshold;
vel_before_crossing = x_vel_ps(indices_before_crossing);
zero_cross_before_indices = find(vel_before_crossing(1:end-1) .* vel_before_crossing(2:end) <= 0);

if isempty(zero_cross_before_indices)
    t_s_auto = time_ps(1); 
    fprintf('  RT_FUNC_INFO: No zero-crossing found before T_cross; using stimulus onset as MJT start.\n');
else
    t_s_auto = time_ps(indices_before_crossing(zero_cross_before_indices(end)));
end

% Find first zero-crossing AFTER the threshold crossing (for te)
indices_after_crossing = idx_cross_threshold:length(time_ps);
vel_after_crossing = x_vel_ps(indices_after_crossing);
zero_cross_after_indices_rel = find(vel_after_crossing(1:end-1) .* vel_after_crossing(2:end) <= 0);

if isempty(zero_cross_after_indices_rel)
    [~, idx_P_peak_overall_ps] = max(abs(x_pos_ps - x_pos_ps(1)));
    t_e_auto = time_ps(idx_P_peak_overall_ps);
    fprintf('  RT_FUNC_INFO: No zero-crossing found after T_cross; using peak position time as MJT end.\n');
else
    t_e_auto = time_ps(indices_after_crossing(zero_cross_after_indices_rel(1)));
end
fprintf('  RT_FUNC_INFO: Automatic MJT segment bounds: ts=%.3f, te=%.3f\n', t_s_auto, t_e_auto);

% --- 3. Fit MJT to the [ts, te] segment ---
if t_e_auto <= t_s_auto
    fprintf('  RT_FUNC_ERROR: MJT segment invalid (te <= ts). Cannot calculate RT.\n');
    return;
end

idx_ts = find(time_ps == t_s_auto, 1);
idx_te = find(time_ps == t_e_auto, 1);
P0_auto = x_pos_ps(idx_ts);
Pf_auto = x_pos_ps(idx_te);
Tm_auto = t_e_auto - t_s_auto;

% Generate finely sampled MJT for accurate peak/threshold finding
sampling_rate = 1/mean(diff(time_vector));
if isnan(sampling_rate) || isinf(sampling_rate) || sampling_rate < 1, sampling_rate = 200; end 
num_fine_points = max(100, round(Tm_auto * sampling_rate));
mjt_time_vector_fine = linspace(0, Tm_auto, num_fine_points)';

mjt_V_auto_seg = generate_mjt_velocity(P0_auto, Pf_auto, Tm_auto, mjt_time_vector_fine);

if isempty(mjt_V_auto_seg) || all(mjt_V_auto_seg == 0)
    fprintf('  RT_FUNC_ERROR: Generated MJT velocity is empty or zero. Cannot calculate RT.\n');
    return;
end

% --- 4. Calculate RT from the Fitted MJT ---
[V_mjt_peak_val, ~] = max(abs(mjt_V_auto_seg));
mjt_V_thresh = rt_velocity_fraction * V_mjt_peak_val;

idx_mjt_rt_init_rel = find(abs(mjt_V_auto_seg) >= mjt_V_thresh, 1, 'first');

if ~isempty(idx_mjt_rt_init_rel)
    t_prime_mjt_initiation = mjt_time_vector_fine(idx_mjt_rt_init_rel);
    reaction_time = t_s_auto + t_prime_mjt_initiation; 
    fprintf('  RT_FUNC_SUCCESS: Automatic RT found: %.3f s\n', reaction_time);
else
    fprintf('  RT_FUNC_NOTE: Did not cross %.1f%% MJT peak velocity threshold. RT is NaN.\n', rt_velocity_fraction*100);
    reaction_time = NaN;
end

end

function V_mjt = generate_mjt_velocity(P0, Pf, Tm, t_relative)
    if Tm <= 1e-9 
        V_mjt = zeros(size(t_relative)); 
        return;
    end
    tau = t_relative ./ Tm;
    tau_calc = min(max(tau, 0), 1);
    V_mjt = ((Pf - P0) / Tm) * (30*(tau_calc.^2) - 60*(tau_calc.^3) + 30*(tau_calc.^4));
    V_mjt(t_relative > Tm | t_relative < 0) = 0; 
end