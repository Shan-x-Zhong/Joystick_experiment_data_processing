function reaction_time = kinematicOnsetMJT_v2(time_vector, position_vector, x_velocity_vector, displacement_threshold, rt_velocity_fraction)
% Calculates Reaction Time using an MJT method based on DISPLACEMENT.
% This version is adapted for data where the starting baseline may vary (e.g., PD tremor).
%
%   ALGORITHM LOGIC:
%   1.  Identify Target Crossing: Calculates displacement from the trial's starting
%       position. Finds the first time this DISPLACEMENT crosses a specified threshold.
%   2.  Determine Movement Bounds: (Same as before) Finds velocity zero-crossings
%       around the threshold crossing to define the start (ts) and end (te) of the movement.
%   3.  Fit MJT: (Same as before) Fits an idealized Minimum Jerk Trajectory to the segment.
%   4.  Calculate RT: (Same as before) RT is when the fitted MJT's velocity
%       first reaches a fraction of its peak.

reaction_time = NaN;

% --- Input Validation ---
if isempty(time_vector) || isempty(position_vector) || isempty(x_velocity_vector) || length(time_vector) < 10
    return;
end
if ~(length(time_vector) == length(position_vector) && length(time_vector) == length(x_velocity_vector))
    return;
end

% --- 1. Find Target Crossing based on Displacement ---
post_stim_indices = find(time_vector >= 0);
if length(post_stim_indices) < 10
    return;
end

time_ps = time_vector(post_stim_indices);
x_pos_ps = position_vector(post_stim_indices);
x_vel_ps = x_velocity_vector(post_stim_indices);

displacement_from_start = abs(x_pos_ps - x_pos_ps(1));
idx_cross_threshold = find(displacement_from_start >= displacement_threshold, 1, 'first');

if isempty(idx_cross_threshold)
    return; % Did not cross the displacement threshold
end

% --- 2. Determine Movement Bounds (ts, te) using Velocity Zero-Crossings ---
% This section remains the same as it operates relative to the threshold crossing event.
indices_before_crossing = 1:idx_cross_threshold;
vel_before_crossing = x_vel_ps(indices_before_crossing);
zero_cross_before_indices = find(vel_before_crossing(1:end-1) .* vel_before_crossing(2:end) <= 0);
if isempty(zero_cross_before_indices)
    t_s_auto = time_ps(1); 
else
    t_s_auto = time_ps(indices_before_crossing(zero_cross_before_indices(end)));
end

indices_after_crossing = idx_cross_threshold:length(time_ps);
vel_after_crossing = x_vel_ps(indices_after_crossing);
zero_cross_after_indices_rel = find(vel_after_crossing(1:end-1) .* vel_after_crossing(2:end) <= 0);
if isempty(zero_cross_after_indices_rel)
    [~, idx_P_peak_overall_ps] = max(displacement_from_start);
    t_e_auto = time_ps(idx_P_peak_overall_ps);
else
    t_e_auto = time_ps(indices_after_crossing(zero_cross_after_indices_rel(1)));
end

% --- 3. Fit MJT to the [ts, te] segment ---
if t_e_auto <= t_s_auto
    return;
end

idx_ts = find(time_ps >= t_s_auto, 1, 'first');
idx_te = find(time_ps >= t_e_auto, 1, 'first');
if isempty(idx_ts) || isempty(idx_te), return; end

P0_auto = x_pos_ps(idx_ts);
Pf_auto = x_pos_ps(idx_te);
Tm_auto = t_e_auto - t_s_auto;
if Tm_auto <= 0, return; end

% --- 4. Calculate RT from the Fitted MJT ---
num_fine_points = 500;
mjt_time_vector_fine = linspace(0, Tm_auto, num_fine_points)';

% Use a helper sub-function or inline calculation for MJT velocity
tau = mjt_time_vector_fine ./ Tm_auto;
V_mjt = ((Pf_auto - P0_auto) / Tm_auto) * (30*(tau.^2) - 60*(tau.^3) + 30*(tau.^4));
if all(V_mjt == 0), return; end

[V_mjt_peak_val, ~] = max(abs(V_mjt));
mjt_V_thresh = rt_velocity_fraction * V_mjt_peak_val;

idx_mjt_rt_init_rel = find(abs(V_mjt) >= mjt_V_thresh, 1, 'first');

if ~isempty(idx_mjt_rt_init_rel)
    t_prime_mjt_initiation = mjt_time_vector_fine(idx_mjt_rt_init_rel);
    reaction_time = t_s_auto + t_prime_mjt_initiation; 
end

end