# Shared helpers for scripts/plot_spindown.jl and scripts/plot_spindown_overlay.jl.

"""
    log_axis_resample(times, series...; logarithmic, time_min, time_max, points=1500)

Extra plotting points for a log time axis. The integrator's accepted states are sparse in the early exponential
growth (steps of ~5e7 M₀ from t ~ 1e7), and straight segments between them show as corners on a log axis. This returns
a geometric time grid on [time_min, time_max] merged with the accepted times, and each series interpolated onto it:
log-linear in t where `logarithmic[k]` is true (occupations and amplitudes, exact for exponential growth between two
states), linear in t otherwise (spin, mass, threshold). Display only: the data are the accepted states.
"""
function log_axis_resample(times, series...; logarithmic, time_min, time_max, points=1500)
    inside = filter(t -> time_min ≤ t ≤ time_max, times)
    new_times = sort(unique(vcat(exp.(range(log(max(time_min, times[1] > 0 ? times[1] : times[2])), log(min(time_max, times[end]));
                                            length=points)), inside)))
    resampled = map(enumerate(series)) do (series_index, values)
        map(new_times) do t
            upper = clamp(searchsortedfirst(times, t), 2, length(times))
            lower = upper - 1
            fraction = times[upper] == times[lower] ? 0.0 : (t - times[lower]) / (times[upper] - times[lower])
            if logarithmic[series_index] && values[lower] > 0 && values[upper] > 0
                exp((1 - fraction) * log(values[lower]) + fraction * log(values[upper]))
            else
                (1 - fraction) * values[lower] + fraction * values[upper]
            end
        end
    end
    new_times, resampled...
end
