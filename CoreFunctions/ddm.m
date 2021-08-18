function [Evidence, rt, choice] = ddm(Vinput, z, a, v, t0, stimdur, dur, dt, stoprule)
total_time_steps = round(dur/dt);
stim_duration = round(stimdur/dt);
dV = Vinput(2) - Vinput(1);
Evidence = z;
rt = NaN;
for ti = 2:total_time_steps
    if ti == stim_duration
        dV = 0;
    end
    Evidence(ti) = Evidence(ti-1) + dV*v*dt + randn(1)*sqrt(dt);
    inside = (Evidence(ti) >= a)*2 + (Evidence(ti) <= 0);
    if inside > 0 && isnan(rt)
        rt = ti*dt + t0;
        choice = inside;
        if stoprule == 1
            break;
        end
    end
end