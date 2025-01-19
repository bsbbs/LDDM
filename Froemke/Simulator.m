%% simulation begin ... excitatory activities with Background white noise (OU process) input to the excitatory cells
if ~exist(SavedResult, 'file')
    R = zeros(Ntwk.Exct.N,numel(time));
    G = zeros(Ntwk.Inhbt.N,numel(time));
    noise = zeros(Ntwk.Exct.N, 1);
    Ntwk.wEE = Ntwk.wEE_initial;
    Ntwk.wEI = Ntwk.wEI_initial;
    Ntwk.wIE = Ntwk.wIE_initial;
    IdxE = Ntwk.Visualization.IdxE;
    IdxI = Ntwk.Visualization.IdxI;
    EEdynamc = Ntwk.wEE(IdxE,IdxE);
    EIdynamc = Ntwk.wEI(IdxI,IdxE);
    IEdynamc = Ntwk.wIE(IdxE,IdxI);
    a_pre_EE = zeros(1, Ntwk.Exct.N);
    a_post_EE = zeros(Ntwk.Exct.N, 1);
    a_pre_EI = zeros(1, Ntwk.Exct.N);
    a_post_EI = zeros(Ntwk.Inhbt.N, 1);
    a_pre_IE = zeros(1, Ntwk.Inhbt.N);
    a_post_IE = zeros(Ntwk.Exct.N, 1);
    for ti = 1:(length(time)-1)
        if (mod(ti*dt, 1)==0)
            fprintf('%is.',ti*dt*1);
        end
        if (mod(ti*dt, 60)==0)
            fprintf('\n');
        end
        % update R and G activities
        dR = (-R(:,ti) + (Ntwk.wEE*R(:,ti) + Ntwk.Exct.tuning*(Seq(ti,:).*value)' + noise)./(1 + Ntwk.wIE*G(:,ti)))*dt/Ntwk.tauR;
        dG = (-G(:,ti) + Ntwk.wEI*R(:,ti))*dt/Ntwk.tauG;
        R(:,ti+1) = R(:,ti) + dR;
        G(:,ti+1) = G(:,ti) + dG;
        R(R(:, ti+1) < 0, ti+1) = 0;
        G(G(:, ti+1) < 0, ti+1) = 0;
        % update OU noise
        noise = noise + (-noise + randn(size(noise))*sqrt(dt)*Ntwk.Noise.sgm)*dt/Ntwk.Noise.tauN;
        % apply E to E STDP kernel
        preAP = R(:,ti)';
        postAP = R(:,ti);
        a_pre_EE = a_pre_EE + (-a_pre_EE/Ntwk.tau_prepost)*dt + Ntwk.A*100*preAP*dt; % the intermiediate variables for pre
        a_post_EE = a_post_EE + (-a_post_EE/Ntwk.tau_postpre)*dt + Ntwk.A*100*postAP*dt; % and post action potential induced plasticity
        a = postAP*dt*a_pre_EE - a_post_EE*preAP*dt; % the AP induced overall plasticity change
        ds = Ntwk.gamma*(1 - Ntwk.wEE).*a.*Ntwk.Cnnct_EE; % a already contains the dt information, only applied to physical connections
        % (-Ntwk.wEE/Ntwk.taus)*dt +
        Ntwk.wEE = Ntwk.wEE + ds;
        Ntwk.wEE(Ntwk.wEE<0) = 0;
        EEdynamc(ti+1) = Ntwk.wEE(IdxE,IdxE);
        % apply E to I STDP kernel
        preAP = R(:,ti)';
        postAP = G(:,ti);
        a_pre_EI = a_pre_EI + (-a_pre_EI/Ntwk.tau_prepost)*dt + Ntwk.A*100*preAP*dt; % the intermiediate variables for pre
        a_post_EI = a_post_EI + (-a_post_EI/Ntwk.tau_postpre)*dt + Ntwk.A*100*postAP*dt; % and post action potential induced plasticity
        a = postAP*dt*a_pre_EI + iSTDPsign*a_post_EI*preAP*dt; % the AP induced overall plasticity change
        ds = Ntwk.gamma*(1 - Ntwk.wEI).*a.*Ntwk.Cnnct_EI; % a already contains the dt information, only applied to physical connections
        % (-Ntwk.wEI/Ntwk.tauhet_e)*dt +
        Ntwk.wEI = Ntwk.wEI + ds;
        Ntwk.wEI(Ntwk.wEI<0) = 0;
        EIdynamc(ti+1) = Ntwk.wEI(IdxI,IdxE);
        % apply I to E STDP kernel to the course  of firing rates
        preAP = G(:,ti)';
        postAP = R(:,ti);
        a_pre_IE = a_pre_IE + (-a_pre_IE/Ntwk.tau_prepost)*dt + Ntwk.A*preAP*dt; % the intermiediate variables for pre
        a_post_IE = a_post_IE + (-a_post_IE/Ntwk.tau_postpre)*dt + Ntwk.A*postAP*dt; % and post action potential induced plasticity
        a = postAP*dt*a_pre_IE - a_post_IE*preAP*dt; % the AP induced overall plasticity change
        ds = Ntwk.gamma*(1 - Ntwk.wIE).*a.*Ntwk.Cnnct_IE; % a already contains the dt information, only applied to physical connections
        % (-Ntwk.wIE/Ntwk.taus)*dt +
        Ntwk.wIE = Ntwk.wIE + ds;
        Ntwk.wIE(Ntwk.wIE<0) = 0;
        IEdynamc(ti+1) = Ntwk.wIE(IdxE,IdxI);
    end
    fprintf('\n');
    segment = round(smpl_time/dt);
    Rsample = R(IdxE, :);
    Gsample = G(IdxI, :);
    Rsection = R(:, segment);
    Gsection = G(:, segment);
    save(SavedResult, 'Ntwk', 'EEdynamc', 'EIdynamc', 'IEdynamc', ...
        'Rsample', 'Gsample', 'Rsection', 'Gsection', 'Seq', 'dt', 'time', 'smpl_time', 'ti');
else
    load(SavedResult);
end