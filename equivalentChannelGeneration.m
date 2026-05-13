function r = equivalentChannelGeneration(Ts,DS,sel,Graphix)
%% --- Modelos TDL ---
    tdlA_delay   = [0.0000 0.3819 0.4025 0.5868 0.4610 0.5375 0.6708 0.5750 0.7618 1.5375 1.8978 2.2242 2.1718 2.4942 2.5119 3.0582 4.0810 4.4579 4.5695 4.7966 5.0066 5.3043 9.6586];
    tdlA_power_dB = [-13.4 0 -2.2 -4 -6 -8.2 -9.9 -10.5 -7.5 -15.9 -6.6 -16.7 -12.4 -15.2 -10.8 -11.3 -12.7 -16.2 -18.3 -18.9 -16.6 -19.9 -29.7];
    tdlB_delay   = [0.0000 0.1072 0.2155 0.2095 0.2870 0.2986 0.3752 0.5055 0.3681 0.3697 0.5700 0.5283 1.1021 1.2756 1.5474 1.7842 2.0169 2.8294 3.0219 3.6187 4.1067 4.2790 4.7834];
    tdlB_power_dB = [0 -2.2 -4 -3.2 -9.8 -1.2 -3.4 -5.2 -7.6 -3 -8.9 -9 -4.8 -5.7 -7.5 -1.9 -7.6 -12.2 -9.8 -11.4 -14.9 -9.2 -11.3];
    tdlC_delay   = [0 0.2099 0.2219 0.2329 0.2176 0.6366 0.6448 0.6560 0.6584 0.7935 0.8213 0.9336 1.2285 1.3083 2.1704 2.7105 4.2589 4.6003 5.4902 5.6077 6.3065 6.6374 7.0427 8.6523];
    tdlC_power_dB = [-4.4 -1.2 -3.5 -5.2 -2.5 0 -2.2 -3.9 -7.4 -7.1 -10.7 -11.1 -5.1 -6.8 -8.7 -13.2 -13.9 -13.9 -15.8 -17.1 -16 -15.7 -21.6 -22.8];
    tdlD_delay   = [0 0.035 0.612 1.363 1.405 1.804 2.596 1.775 4.042 7.937 9.424 9.708 12.525];
    tdlD_power_dB = [-0.2 -18.8 -21 -22.8 -17.9 -20.1 -21.9 -22.9 -27.8 -23.6 -24.8 -30.0 -27.7];
    tdlE_delay   =[0 0.035 0.612 1.363 1.405 1.775 1.804 2.596 4.042 7.937 9.424 9.708 12.525];
    tdlE_power_dB =[-0.2 -18.8 -21 -22.8 -17.9 -22.9 -20.1 -21.9 -27.8 -23.6 -24.8 -30 -27.7];
    switch sel
        case 'A'
        power_dB  = tdlA_power_dB; 
        normDelay = tdlA_delay;
        case 'B'
        power_dB  = tdlB_power_dB; 
        normDelay = tdlB_delay;
        case 'C'
        power_dB  = tdlC_power_dB; 
        normDelay = tdlC_delay;
        case 'D'
        power_dB  = tdlD_power_dB; 
        normDelay = tdlD_delay;
        case 'E'
        power_dB  = tdlE_power_dB; 
        normDelay = tdlE_delay;
    end
    P = 10.^(power_dB/10);
    L = length(P);
    ScalDel = normDelay.*DS; 
    % --- Definición temporal ---
    dt = Ts/100;
    % t = -40*Ts:dt:ScalDel(end)+10*Ts;
    t = -40*Ts:dt:10000*Ts;
    % --- Pulso equivalente g(t) ---
    g = sinc(t/Ts);
    % --- Coeficientes complejos del TDL ---
    if sel=='D'||sel=='E'
        % Para el caso de LOS
        K_rice = 13.3;
        K_lin = 10^(K_rice/10);
        s = sqrt(P(1) * K_lin/(K_lin+1));
        scatter_power = P(1)/(K_lin+1);
        gamma(1) = s + sqrt(scatter_power/2)*(randn+1j*randn);
        gamma(2:L) = sqrt(P(2:L)/2).*(randn(1,L-1)+1j*randn(1,L-1));
    else
        % Para el caso de NLOS
        gamma = sqrt(P/2).*(randn(1,L)+1j*randn(1,L));
    end
    % --- Construcción de r(t) ---
    rt = zeros(size(t));
    for i = 1:L
        x=gamma(i) * interp1(t, g, t-ScalDel(i), 'linear', 0);
        rt = rt + x;
    end
    % --- Muestreo de r(t) ---
    n = 0:floor(max(t)/Ts);
    tn = n*Ts;
    r = interp1(t, rt, tn, 'linear', 0);

    %% --- Gráficas ---
    if Graphix==1
    t_lim = [0, t(end)];
    % --- Subplot 1: Coeficientes gamma ---
    figure('Position', [0,200, 1000, 250]) 
    stem(ScalDel, real(gamma), 'Marker', '.', 'MarkerSize', 24), hold on
    stem(ScalDel, imag(gamma), 'Marker', '.', 'MarkerSize', 24)
    grid on
    xlim(t_lim)
    xlabel('$\mathrm{Tiempo\ [ns]}$', 'Interpreter', 'latex', 'FontSize', 20)
    lg=legend('$\Re\{\gamma_i\}$', '$\Im\{\gamma_i\}$', 'Interpreter', 'latex', 'Location', 'best');
    lg.FontSize = 20;
    ylabel('$\sum_{i=0}^{N_p-1}\limits \gamma_i \, \delta(t - \tau_i)$', 'Interpreter', 'latex')
    ax = gca; 
    ax.FontSize = 18; 
    ax.TickLabelInterpreter = 'latex';
    ax.YLabel.FontSize = 25;
    box on; hold off
    
    % --- Subplot 2: r(t) continua ---
    figure('Position', [0,450, 1000, 250])
    hold on
    plot(t, real(rt), '-',  'LineWidth', 2)
    plot(t, imag(rt), '-.', 'LineWidth', 2)
    grid on
    xlim(t_lim)
    xlabel('$\mathrm{Tiempo\ [ns]}$', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('$r(t)$', 'Interpreter', 'latex')
    lg=legend({'$\Re\{r(t)\}$','$\Im\{r(t)\}$'}, 'Interpreter', 'latex', 'Location', 'best');
    lg.FontSize = 20;
    ax = gca; 
    ax.FontSize = 18; 
    ax.TickLabelInterpreter = 'latex';
    ax.YLabel.FontSize = 30;
    box on; hold off
    
    % --- Subplot 3: r[n] muestreada ---
    figure('Position', [0,100, 1000, 250]) 
    stem(real(r), 'Marker', '.', 'MarkerSize', 20), hold on
    stem(imag(r), 'Marker', '.', 'MarkerSize', 20),hold on
    xline(idx99, 'B','LineStyle','--', 'LineWidth', 1.8)
    grid on
    xlim([1 length(r)])
    xlabel('$\mathrm{n}$', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('$r[n]$', 'Interpreter', 'latex')
    lg=legend({'$\Re\{r[n]\}$','$\Im\{r[n]\}$'}, 'Interpreter', 'latex', 'Location', 'best');
    lg.FontSize = 20;
    ax = gca; 
    ax.FontSize = 18; 
    ax.TickLabelInterpreter = 'latex';
    ax.YLabel.FontSize = 30;
    box on; hold off
    %% --- Respuesta en Frecuencia del Canal ---
    fs   = 1/Ts;                                        % [GHz]
    Nfft = 2048;
    f    = (-Nfft/2:Nfft/2-1) * (fs/Nfft);             % eje en GHz, dentro de Nyquist

    H_shift = fftshift(fft(r, Nfft));

    figure('Position', [1280, 580, 640, 420])
    plot(f, 20*log10(abs(H_shift)), 'LineWidth', 2)
    grid on
    xlabel('$f$ [GHz]',                'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$|R(f)|$ [dB]',       'Interpreter', 'latex', 'FontSize', 14)
    title('Respuesta en frecuencia de $r[n]$', 'Interpreter', 'latex')
    ax = gca; ax.FontSize = 12; ax.TickLabelInterpreter = 'latex'; box on
    H=(20*log10(abs(H_shift)));
    ylim([-20 2*max(H)])
    
    % --- Respuesta analítica ---
    f_eval      = linspace(-fs*3, fs*3, 10000);         % GHz, mismo rango que FFT
    H_analitico = zeros(size(f_eval));
    for i = 1:L
        H_analitico = H_analitico + gamma(i) .* exp(-1j*2*pi*f_eval*ScalDel(i));
    end
    % f_eval [GHz] x ScalDel [ns] = adimensional ✅

    figure('Position', [1280, 75, 640, 420])
    plot(f_eval, 20*log10(abs(H_analitico)), 'LineWidth', 1.5)
    grid on
    xlabel('$f$ [GHz]',          'Interpreter', 'latex', 'FontSize', 14)
    ylabel('$|H(f)|$ [dB]',      'Interpreter', 'latex', 'FontSize', 14)
    % title('Respuesta analítica del canal TDL', 'Interpreter', 'latex')
    ax = gca; ax.FontSize = 12; ax.TickLabelInterpreter = 'latex'; box on
    xlim([-fs*3 fs*3])
    end
end