
% % % % % % % % % % % %  OFDM System  % % % % % % % % % % % % % % % %
%                                                                   %
% Include TDL Channel Model from TR 138 901                         %
%                                                                   %    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear
close all
format short

%% Setup

% Number of Transmitted Symbols
M=1e6;
% Modulation Order and Scaling Factor (λ)
Mod=16; A=1; 
% Symbol Rate
Rs=100e6;
% Symbol Duration
Ts=1/Rs;
Rb=Rs*log2(Mod);
B=1/Ts;

% Constellation
constell=generateConstellation(Mod,A);

% Noise power
EsNo_db=30;

% Multipath Channel 

DelaySpread=100;         %[ns]
Bc=1/(DelaySpread*1e-9); % Coherence BW
Channel_Ts=Ts/1e-9;      %[ns]

channellength=4000;

% New Channel Generarion
if 1
    req=equivalentChannelGeneration(Channel_Ts,DelaySpread,'D',0);
end    

% Effective Channel

r=req(1:channellength); 

% % Channel Selection from Preloaded Dataset

% rTest=load("ChannelData\canales_DS10_Ts10ns_TDLD_Len10000.mat")
% req=rTest.r_reali;
% r=req(2,1:channellength)

% Number of Subcarriers
N=1024; 

% Cyclic Prefix Length

% Energy accumulation criterion
powr = abs(r).^2;
pr = sum(powr);
cumr = cumsum(powr)/pr;
idx = find(cumr >= 0.999,1);

if(idx>=N)
    L=N;
else
    L=idx;
end
% % Fixed CP Length
% L=512; 

% Estimation Capability
mu=idx;

displaySystemParameters(M,Mod,Rs,Ts,DelaySpread,N,L,mu,Bc,B,(B/N));

% OFDM Blocks per Pilot
bspp=0;% 0: add pilot only at the beginning of the transmission

% Channel Impulse Response Plots

figure('Position', [0, 550, 800, 400])
stem(0:length(r)-1,real(r),'Marker','.','MarkerSize',15, 'LineWidth', 2.8),hold on
stem(0:length(r)-1,imag(r),'Marker','.','MarkerSize',15, 'LineWidth', 2.8),hold on
xline(L,  'g','LineStyle','-', 'LineWidth', 1.8),grid on
xline(idx,'b','LineStyle','--' , 'LineWidth', 1.8),hold on
xline(mu, 'r','LineStyle',':', 'LineWidth', 1.8),hold on
legend({'$\Re\{r[n]\}$','$\Im\{r[n]\}$','$L$','$L_{e}$','$\mu$'},'Interpreter','latex', 'Location','best', 'FontSize', 20)
xlim([0 1000])
set(gca, 'FontSize', 16, 'TickLabelInterpreter', 'latex')
ylabel('${r}[n]$', 'Interpreter', 'latex', 'FontSize', 20)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 20)
title('Channel Impulse Response $r[n]$', 'Interpreter', 'latex','FontSize', 14)
grid minor
hold off

%% Symbol generation
s=generateSymbols(constell,M);
Es = mean(abs(constell).^2);
% Pilot Sequence
sp=(1+1j)*[1 zeros(1,mu-1)];

%% Transmitter
% Frame Length Adjustment: if M/N is not an integer, zero-pad the input
% sequence so that its length becomes an integer multiple of N.
if mod(M,N)~=0
    ext=N-mod(M,N);
    sk=[s zeros(1,ext)];
else
    sk=s;
end
% Pilot sequence length adjustment: if the pilot sequence length is
% smaller than N, apply zero padding so that it occupies exactly one block.
if mu<N
    sptx=[sp zeros(1,N-mu)];
    pilBlks=1;
elseif mu>=N
    if mod(mu,N)~=0
        sptx=[sp zeros(1,N-mod(mu,N))];
    else
        sptx=sp;
    end
    pilBlks=length(sptx)/N;
end
numBlk = length(sk)/N;
if bspp>numBlk||bspp==0
    bspp=numBlk;
end
if mod(numBlk,bspp)==0 
    voidBlk = 0;
else
    voidBlk = bspp-mod(numBlk,bspp);
end
voidBlk;
numIns = ceil(numBlk/bspp);
% Normalizacion de la potencia de la secuencia piloto
P_sptx = mean(abs(sptx).^2);
sptx = sptx * sqrt(Es / P_sptx)*(1/sqrt(N));
Esec=sum(abs(sptx*sqrt(N)).^2);
% OFDM Modulator

SPconv = reshape(sk,N,[]);
Sm=ifft(SPconv,N);
% [sptx, sn(1:K),sptx ,sn(K+1:2K), ... ,relleno]
K = N*bspp;
% Zero padding
Sm_aux = [Sm zeros(N,voidBlk)];
Sm_mat = reshape(Sm_aux, K, numIns);
Stx = [repmat(sptx(:),1,numIns); Sm_mat];
Stx = reshape(Stx,1,[])*sqrt(N);
% Remove padding before transmission
if voidBlk > 0
    Stx = Stx(1:end-N*voidBlk);
end
Sm_p=reshape(Stx,N,[]);
% Cyclic Prefix Addition
i=1:size(Sm_p,2);
CP=[Sm_p((N-(L-1):N),i);Sm_p(:,i)];
PSconv = reshape(CP,1,[]);
xn=PSconv;

if M<1000
    PlotData(Sm,L,N,sptx,pilBlks,numIns,bspp,xn)
end
%% Channel
SNR=10^(EsNo_db/10);
Px=mean(abs(xn).^2);  
No=Px/10^(EsNo_db/10);
sigma2=No/2;

xr=conv(xn,r);
z = sqrt(sigma2) * (randn(1,length(xr)) + 1i*randn(1,length(xr)));
y=xr+z;


%% Receiver
u=y;
u(length(xn)+1:length(y))=[]; % Removal of the Convolution Extension
un=u;

% OFDM Demodulator

Ok=zeros(1,length(sk));
Um=zeros(1,length(sk));

num_rx_blocks = length(u)/(N+L);

% Pilot Sequence and OFDM Block Detection

pos_in_period = mod((1:num_rx_blocks) - 1, pilBlks + bspp);
isPilot = pos_in_period < pilBlks;
blk=find(~isPilot);
pil=find(isPilot);

um_remCP=zeros(1,num_rx_blocks*N);
rEst=zeros(length(pil),mu);

% Cyclic Prefix Removal 

for m = 1:num_rx_blocks
    % Block Processing 
    um=un(((m-1)*(N+L))+1:((m)*(N+L)));
    um_remCP(((m-1)*(N))+1:((m)*(N)))=um(L+1:end)/sqrt(N);
end

% Pilot extraction
pil_aux=zeros(length(pil),N);
for b=1:length(pil)
    pil_aux(b,:)=um_remCP((pil(b)-1)*(N)+1:pil(b)*(N));
end
pil_mat=(reshape(pil_aux.',N*pilBlks,[]).');
Esptx=sum(abs(sptx*sqrt(N)).^2);

for k = 1:length(pil)/pilBlks
    pilot=pil_mat(k,:)*sqrt(N);

    % Channel Estimation 

    rML= ML_Estimate(pilot,sptx*sqrt(N),mu,Esec);
    rEst(k,:)=rML(1:mu);
end
est=repelem(1:size(rEst,1),1,bspp);
est(:,length(blk)+1:end)=[];

% Equalization

for m=1:length(blk)

    i=blk(m);
    p=est(m);
    % FFT Block Output
    Um(((m-1)*(N))+1:((m)*(N)))=fft(um_remCP(((i-1)*(N))+1:((i)*(N))),N);
    % Equalizer Output
    Om=ZF_Equalize(N,rEst(p,:),Um(((m-1)*(N))+1:((m)*(N))));
    % P/S output
    Ok(((m-1)*(N))+1:((m)*(N)))=Om;
end
% Decition
eS=Decision(Ok(1:length(s)),constell);

re=r(1:mu);
rML=rML(1:mu);
mseML = mean((abs(rML-re)).^2)

% Symbol Error Rate
[numErr,SER]=symerr(s,eS);
fprintf('  %-10s %18d\n','Symbol Error Rate:',SER);

%% Graphics
% Constellation Diagram

figure('Position', [1220, 40, 700, 450])
% titulo=sprintf('SNR=%d',SNR(a));
h1 = scatter(real(xn), imag(xn), 'MarkerEdgeColor','r', 'Marker','.'); grid on; hold on
h2 = scatter(real(Um), imag(Um), 'MarkerEdgeColor','#e7e01f', 'Marker','.'); hold on
h3 = scatter(real(Ok(1:length(s))), imag(Ok(1:length(s))), 'MarkerEdgeColor','#00ca18', 'Marker','.'); hold on
h4 = scatter(real(constell), imag(constell), 'MarkerEdgeColor','#008fff', 'Marker','x', 'LineWidth',1.5); hold on
box on

leg = legend([h1, h2,h3,h4], ...
    {'Transmitted symbols', 'Demodulated symbols', 'Equalized symbols','Reference constellation'});
set(leg, 'Color','k', 'TextColor','w', 'Location','northeast', ...
         'FontSize',14, 'Interpreter','latex')
xlabel('In-Phase (I)',    'FontSize',14,'Color','k','Interpreter','latex')
ylabel('Quadrature (Q)', 'FontSize',14,'Color','k','Interpreter','latex')
ax = gca;
set(gcf,'Color','w')
set(ax,'Color','k')
set(ax,'GridColor','w')
set(ax,'XColor','k')
set(ax,'YColor','k')
title('Constellation Diagram', 'Interpreter', 'latex','FontSize', 14)
% title(titulo)
xlim([-11.5,11.5]) 
ylim([-12,12]) 
% Spectral analysis 
fs   = 1/Ts;                                        % [GHz]
Nfft = 2^16;
f    = (-Nfft/2:Nfft/2-1) * (fs/Nfft);             % eje en GHz, dentro de Nyquist
R = fftshift(fft(r, Nfft));
X = fftshift(fft(xn, Nfft));
Y = fftshift(fft(un, Nfft));

figure('Position', [800, 550, 1120, 400])
plot(f, 20*log10(abs(R)), 'LineWidth', 2), hold on
plot(f, 20*log10(abs(X)), 'LineWidth', 2), hold on
plot(f, 20*log10(abs(Y)), 'LineWidth', 2)
grid on
xlabel('$f$ [GHz]',        'Interpreter', 'latex', 'FontSize', 14)
ylabel('$|R(f)|$ [dB]',    'Interpreter', 'latex', 'FontSize', 14)
title('Channel Frequency Response $R(f)$', 'Interpreter', 'latex','FontSize', 14)

lgd = legend( ...
    '$R(f)$', ...
    '$X(f)$', ...
    '$Y(f)$', ...
    'Interpreter', 'latex', ...
    'FontSize',    12, ...
    'Location',    'northeast');
lgd.Color       = [0.15 0.15 0.15];   % fondo oscuro (coherente con ax.Color='k')
lgd.TextColor   = 'w';                % texto blanco
lgd.EdgeColor   = 'w';                % borde blanco

ax = gca; ax.FontSize = 12; ax.TickLabelInterpreter = 'latex'; box on
ax.Color     = 'k';
ax.GridColor = 'w';

%% Local Functions


function displaySystemParameters(N_symbols, mod_order, symbol_rate, ...
                                  symbol_duration, delay_spread, ...
                                  N_subcarriers, cp_length, est_capability,CohBW,txBW,scBW)
% =========================================================
%                   SYSTEM INPUT PARAMETERS
% =========================================================
    
        line = repmat('=', 1, 55);
        thin_line = repmat('-', 1, 55);
    
        fprintf('\n%s\n', line);
        fprintf('          %s\n', 'S Y S T E M   I N P U T   P A R A M E T E R S');
        fprintf('%s\n', line);
        fprintf('  %-30s %20s\n', 'Parameter', 'Value');
        fprintf('%s\n', thin_line);
        fprintf('  %-30s %18d\n',    'Number of Transmitted Symbols',  N_symbols);
        fprintf('  %-30s %18d\n',    'Modulation Order',               mod_order);
        fprintf('  %-30s %15.2f Mbaud/s\n', 'Symbol Rate',             symbol_rate/1e6);
        fprintf('  %-30s %15.2f ns\n',  'Symbol Duration',             symbol_duration*1e9);
        fprintf('  %-30s %15.2f Mbps\n',    'Bit Rate',                symbol_rate*log2(mod_order)/1e6);
        fprintf('  %-30s %15.2f ns\n',  'Delay Spread',                delay_spread);
        fprintf('  %-30s %18d\n',    'Number of Subcarriers',          N_subcarriers);
        fprintf('  %-30s %18d\n',    'Cyclic Prefix Length',           cp_length);
        fprintf('  %-30s %18d\n',    'Estimation Capability',          est_capability);
        fprintf('  %-30s %15.4f MHz\n','Coherence Bandwidth',      CohBW/1e6);
        fprintf('  %-30s %15.2f MHz\n','Tx Bandwidth',             txBW/1e6);
        fprintf('  %-30s %15.2f kHz\n','SubCarrier Bandwidth',             scBW/1e3);
    
        fprintf('%s\n\n', line);
    
    end
    
    function sy = generateSymbols(constelation, n)
        sy=[datasample(constelation,n, 'Replace', true)];
    end
    
    function output = generateConstellation(M, A)
    posibleConstellations=[4 16 64 256 1024 4096];
    check=sum(ismember(M,posibleConstellations));
    assert(A > 0,'Zero Symbol Power error')
    if check>0
        
        values=-(sqrt(M)-1):2:(sqrt(M)-1);
        [Re, Im] =meshgrid(values, values);
        constell=A*(Re+1i*Im);
        constell=reshape(constell,1,numel(constell));
        output=constell;
       
        
    else
        assert(check<=0)
    end

end

function output = Decision(input, constellation)
    u = input(:);          % Nx1
    c = constellation(:).';% 1xM
    
    [~,I] = min(abs(u - c),[],2);
    output = c(I);

end


function [] = PlotData(Sm,L,N,sptx,pilBlks,numIns,bspp,xn)
    % --- Gráfica stem con colores diferenciados ---
    numSym   = size(Sm, 2);
    cpLen    = L;
    blockLen = N + cpLen;
    pilotLen = length(sptx);          % longitud de la secuencia piloto (con zero-pad)
    pilotBlockLen = pilBlks * (N + cpLen);

    figure('Position', [100, 100, 1500, 450]); hold on;
    hCP    = [];
    hOFDM  = [];
    hPilot = [];

    % ── Bloques piloto (uno cada bspp+1 bloques) ──────────────────────────
    for ins = 1:numIns
        % Posición de inicio del bloque piloto dentro de xn
        % Estructura: [CP+piloto | CP+datos x bspp | CP+piloto | ...]
        blkStart = (ins - 1) * (pilotBlockLen + blockLen * bspp) + 1;

        idxPilotCP   = blkStart : blkStart + cpLen - 1;
        idxPilotSym  = blkStart + cpLen : blkStart + pilotBlockLen - 1;

        % CP del piloto — color verde oscuro
        h3cp = stem(idxPilotCP, real(xn(idxPilotCP)), ...
            'Color',           [0.10 0.60 0.30], ...
            'MarkerFaceColor', [0.10 0.60 0.30], ...
            'Marker',          '.', ...
            'MarkerSize',      12, ...
            'LineWidth',       1.4, ...
            'HandleVisibility','off');

        % Símbolo piloto — verde lima
        h3 = stem(idxPilotSym, real(xn(idxPilotSym)), ...
            'Color',           [0.18 0.80 0.44], ...
            'MarkerFaceColor', [0.18 0.80 0.44], ...
            'Marker',          '.', ...
            'MarkerSize',      12, ...
            'LineWidth',       1.4, ...
            'HandleVisibility','off');

        if ins == 1
            hPilot = h3;
            set(hPilot, 'HandleVisibility', 'on');
        end
    end

    % ── Bloques de datos ───────────────────────────────────────────────────
    for ins = 1:numIns
        pilotOffset = (ins - 1) * (pilotBlockLen + blockLen * bspp) + pilotBlockLen;

        for blk = 1:bspp
            idxStart = pilotOffset + (blk - 1) * blockLen + 1;
            if idxStart > length(xn), break; end

            idxCP   = idxStart : min(idxStart + cpLen - 1,        length(xn));
            idxOFDM = idxStart + cpLen : min(idxStart + blockLen - 1, length(xn));

            % Prefijo cíclico — azul
            h1 = stem(idxCP, real(xn(idxCP)), ...
                'Color',           [0.20 0.45 0.85], ...
                'MarkerFaceColor', [0.20 0.45 0.85], ...
                'Marker',          '.', ...
                'MarkerSize',      12, ...
                'LineWidth',       1.4, ...
                'HandleVisibility','off');

            % Símbolo OFDM — naranja
            h2 = stem(idxOFDM, real(xn(idxOFDM)), ...
                'Color',           [0.90 0.45 0.10], ...
                'MarkerFaceColor', [0.90 0.45 0.10], ...
                'Marker',          '.', ...
                'MarkerSize',      12, ...
                'LineWidth',       1.4, ...
                'HandleVisibility','off');

            if ins == 1 && blk == 1
                hCP   = h1;
                hOFDM = h2;
                set(hCP,   'HandleVisibility', 'on');
                set(hOFDM, 'HandleVisibility', 'on');
            end
        end
    end

    leg = legend([hCP, hOFDM, hPilot], ...
        {'Prefijo c\''iclico', 'S\''imbolos OFDM', 'Secuencia piloto'});
    set(leg, 'Interpreter','latex', 'FontSize',16, 'Location','northeast')
    set(gca, ...
        'FontSize',             14, ...
        'FontName',             'Times New Roman', ...
        'TickLabelInterpreter', 'latex', ...
        'LineWidth',            0.8);
    xlabel('Muestra',                'Interpreter','latex', 'FontSize',15)
    ylabel('Amplitud $\Re\{\hat{x}[n]\}$', 'Interpreter','latex', 'FontSize',18)
    xlim([1, length(xn)])
    grid on
    hold off
end
