function result = A_insilico_1_eisim(EI_ratio, MAKE_PLOT)
% A_insilico_1_eisim - simulate LFP based on model from Gao, Peterson, & Voytek
%
% INPUT
%   EI_ratio = [2:0.2:6];
%   MAKE_PLOT = set to 1 to make plot
%
% Usage 
% EI_ratio = [2:0.2:6];
% MAKE_PLOT = 1;
% result = A_insilico_1_eisim(EI_ratio, MAKE_PLOT);
%

%%
outpath = '/Users/mvlombardo/projects/EI_hurst/reproAnalysis/ephys_sim';

%% set up simulation
dt = 0.001; %simulation time step
fs = 1/dt; %sampling rate
tk = 0:dt:1; %PSC kernel time vector
t = 0:dt:60*2; %simulation time vector

%spike train parameters
FR_E = 2;
FR_I = 5;
N_E = 8000;
N_I = 2000;

%ampa/gaba PSC kernels
Vr = -65;
Ee = 0;
Ei = -80;
AMPA_tau = [0.1 2]/1000;
GABA_tau = [0.5 10]/1000;
kA = syn_kernel(tk,AMPA_tau);
kG = syn_kernel(tk,GABA_tau);

%% simulate with different EI ratios
% simulate 128 regions
num_trs = 128;  
for i = 1:length(EI_ratio)
    
    leg_labels{i} = sprintf('1:%0.1f',EI_ratio(i));
    boost = EI_ratio(i)./((N_I*FR_I*sum(kG))/(N_E*FR_E*sum(kA)));
    disp(EI_ratio(i))   
    
    for tr = 1:num_trs        
        spk_E = pois_spikes(t(end)+tk(end), dt,N_E,FR_E);
        spk_I = pois_spikes(t(end)+tk(end), dt,N_I,FR_I);

        GE(:,i,tr) = conv(spk_E,kA, 'valid');
        GI(:,i,tr) = conv(spk_I,kG, 'valid')*boost;
        LFP_E(:,i,tr) = detrend(GE(:,i,tr),'constant')*(Ee-Vr);
        LFP_I(:,i,tr) = detrend(GI(:,i,tr),'constant')*(Ei-Vr);
        LFP(:,i,tr) = LFP_E(:,i,tr) + LFP_I(:,i,tr); 
    end % for tr
    
    % grab data to throw into function to compute H
    data2use = squeeze(LFP(:,i,:)); 
    % compute Hurst exponent 
    [H(:,i), NFC, FC] = bfn_mfin_ml(data2use,'filter','haar','lb',[-0.5,0],'ub',[1.5,10],'verbose',false);

end % for i


if MAKE_PLOT
    % plot relationship between mean H across all regions and E:I ratio
    figure; set(gcf,'color','white');
    scatter(log2([1./EI_ratio])',fliplr(nanmean(H(:,:),1))');
    set(gca,'XTick',log2(fliplr(1./EI_ratio)),'XTickLabel',leg_labels);
    ylabel('Hurst exponent (H)'); xlabel('Excitation:Inhibition Ratio');
    rotateXLabels(gca,90);
    grid on;
    
    fname2save = 'EIratio_1to6.pdf';
    print(gcf,'-dpdf','-opengl','-r300',fname2save);
end % if MAKE_PLOT

%% write out result
for i = 1:num_trs
    varnames{i} = sprintf('channel_%03d',i);
end % for i
varnames = [{'EIlabels','EIratio'}, varnames];

tab2write = cell2table([leg_labels', num2cell([EI_ratio', H'])],'VariableNames',varnames);
fname2save = fullfile(outpath,'EIsim_H.csv');
writetable(tab2write, fname2save, 'FileType','text','delimiter',',');
result = tab2write;

end % function A_insilico_1_eisim


%%
function kernel = syn_kernel(t, tau, type)
%function kernel = syn_kernel(t, tau, type)
% given a specific synaptic kernel type and time constant, this returns a
% time series of the kernel that spans the time defined (t) in seconds
%
% t: time vector in seconds (e.g. t=0:0.001:5)
% tau: t_decay or [t_rise t_decay] in seconds
% type: single- or double-exponential, or alpha synapse

if length(tau)==2
    type='2exp';
end % if length(tau)==2
switch type
    case 'exp'
        %single decay exponential (e^-t)
        kernel = exp(-t./tau(1));
        
    case '2exp'
        %double exponential -(e^t/t_rise) + e^t/t_decay
        if length(tau)==2
            % compute the normalization factor
            tpeak = tau(2)*tau(1)/(tau(2)-tau(1))*log(tau(2)/tau(1));
            normf = 1/(-exp(-tpeak./tau(1))+exp(-tpeak./tau(2)));            
            kernel = (-exp(-t./tau(1))+exp(-t./tau(2)))*normf;
        else
            kernel = [];
            disp('Need two time constants for double exponential.')
        end
                
    case 'alpha'
        %alpha kernel (t*e^t)
        kernel = (t/tau).*exp(-t./tau);        
end % switch type

end % function syn_kernel



%%
function discretized = pois_spikes(sim_t, dt, N_neu, FR)
%discretized = pois_spikes(sim_t, dt, N_neu, FR)
%   simulate population spiking of N neurons firing at FR each, return a
%   single spike train that is the total spiking

%mu parameter for exponential distribution
MU = 1./(N_neu*FR);

%draw ISI from exp RV 
ISI = exprnd(MU, [((sim_t+2)/MU) 1]); %pad 2s of sim time
spk_times = cumsum(ISI);
spk_times(spk_times>sim_t)=[];

%discretize
bins = (0:dt:sim_t)+dt/2; %make discretizing bins
discretized = hist(spk_times,bins);
end % function pois_spikes

