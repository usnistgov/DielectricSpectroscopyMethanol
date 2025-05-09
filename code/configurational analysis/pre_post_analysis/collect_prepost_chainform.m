% Collect per-molecule data pre and post chain formation
nchain = 100000; % Number of chains to be analyzed
folder = './'; % Folder in which chain files are located

% Load molecular data
load('./static_dipole.mat')
load('./induced_dipole.mat')
load('./oh_zangle.mat','zangle_all')
oh = zangle_all;
load('./co_zangle.mat','zangle_all')
co = zangle_all;
clear zangle_all

% Data structures
oh_pre = [];
oh_post = [];
co_pre = [];
co_post = [];
static_pre = [];
static_post = [];
induced_pre = [];
induced_post = [];
times = [];

% Collection
for i = 1:1:nchain % Loop through all chains specified
    a = importdata([folder,'chain',num2str(i),'_start']); % Get chain start time
    b = importdata([folder,'chain',num2str(i),'_',num2str(a),'mol']); % Get molecule ID's at chain start time
    % Zero out values for determination
    zan_co_pre = 0;
    zan_co_post = 0;
    zan_oh_pre = 0;
    zan_oh_post = 0;
    sta_pre = 0;
    sta_post = 0;
    ind_pre = 0;
    ind_post = 0;
    % Collect values for start and just before start of current chain
    for j = 1:1:size(b,2)
        zan_co_pre = zan_co_pre + co(b(1,j),a-1);
        zan_co_post = zan_co_post + co(b(1,j),a);
        zan_oh_pre = zan_oh_pre + oh(b(1,j),a-1);
        zan_oh_post = zan_oh_post + oh(b(1,j),a);
        sta_pre = sta_pre + static(b(1,j),a-1);
        sta_post = sta_post + static(b(1,j),a);
        ind_pre = ind_pre + induced(b(1,j),a-1);
        ind_post = ind_post + induced(b(1,j),a);
    end
    co_pre = [co_pre;zan_co_pre/size(b,2)];
    co_post = [co_post;zan_co_post/size(b,2)];
    oh_pre = [oh_pre;zan_oh_pre/size(b,2)];
    oh_post = [oh_post;zan_oh_post/size(b,2)];
    static_pre = [static_pre;sta_pre];
    static_post = [static_post;sta_post];
    induced_pre = [induced_pre;ind_pre];
    induced_post = [induced_post;ind_post];
    times = [times;a];
end

% Analysis
% All time sums of values data structures
start_co_all = zeros(1,max(times));
prestart_co_all = zeros(1,max(times));
start_oh_all = zeros(1,max(times));
prestart_oh_all = zeros(1,max(times));
start_static_all = zeros(1,max(times));
prestart_static_all = zeros(1,max(times));
start_induced_all = zeros(1,max(times));
prestart_induced_all = zeros(1,max(times));
start_count_all = zeros(1,max(times));
% Collect sums at each time of all values
for i = 1:1:size(times,1)
    start_co_all(1,times(i,1)) = start_co_all(1,times(i,1)) + co_post(i,1);
    prestart_co_all(1,times(i,1)) = prestart_co_all(1,times(i,1)) + co_pre(i,1);
    start_oh_all(1,times(i,1)) = start_oh_all(1,times(i,1)) + oh_post(i,1);
    prestart_oh_all(1,times(i,1)) = prestart_oh_all(1,times(i,1)) + oh_pre(i,1);
    start_static_all(1,times(i,1)) = start_static_all(1,times(i,1)) + static_post(i,1);
    prestart_static_all(1,times(i,1)) = prestart_static_all(1,times(i,1)) + static_pre(i,1);
    start_induced_all(1,times(i,1)) = start_induced_all(1,times(i,1)) + induced_post(i,1);
    prestart_induced_all(1,times(i,1)) = prestart_induced_all(1,times(i,1)) + induced_pre(i,1);
    start_count_all(1,times(i,1)) = start_count_all(1,times(i,1)) + 1;
end
%Averaged over cycles data structures
start_co_cycle = zeros(1,1000);
prestart_co_cycle = zeros(1,1000);
start_oh_cycle = zeros(1,1000);
prestart_oh_cycle = zeros(1,1000);
start_static_cycle = zeros(1,1000);
prestart_static_cycle = zeros(1,1000);
start_induced_cycle = zeros(1,1000);
prestart_induced_cycle = zeros(1,1000);
start_count_cycle = zeros(1,1000);
% Average across cycles
ct = 0;
for i = 1:1:size(start_count_all,2)
    ct = ct + 1;
    if ct > 1000
        ct = 1;
    end
    start_co_cycle(1,ct) = start_co_cycle(1,ct) + start_co_all(1,i);
    prestart_co_cycle(1,ct) = prestart_co_cycle(1,ct) + prestart_co_all(1,i);
    start_oh_cycle(1,ct) = start_oh_cycle(1,ct) + start_oh_all(1,i);
    prestart_oh_cycle(1,ct) = prestart_oh_cycle(1,ct) + prestart_oh_all(1,i);
    start_static_cycle(1,ct) = start_static_cycle(1,ct) + start_static_all(1,i);
    prestart_static_cycle(1,ct) = prestart_static_cycle(1,ct) + prestart_static_all(1,i);
    start_induced_cycle(1,ct) = start_induced_cycle(1,ct) + start_induced_all(1,i);
    prestart_induced_cycle(1,ct) = prestart_induced_cycle(1,ct) + prestart_induced_all(1,i);
    start_count_cycle(1,ct) = start_count_cycle(1,ct) + start_count_all(1,i);
end
% Normalize to per-molecule values in cycle
start_co_cycle = start_co_cycle./start_count_cycle;
prestart_co_cycle = prestart_co_cycle./start_count_cycle;
start_oh_cycle = start_oh_cycle./start_count_cycle;
prestart_oh_cycle = prestart_oh_cycle./start_count_cycle;
start_static_cycle = start_static_cycle./start_count_cycle;
prestart_static_cycle = prestart_static_cycle./start_count_cycle;
start_induced_cycle = start_induced_cycle./start_count_cycle;
prestart_induced_cycle = prestart_induced_cycloe./start_count_cycle;

% Save values
save('./raw_co_prepost_chainform.mat','co_pre','co_post')
save('./raw_oh_prepost_chainform.mat','oh_pre','oh_post')
save('./raw_static_prepost_chainform.mat','static_pre','static_post')
save('./raw_induced_prepost_chainform.mat','induced_pre','induced_post')
save('./times_chainform.mat','times')
save('./average_co_prepost_chainform.mat','start_co_cycle','prestart_co_cycle')
save('./average_oh_prepost_chainform.mat','start_oh_cycle','prestart_oh_cycle')
save('./average_static_prepost_chainform.mat','start_static_cycle','prestart_static_cycle')
save('./average_induced_prepost_chainform.mat','start_induced_cycle','prestart_induced_cycle')