% Collect per-molecule data pre and post chain death
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
load('./lone_ids.mat')

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
    a = importdata([folder,'chain',num2str(i),'_end']); % Get chain end time
    b = importdata([folder,'chain',num2str(i),'_',num2str(a-1),'mol']); % Get molecule ID's at chain start time
    % Check that all molecules are alone after chain death
    alone = zeros(1,size(b,2));
    for q = 1:1:size(b,2)
        if lone(b(1,q),a) == 0
            alone(1,q) = 1;
        end
    end
    if sum(alone) == size(b,2)
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
end

% Analysis
% All time sums of values data structures
die_co_all = zeros(1,max(times));
predie_co_all = zeros(1,max(times));
die_oh_all = zeros(1,max(times));
predie_oh_all = zeros(1,max(times));
die_static_all = zeros(1,max(times));
predie_static_all = zeros(1,max(times));
die_induced_all = zeros(1,max(times));
predie_induced_all = zeros(1,max(times));
die_count_all = zeros(1,max(times));
% Collect sums at each time of all values
for i = 1:1:size(times,1)
    die_co_all(1,times(i,1)) = die_co_all(1,times(i,1)) + co_post(i,1);
    predie_co_all(1,times(i,1)) = predie_co_all(1,times(i,1)) + co_pre(i,1);
    die_oh_all(1,times(i,1)) = die_oh_all(1,times(i,1)) + oh_post(i,1);
    predie_oh_all(1,times(i,1)) = predie_oh_all(1,times(i,1)) + oh_pre(i,1);
    die_static_all(1,times(i,1)) = die_static_all(1,times(i,1)) + static_post(i,1);
    predie_static_all(1,times(i,1)) = predie_static_all(1,times(i,1)) + static_pre(i,1);
    die_induced_all(1,times(i,1)) = die_induced_all(1,times(i,1)) + induced_post(i,1);
    predie_induced_all(1,times(i,1)) = predie_induced_all(1,times(i,1)) + induced_pre(i,1);
    die_count_all(1,times(i,1)) = die_count_all(1,times(i,1)) + 1;
end
%Averaged over cycles data structures
die_co_cycle = zeros(1,1000);
predie_co_cycle = zeros(1,1000);
die_oh_cycle = zeros(1,1000);
predie_oh_cycle = zeros(1,1000);
die_static_cycle = zeros(1,1000);
predie_static_cycle = zeros(1,1000);
die_induced_cycle = zeros(1,1000);
predie_induced_cycle = zeros(1,1000);
die_count_cycle = zeros(1,1000);
% Average across cycles
ct = 0;
for i = 1:1:size(die_count_all,2)
    ct = ct + 1;
    if ct > 1000
        ct = 1;
    end
    die_co_cycle(1,ct) = die_co_cycle(1,ct) + die_co_all(1,i);
    predie_co_cycle(1,ct) = predie_co_cycle(1,ct) + predie_co_all(1,i);
    die_oh_cycle(1,ct) = die_oh_cycle(1,ct) + die_oh_all(1,i);
    predie_oh_cycle(1,ct) = predie_oh_cycle(1,ct) + predie_oh_all(1,i);
    die_static_cycle(1,ct) = die_static_cycle(1,ct) + die_static_all(1,i);
    predie_static_cycle(1,ct) = predie_static_cycle(1,ct) + predie_static_all(1,i);
    die_induced_cycle(1,ct) = die_induced_cycle(1,ct) + die_induced_all(1,i);
    predie_induced_cycle(1,ct) = predie_induced_cycle(1,ct) + predie_induced_all(1,i);
    die_count_cycle(1,ct) = die_count_cycle(1,ct) + die_count_all(1,i);
end
% Normalize to per-molecule values in cycle
die_co_cycle = die_co_cycle./die_count_cycle;
predie_co_cycle = predie_co_cycle./die_count_cycle;
die_oh_cycle = die_oh_cycle./die_count_cycle;
predie_oh_cycle = predie_oh_cycle./die_count_cycle;
die_static_cycle = die_static_cycle./die_count_cycle;
predie_static_cycle = predie_static_cycle./die_count_cycle;
die_induced_cycle = die_induced_cycle./die_count_cycle;
predie_induced_cycle = prestart_induced_cycloe./die_count_cycle;

% Save values
save('./raw_co_prepost_chaindeath.mat','co_pre','co_post')
save('./raw_oh_prepost_chaindeath.mat','oh_pre','oh_post')
save('./raw_static_prepost_chaindeath.mat','static_pre','static_post')
save('./raw_induced_prepost_chaindeath.mat','induced_pre','induced_post')
save('./times_chaindeath.mat','times')
save('./average_co_prepost_chaindeath.mat','die_co_cycle','predie_co_cycle')
save('./average_oh_prepost_chaindeath.mat','die_oh_cycle','predie_oh_cycle')
save('./average_static_prepost_chaindeath.mat','die_static_cycle','predie_static_cycle')
save('./average_induced_prepost_chaindeath.mat','die_induced_cycle','predie_induced_cycle')
