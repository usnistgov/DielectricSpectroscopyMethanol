% Collect per-molecule data pre and post addition to extant chain
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
    a = importdata([folder,'chain',num2str(i),'_start']); % Get chain start time
    b = importdata([folder,'chain',num2str(i),'_end']); % Get chain end time
    if a ~= b-1 % Chain must last long enough to have time steps to compare
        for j = a+1:1:b-1 % Loop through n-1 time steps of chain existence
            c = importdata([folder,'chain',num2str(i),'_',num2str(j-1),'mol']); % Load j-1 molecules in chain
            d = importdata([folder,'chain',num2str(i),'_',num2str(j),'mol']); % Load j molecules in chain
            for k = 1:1:size(d,2) % Loop through molecules at time j in chain
                fl = 0; % Flag of if molecule k if new
                for m = 1:1:size(c,2) % Loop through molecules at time j-1 in chain
                    if d(1,k) == c(1,m) % If molecule IDs match
                        fl = 1; % Flip flag
                        break; % Break from checking
                    end
                end
                if fl == 0 % If molecule k is new
                    if lone(d(1,k),j-1) == 0
                        co_pre = [co_pre;co(d(1,k),j-1)];
                        co_post = [co_post;co(d(1,k),j)];
                        oh_pre = [oh_pre;oh(d(1,k),j-1)];
                        oh_post = [oh_post;oh(d(1,k),j)];
                        static_pre = [static_pre;static(d(1,k),j-1)];
                        static_post = [static_post;static(d(1,k),j)];
                        induced_pre = [induced_pre;induced(d(1,k),j-1)];
                        induced_post = [induced_post;induced(d(1,k),j)];
                        times = [times;j];
                    end
                end
            end
        end
    end
end

% Analysis
% All time sums of values data structures
add_co_all = zeros(1,max(times));
preadd_co_all = zeros(1,max(times));
add_oh_all = zeros(1,max(times));
preadd_oh_all = zeros(1,max(times));
add_static_all = zeros(1,max(times));
preadd_static_all = zeros(1,max(times));
add_induced_all = zeros(1,max(times));
preadd_induced_all = zeros(1,max(times));
add_count_all = zeros(1,max(times));
% Collect sums at each time of all values
for i = 1:1:size(times,1)
    add_co_all(1,times(i,1)) = add_co_all(1,times(i,1)) + co_post(i,1);
    preadd_co_all(1,times(i,1)) = preadd_co_all(1,times(i,1)) + co_pre(i,1);
    add_oh_all(1,times(i,1)) = add_oh_all(1,times(i,1)) + oh_post(i,1);
    preadd_oh_all(1,times(i,1)) = preadd_oh_all(1,times(i,1)) + oh_pre(i,1);
    add_static_all(1,times(i,1)) = add_static_all(1,times(i,1)) + static_post(i,1);
    preadd_static_all(1,times(i,1)) = preadd_static_all(1,times(i,1)) + static_pre(i,1);
    add_induced_all(1,times(i,1)) = add_induced_all(1,times(i,1)) + induced_post(i,1);
    preadd_induced_all(1,times(i,1)) = preadd_induced_all(1,times(i,1)) + induced_pre(i,1);
    add_count_all(1,times(i,1)) = add_count_all(1,times(i,1)) + 1;
end
%Averaged over cycles data structures
add_co_cycle = zeros(1,1000);
preadd_co_cycle = zeros(1,1000);
add_oh_cycle = zeros(1,1000);
preadd_oh_cycle = zeros(1,1000);
add_static_cycle = zeros(1,1000);
preadd_static_cycle = zeros(1,1000);
add_induced_cycle = zeros(1,1000);
preadd_induced_cycle = zeros(1,1000);
add_count_cycle = zeros(1,1000);
% Average across cycles
ct = 0;
for i = 1:1:size(add_count_all,2)
    ct = ct + 1;
    if ct > 1000
        ct = 1;
    end
    add_co_cycle(1,ct) = add_co_cycle(1,ct) + add_co_all(1,i);
    preadd_co_cycle(1,ct) = preadd_co_cycle(1,ct) + preadd_co_all(1,i);
    add_oh_cycle(1,ct) = add_oh_cycle(1,ct) + add_oh_all(1,i);
    preadd_oh_cycle(1,ct) = preadd_oh_cycle(1,ct) + preadd_oh_all(1,i);
    add_static_cycle(1,ct) = add_static_cycle(1,ct) + add_static_all(1,i);
    preadd_static_cycle(1,ct) = preadd_static_cycle(1,ct) + preadd_static_all(1,i);
    add_induced_cycle(1,ct) = add_induced_cycle(1,ct) + add_induced_all(1,i);
    preadd_induced_cycle(1,ct) = preadd_induced_cycle(1,ct) + preadd_induced_all(1,i);
    add_count_cycle(1,ct) = add_count_cycle(1,ct) + add_count_all(1,i);
end
% Normalize to per-molecule values in cycle
add_co_cycle = add_co_cycle./add_count_cycle;
preadd_co_cycle = preadd_co_cycle./add_count_cycle;
add_oh_cycle = add_oh_cycle./add_count_cycle;
preadd_oh_cycle = preadd_oh_cycle./add_count_cycle;
add_static_cycle = add_static_cycle./add_count_cycle;
preadd_static_cycle = preadd_static_cycle./add_count_cycle;
add_induced_cycle = add_induced_cycle./add_count_cycle;
preadd_induced_cycle = prestart_induced_cycloe./add_count_cycle;

% Save values
save('./raw_co_prepost_addition.mat','co_pre','co_post')
save('./raw_oh_prepost_addition.mat','oh_pre','oh_post')
save('./raw_static_prepost_addition.mat','static_pre','static_post')
save('./raw_induced_prepost_addition.mat','induced_pre','induced_post')
save('./times_addition.mat','times')
save('./average_co_prepost_addition.mat','add_co_cycle','preadd_co_cycle')
save('./average_oh_prepost_addition.mat','add_oh_cycle','preadd_oh_cycle')
save('./average_static_prepost_addition.mat','add_static_cycle','preadd_static_cycle')
save('./average_induced_prepost_addition.mat','add_induced_cycle','preadd_induced_cycle')
