% Collect per-molecule data pre and post removal from chain not dying
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
            for k = 1:1:size(c,2) % Loop through molecules at time j in chain
                fl = 0; % Flag of if molecule k if new
                for m = 1:1:size(d,2) % Loop through molecules at time j-1 in chain
                    if c(1,k) == d(1,m) % If molecule IDs match
                        fl = 1; % Flip flag
                        break; % Break from checking
                    end
                end
                if fl == 0 % If molecule k is new
                    if lone(c(1,k),j) == 0
                        co_pre = [co_pre;co(c(1,k),j-1)];
                        co_post = [co_post;co(c(1,k),j)];
                        oh_pre = [oh_pre;oh(c(1,k),j-1)];
                        oh_post = [oh_post;oh(c(1,k),j)];
                        static_pre = [static_pre;static(c(1,k),j-1)];
                        static_post = [static_post;static(c(1,k),j)];
                        induced_pre = [induced_pre;induced(c(1,k),j-1)];
                        induced_post = [induced_post;induced(c(1,k),j)];
                        times = [times;j];
                    end
                end
            end
        end
    end
end

% Analysis
% All time sums of values data structures
rem_co_all = zeros(1,max(times));
prerem_co_all = zeros(1,max(times));
rem_oh_all = zeros(1,max(times));
prerem_oh_all = zeros(1,max(times));
rem_static_all = zeros(1,max(times));
prerem_static_all = zeros(1,max(times));
rem_induced_all = zeros(1,max(times));
prerem_induced_all = zeros(1,max(times));
rem_count_all = zeros(1,max(times));
% Collect sums at each time of all values
for i = 1:1:size(times,1)
    rem_co_all(1,times(i,1)) = rem_co_all(1,times(i,1)) + co_post(i,1);
    prerem_co_all(1,times(i,1)) = prerem_co_all(1,times(i,1)) + co_pre(i,1);
    rem_oh_all(1,times(i,1)) = rem_oh_all(1,times(i,1)) + oh_post(i,1);
    prerem_oh_all(1,times(i,1)) = prerem_oh_all(1,times(i,1)) + oh_pre(i,1);
    rem_static_all(1,times(i,1)) = rem_static_all(1,times(i,1)) + static_post(i,1);
    prerem_static_all(1,times(i,1)) = prerem_static_all(1,times(i,1)) + static_pre(i,1);
    rem_induced_all(1,times(i,1)) = rem_induced_all(1,times(i,1)) + induced_post(i,1);
    prerem_induced_all(1,times(i,1)) = prerem_induced_all(1,times(i,1)) + induced_pre(i,1);
    rem_count_all(1,times(i,1)) = rem_count_all(1,times(i,1)) + 1;
end
%Averaged over cycles data structures
rem_co_cycle = zeros(1,1000);
prerem_co_cycle = zeros(1,1000);
rem_oh_cycle = zeros(1,1000);
prerem_oh_cycle = zeros(1,1000);
rem_static_cycle = zeros(1,1000);
prerem_static_cycle = zeros(1,1000);
rem_induced_cycle = zeros(1,1000);
prerem_induced_cycle = zeros(1,1000);
rem_count_cycle = zeros(1,1000);
% Average across cycles
ct = 0;
for i = 1:1:size(rem_count_all,2)
    ct = ct + 1;
    if ct > 1000
        ct = 1;
    end
    rem_co_cycle(1,ct) = rem_co_cycle(1,ct) + rem_co_all(1,i);
    prerem_co_cycle(1,ct) = prerem_co_cycle(1,ct) + prerem_co_all(1,i);
    rem_oh_cycle(1,ct) = rem_oh_cycle(1,ct) + rem_oh_all(1,i);
    prerem_oh_cycle(1,ct) = prerem_oh_cycle(1,ct) + prerem_oh_all(1,i);
    rem_static_cycle(1,ct) = rem_static_cycle(1,ct) + rem_static_all(1,i);
    prerem_static_cycle(1,ct) = prerem_static_cycle(1,ct) + prerem_static_all(1,i);
    rem_induced_cycle(1,ct) = rem_induced_cycle(1,ct) + rem_induced_all(1,i);
    prerem_induced_cycle(1,ct) = prerem_induced_cycle(1,ct) + prerem_induced_all(1,i);
    rem_count_cycle(1,ct) = rem_count_cycle(1,ct) + rem_count_all(1,i);
end
% Normalize to per-molecule values in cycle
rem_co_cycle = rem_co_cycle./rem_count_cycle;
prerem_co_cycle = prerem_co_cycle./rem_count_cycle;
rem_oh_cycle = rem_oh_cycle./rem_count_cycle;
prerem_oh_cycle = prerem_oh_cycle./rem_count_cycle;
rem_static_cycle = rem_static_cycle./rem_count_cycle;
prerem_static_cycle = prerem_static_cycle./rem_count_cycle;
rem_induced_cycle = rem_induced_cycle./rem_count_cycle;
prerem_induced_cycle = prestart_induced_cycloe./rem_count_cycle;

% Save values
save('./raw_co_prepost_removal.mat','co_pre','co_post')
save('./raw_oh_prepost_removal.mat','oh_pre','oh_post')
save('./raw_static_prepost_removal.mat','static_pre','static_post')
save('./raw_induced_prepost_removal.mat','induced_pre','induced_post')
save('./times_removal.mat','times')
save('./average_co_prepost_removal.mat','rem_co_cycle','prerem_co_cycle')
save('./average_oh_prepost_removal.mat','rem_oh_cycle','prerem_oh_cycle')
save('./average_static_prepost_removal.mat','rem_static_cycle','prerem_static_cycle')
save('./average_induced_prepost_removal.mat','rem_induced_cycle','prerem_induced_cycle')
