% Values over chain lifetimes analysis
folder = './'; % Folder containing chain data
nchain = 100000; % Number of chains to be analyzed

% Collect molecular values
load('./oh_zangle.mat','zangle_all')
oh = zangle_all;
load('./co_zangle.mat','zangle_all')
co = zangle_all;
clear zangle_all
load('./static_dipole.mat')
load('./induced_dipole.mat')

% Data structures
key = {'times','lengths','oh zangle','co zangle','static dipole','induced dipole'};
chain_history = cell(nchain,6);

% Collect chain values
for i = 1:1:nchain % Loop through chains to analyze
    a = importdata([folder,'chain',num2str(i),'_start']); % Get chain start time
    b = importdata([folder,'chain',num2str(i),'_end']); % Get chain end time
    for j = a:1:b-1 % Loop through times chain exists
        c = importdata([folder,'chain',num2str(i),'_',num2str(j),'mol']); % Get molecule IDs of chain at time j
        chain_history{i,1} = [chain_history{i,1},j]; % Time index
        chain_history{i,2} = [chain_history{i,2},size(c,2)]; % Length at time j
        % Zero values for molecular collection
        zan_oh = 0;
        zan_co = 0;
        sta = 0;
        ind = 0;
        % Molecular collection
        for k = 1:1:size(c,2) % Loop through molecules in chain i at time j
            zan_oh = zan_oh + oh(c(1,k),j);
            zan_co = zan_co + co(c(1,k),j);
            sta = sta + static(c(1,k),j);
            ind = ind + induced(c(1,k),j);
        end
        % Add molecular collections to array
        chain_history{i,3} = [chain_history{i,3},zan_oh/size(c,2)];
        chain_history{i,4} = [chain_history{i,4},zan_co/size(c,2)];
        chain_history{i,5} = [chain_history{i,5},sta];
        chain_history{i,6} = [chain_history{i,6},ind];
    end
end

% Save chain history
save('./chain_history.mat','chain_history','key')

% Analysis

% Initial data structures for all chains
start_co = zeros(nchain,1);
end_co = zeros(nchain,1);
max_co = zeros(nchain,1);
min_co = zeros(nchain,1);
start_oh = zeros(nchain,1);
end_oh = zeros(nchain,1);
max_oh = zeros(nchain,1);
min_oh = zeros(nchain,1);
start_static = zeros(nchain,1);
end_static = zeros(nchain,1);
max_static = zeros(nchain,1);
min_static = zeros(nchain,1);
start_induced = zeros(nchain,1);
end_induced = zeros(nchain,1);
max_induced = zeros(nchain,1);
min_induced = zeros(nchain,1);
start_length = zeros(nchain,1);
end_length = zeros(nchain,1);
max_length = zeros(nchain,1);
min_length = zeros(nchain,1);
start_time = zeros(nchain,1);
end_time = zeros(nchain,1);

% Collect each chain information
for i = 1:1:nchain
    start_co(i,1) = chain_history{i,4}(1);
    end_co(i,1) = chain_history{i,4}(end);
    max_co(i,1) = max(chain_history{i,4});
    min_co(i,1) = min(chain_history{i,4});
    start_oh(i,1) = chain_history{i,3}(1);
    end_oh(i,1) = chain_history{i,3}(end);
    max_oh(i,1) = max(chain_history{i,3});
    min_oh(i,1) = min(chain_history{i,3});
    start_static(i,1) = chain_history{i,5}(1);
    end_static(i,1) = chain_history{i,5}(end);
    max_static(i,1) = max(chain_history{i,5});
    min_static(i,1) = min(chain_history{i,5});
    start_induced(i,1) = chain_history{i,6}(1);
    end_induced(i,1) = chain_history{i,6}(end);
    max_induced(i,1) = max(chain_history{i,6});
    min_induced(i,1) = min(chain_history{i,6});
    start_length(i,1) = chain_history{i,2}(1);
    end_length(i,1) = chain_history{i,2}(end);
    max_length(i,1) = max(chain_history{i,2});
    min_length(i,1) = min(chain_history{i,2});
    start_time(i,1) = chain_history{i,1}(1);
    end_time(i,1) = chain_history{i,1}(end);
end

% Data structures for sums at each time of all values
start_co_all = zeros(1,max(start_time));
end_co_all = zeros(1,max(end_time));
max_co_all = zeros(1,max(start_time));
min_co_all = zeros(1,max(start_time));
start_oh_all = zeros(1,max(start_time));
end_oh_all = zeros(1,max(start_time));
max_oh_all = zeros(1,max(start_time));
min_oh_all = zeros(1,max(start_time));
start_static_all = zeros(1,max(start_time));
end_static_all = zeros(1,max(end_time));
max_static_all = zeros(1,max(start_time));
min_static_all = zeros(1,max(start_time));
start_induced_all = zeros(1,max(start_time));
end_induced_all = zeros(1,max(end_time));
max_induced_all = zeros(1,max(start_time));
min_induced_all = zeros(1,max(start_time));
start_length_all = zeros(1,max(start_time));
end_length_all = zeros(1,max(end_time));
max_length_all = zeros(1,max(start_time));
min_length_all = zeros(1,max(start_time));
start_count_all = zeros(1,max(start_time));
end_count_all = zeros(1,max(end_time));

% Collect sum at all times of values
for i = 1:1:size(start_time,1)
    start_co_all(1,start_time(i,1)) = start_co_all(1,start_time(i,1)) + start_co(i,1);
    end_co_all(1,end_time(i,1)) = end_co_all(1,end_time(i,1)) + end_co(i,1);
    max_co_all(1,start_time(i,1)) = max_co_all(1,start_time(i,1)) + max_co(i,1);
    min_co_all(1,start_time(i,1)) = min_co_all(1,start_time(i,1)) + min_co(i,1);
    start_oh_all(1,start_time(i,1)) = start_oh_all(1,start_time(i,1)) + start_oh(i,1);
    end_oh_all(1,end_time(i,1)) = end_oh_all(1,end_time(i,1)) + end_oh(i,1);
    max_oh_all(1,start_time(i,1)) = max_oh_all(1,start_time(i,1)) + max_oh(i,1);
    min_oh_all(1,start_time(i,1)) = min_oh_all(1,start_time(i,1)) + min_oh(i,1);
    start_static_all(1,start_time(i,1)) = start_static_all(1,start_time(i,1)) + start_static(i,1);
    end_static_all(1,end_time(i,1)) = end_static_all(1,end_time(i,1)) + end_static(i,1);
    max_static_all(1,start_time(i,1)) = max_static_all(1,start_time(i,1)) + max_static(i,1);
    min_static_all(1,start_time(i,1)) = min_static_all(1,start_time(i,1)) + min_static(i,1);
    start_induced_all(1,start_time(i,1)) = start_induced_all(1,start_time(i,1)) + start_induced(i,1);
    end_induced_all(1,end_time(i,1)) = end_induced_all(1,end_time(i,1)) + end_induced(i,1));
    max_induced_all(1,start_time(i,1)) = max_induced_all(1,start_time(i,1)) + max_induced(i,1);
    min_induced_all(1,start_time(i,1)) = min_induced_all(1,start_time(i,1)) + min_induced(i,1);
    start_length_all(1,start_time(i,1)) = start_length_all(1,start_time(i,1)) + start_length(i,1);
    end_length_all(1,end_time(i,1)) = end_length_all(1,end_time(i,1)) + end_length(i,1);
    max_length_all(1,start_time(i,1)) = max_length_all(1,start_time(i,1)) + max_length(i,1);
    min_length_all(1,start_time(i,1)) = min_length_all(1,start_time(i,1)) + min_length(i,1);
    start_count_all(1,start_time(i,1)) = start_count_all(1,start_time(i,1)) + 1;
    end_count_all(1,end_time(i,1)) = end_count_all(1,end_time(i,1)) + 1;
end

% Data structures for average over cycles
start_co_cycle = zeros(1,1000);
end_co_cycle = zeros(1,1000);
max_co_cycle = zeros(1,1000);
min_co_cycle = zeros(1,1000);
start_oh_cycle = zeros(1,1000);
end_oh_cycle = zeros(1,1000);
max_oh_cycle = zeros(1,1000);
min_oh_cycle = zeros(1,1000);
start_static_cycle = zeros(1,1000);
end_static_cycle = zeros(1,1000);
max_static_cycle = zeros(1,1000);
min_static_cycle = zeros(1,1000);
start_induced_cycle = zeros(1,1000);
end_induced_cycle = zeros(1,1000);
max_induced_cycle = zeros(1,1000);
min_induced_cycle = zeros(1,1000);
start_length_cycle = zeros(1,1000);
end_length_cycle = zeros(1,1000);
max_length_cycle = zeros(1,1000);
min_length_cycle = zeros(1,1000);
start_count_cycle = zeros(1,1000);
end_count_cycle = zeros(1,1000);

% Average across cycles
ct = 0;
for i = 1:1:size(end_count_all,2)
    ct = ct + 1;
    if ct > 1000
        ct = 1;
    end
    if i <= max(start_time)
        start_co_cycle(1,ct) = start_co_cycle(1,ct) + start_co_all(1,i);
        max_co_cycle(1,ct) = max_co_cycle(1,ct) + max_co_all(1,i);
        min_co_cycle(1,ct) = min_co_cycle(1,ct) + min_co_all(1,i);
        start_oh_cycle(1,ct) = start_oh_cycle(1,ct) + start_oh_all(1,i);
        max_oh_cycle(1,ct) = max_oh_cycle(1,ct) + max_oh_all(1,i);
        min_oh_cycle(1,ct) = min_oh_cycle(1,ct) + min_oh_all(1,i);
        start_static_cycle(1,ct) = start_static_cycle(1,ct) + start_static_all(1,i);
        max_static_cycle(1,ct) = max_static_cycle(1,ct) + max_static_all(1,i);
        min_static_cycle(1,ct) = min_static_cycle(1,ct) + min_static_all(1,i);
        start_induced_cycle(1,ct) = start_induced_cycle(1,ct) + start_induced_all(1,i);
        max_induced_cycle(1,ct) = max_induced_cycle(1,ct) + max_induced_all(1,i);
        min_induced_cycle(1,ct) = min_induced_cycle(1,ct) + min_induced_all(1,i);
        start_length_cycle(1,ct) = start_length_cycle(1,ct) + start_length_all(1,i);
        max_length_cycle(1,ct) = max_length_cycle(1,ct) + max_length_all(1,i);
        min_length_cycle(1,ct) = min_length_cycle(1,ct) + min_length_all(1,i);
        start_count_cycle(1,ct) = start_count_cycle(1,ct) + start_count_all(1,i);
    end
    end_co_cycle(1,ct) = end_co_cycle(1,ct) + end_co_all(1,i);
    end_oh_cycle(1,ct) = end_oh_cycle(1,ct) + end_oh_all(1,i);
    end_static_cycle(1,ct) = end_static_cycle(1,ct) + end_static_all(1,i);
    end_induced_cycle(1,ct) = end_induced_cycle(1,ct) + end_induced_all(1,i);
    end_length_cycle(1,ct) = end_length_cycle(1,ct) + end_length_all(1,i);
    end_count_cycle(1,ct) = end_count_cycle(1,ct) + end_count_all(1,i);
end

% Normalize
start_co_cycle = start_co_cycle./start_count_cycle;
end_co_cycle = end_co_cycle./end_count_cycle;
max_co_cycle = max_co_cycle./start_count_cycle;
min_co_cycle = min_co_cycle.start_count_cycle;
start_oh_cycle = start_oh_cycle./start_count_cycle;
end_oh_cycle = end_oh_cycle./end_count_cycle;
max_oh_cycle = max_oh_cycle./start_count_cycle;
min_oh_cycle = min_oh_cycle./start_count_cycle;
start_static_cycle = start_static_cycle./start_count_cycle;
end_static_cycle = end_static_cycle./end_count_cycle;
max_static_cycle = max_static_cycle./start_count_cycle;
min_static_cycle = min_static_cycle./start_count_cycle;
start_induced_cycle = start_induced_cycle./start_count_cycle;
end_induced_cycle = end_induced_cycle./end_count_cycle;
max_induced_cycle = max_induced_cycle./start_count_cycle;
min_induced_cycle = min_induced_cycle./start_count_cycle;
start_length_cycle = start_length_cycle./start_count_cycle;
end_length_cycle = end_length_cycle./end_count_cycle;
max_length_cycle = max_length_cycle./start_count_cycle;
min_length_cycle = min_length_cycle./start_count_cycle;

% Save averaged data
save('./average_co_zangle_chain.mat','start_co_cycle','end_co_cycle','max_co_cycle','min_co_cycle')
save('./average_oh_zangle_chain.mat','start_oh_cycle','end_oh_cycle','max_oh_cycle','min_oh_cycle')
save('./average_static_chain.mat','start_static_cycle','end_static_cycle','max_static_cycle','min_static_cycle')
save('./average_induced_chain.mat','start_induced_cycle','end_induced_cycle','max_induced_cycle','min_induced_cycle')
save('./average_length_chain.mat','start_length_cycle','end_length_cycle','max_length_cycle','min_length_cycle')