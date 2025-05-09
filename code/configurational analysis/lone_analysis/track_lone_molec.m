% Lone molecule lifetime analysis
nmol = 500; % Number of molecules in simulation
tmax = 50000; % Number of time steps in simulation

% Collect molecular values
load('./oh_zangle.mat','zangle_all')
oh = zangle_all;
load('./co_zangle.mat','zangle_all')
co = zangle_all;
clear zangle_all
load('./static_dipole.mat')
load('./induced_dipole.mat')
load('./lone_molecules.mat')

% Data structures
key = {'times','oh zangle','co zangle','static dipole','induced dipole'};
lone_history = cell(1,5);

ct = 0; % Count of lone molecule histories recorded
for i = 1:1:nmol
    if lone(i,1) == 1
        fl = 0;
    else
        fl = 1;
        co_zan = co(i,1);
        oh_zan = oh(i,1);
        sta = static(i,1);
        ind = induced(i,1);
        tm = 1;
    end
    for j = 2:1:tmax
        if fl == 0
            if lone(i,j) == 0
                co_zan = co(i,1);
                oh_zan = oh(i,1);
                sta = static(i,1);
                ind = induced(i,1);
                tm = j;
                fl = 1;
            end
        else
            if lone(i,j) == 0
                co_zan = [co_zan,co(i,j)];
                oh_zan = [oh_zan,oh(i,j)];
                sta = [sta,static(i,j)];
                ind = [ind,induced(i,j)];
                tm = [tm,j];
            else
                ct = ct + 1;
                lone_history{ct,1} = tm;
                lone_history{ct,2} = oh_zan;
                lone_history{ct,3} = co_zan;
                lone_history{ct,4} = sta;
                lone_history{ct,5} = ind;
                co_zan = [];
                oh_zan = [];
                sta = [];
                ind = [];
                tm = [];
                fl = 0;
            end
        end
    end
end

% Save lone molecule histories
save('./lone_history.mat','lone_history','key')

% Data structures for raw values
start_co = zeros(ct,1);
end_co = zeros(ct,1);
max_co = zeros(ct,1);
min_co = zeros(ct,1);
start_oh = zeros(ct,1);
end_oh = zeros(ct,1);
max_oh = zeros(ct,1);
min_oh = zeros(ct,1);
start_static = zeros(ct,1);
end_static = zeros(ct,1);
max_static = zeros(ct,1);
min_static = zeros(ct,1);
start_induced = zeros(ct,1);
end_induced = zeros(ct,1);
max_induced = zeros(ct,1);
min_induced = zeros(ct,1);
start_time = zeros(ct,1);
end_time = zeros(ct,1);

% Grab raw values
for i = 1:1:ct
    start_co(i,1) = lone_history{i,3}(1);
    end_co(i,1) = lone_history{i,3}(end);
    max_co(i,1) = max(lone_history{i,3});
    min_co(i,1) = min(lone_hisory{i,3});
    start_oh(i,1) = lone_history{i,2}(1);
    end_oh(i,1) = lone_history{i,2}(end);
    max_oh(i,1) = max(lone_history{i,2});
    min_oh(i,1) = min(lone_history{i,2});
    start_static(i,1) = lone_history{i,4}(1);
    end_static(i,1) = lone_history{i,4}(end);
    max_static(i,1) = max(lone_history{i,4});
    min_static(i,1) = min(lone_history{i,4});
    start_induced(i,1) = lone_history{i,5}(1);
    end_induced(i,1) = lone_history{i,5}(end);
    max_induced(i,1) = max(lone_history{i,5});
    min_induced(i,1) = min(lone_history{i,5});
    start_time(i,1) = lone_history{i,1}(1);
    end_time(i,1) = lone_history{i,1}(end);
end

% Data structures for running sum over all times
start_co_all = zeros(1,max(start_time));
end_co_all = zeros(1,max(end_time));
max_co_all = zeros(1,max(start_time));
min_co_all = zeros(1,max(start_time));
start_oh_all = zeros(1,max(start_time));
end_oh_all = zeros(1,max(end_time));
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
start_count_all = zeros(1,max(start_time));
end_count_all = zeros(1,max(end_time));

% Collect all times
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
    end_static_all(1,end_time(i,1)) = end_static_all(1,end_time(i,1)) + end_time(i,1);
    max_static_all(1,static_time(i,1)) = max_static_all(1,static_time(i,1)) + max_static(i,1);
    min_static_all(1,static_time(i,1)) = min_static_all(1,static_time(i,1)) + min_static(i,1);
    start_induced_all(1,start_time(i,1)) = start_induced_all(1,start_time(i,1)) + start_induced(i,1);
    end_induced_all(1,end_time(i,1)) = end_induced_all(1,end_time(i,1)) + end_induced(i,1);
    max_induced_all(1,start_time(i,1)) = max_induced_all(1,start_time(i,1)) + max_induced(i,1);
    min_induced_all(1,start_time(i,1)) = min_induced_all(1,start_time(i,1)) + min_induced(i,1);
    start_count_all(1,start_time(i,1)) = start_count_all(1,start_time(i,1)) + 1;
    end_count_all(1,end_time(i,1)) = end_count_all(1,end_time(i,1)) + 1;
end

% Data structures for averages in cycles
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
start_count_cycle = zeros(1,1000);
end_count_cycle = zeros(1,1000);

% Average over cycles
ct = 0;
for i = 1:1:size(start_count_all,2)
    ct = ct + 1;
    if ct > 1000
        ct = 1;
    end
    start_oh_cycle(1,ct) = start_oh_cycle(1,ct) + start_oh_all(1,i);
    end_oh_cycle(1,ct) = end_oh_cycle(1,ct) + end_oh_all(1,i);
    max_oh_cycle(1,ct) = max_oh_cycle(1,ct) + max_oh_all(1,i);
    min_oh_cycle(1,ct) = min_oh_cycle(1,ct) + min_oh_all(1,i);
    start_co_cycle(1,ct) = start_co_cycle(1,ct) + start_co_all(1,i);
    end_co_cycle(1,ct) = end_co_cycle(1,ct) + end_co_all(1,i);
    max_co_cycle(1,ct) = max_co_cycle(1,ct) + max_co_all(1,i);
    min_co_cycle(1,ct) = min_co_cycle(1,ct) + min_co_all(1,i);
    start_static_cycle(1,ct) = start_static_cycle(1,ct) + start_static_all(1,i);
    end_static_cycle(1,ct) = end_static_cycle(1,ct) + end_static_all(1,i);
    max_static_cycle(1,ct) = max_static_cycle(1,ct) + max_static_all(1,i);
    min_static_cycle(1,ct) = min_static_cycle(1,ct) + min_static_all(1,i);
    start_induced_cycle(1,ct) = start_induced_cycle(1,ct) + start_induced_all(1,i);
    end_induced_cycle(1,ct) = end_induced_cycle(1,ct) + end_induced_all(1,i);
    max_induced_cycle(1,ct) = max_induced_cycle(1,ct) + max_induced_all(1,i);
    min_induced_cycle(1,ct) = min_induced_cycle(1,ct) + min_induced_all(1,i);
    start_count_cycle(1,ct) = start_count_cycle(1,ct) + start_count_all(1,i);
    end_count_cycle(1,ct) = end_count_cycle(1,ct) + end_count_all(1,i);
end

% Normalize values
start_co_cycle = start_co_cycle./start_count_cycle;
end_co_cycle = end_co_cycle./end_count_cycle;
max_co_cycle = max_co_cycle./start_count_cycle;
min_co_cycle = min_co_cycle./start_count_cycle;
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

% Save averaged values
save('./average_co_lone_molec.mat','start_co_cycle','end_co_cycle','max_co_cycle','min_c_cycle')
save('./average_oh_lone_molec.mat','start_oh_cycle','end_oh_cycle','max_oh_cycle','min_oh_cycle')
save('./average_static_lone_molec.mat','start_static_cycle','end_static_cycle','max_static_cycle','min_static_cycle')
save('./average_induced_lone_molec.mat','start_induced_cycle','end_induced_cycle','max_induced_cycle','min_induced_cycle')