% Gather per-molecule values
load('./oh_zangle.mat')
oh = zangle_all;
load('./co_zangle.mat')
co = zangle_all;
clear zangle_all
load('./static.mat')
load('./induced.mat')

% Grab lone molecule identities
load('./lone_ids.mat')

% Construct data structures
oh_zangle_chain = zeros(1,size(lone,2));
co_zangle_chain = zeros(1,size(lone,2));
static_dipole_chain = zeros(1,size(lone,2));
induced_dipole_chain = zeros(1,size(lone,2));
oh_zangle_lone = zeros(1,size(lone,2));
co_zangle_lone = zeros(1,size(lone,2));
static_dipole_lone = zeros(1,size(lone,2));
induced_dipole_lone = zeros(1,size(lone,2));
count_chain = zeros(1,size(lone,2));
count_lone = zeros(1,size(lone,2));

% Collect sum of values in chained and lone categories
for i = 1:1:size(lone,2)
    for j = 1:1:size(lone,1)
        if lone(j,i) == 1 % Chained
            oh_zangle_chain(1,i) = oh_zangle_chain(1,i) + oh(j,i);
            co_zangle_chain(1,i) = co_zangle_chain(1,i) + co(j,i);
            static_dipole_chain(1,i) = static_dipole_chain(1,i) + static(j,i);
            induced_dipole_chain(1,i) = induced_dipole_chain(1,i) + induced(j,i);
            count_chain(1,i) = count_chain(1,i) + 1;
        else % Lone
            oh_zangle_lone(1,i) = oh_zangle_lone(1,i) + oh(j,i);
            co_zangle_lone(1,i) = co_zangle_lone(1,i) + co(j,i);
            static_dipole_lone(1,i) = static_dipole_lone(1,i) + static(j,i);
            induced_dipole_lone(1,i) = induced_dipole_lone(1,i) + induced(j,i);
            count_lone(1,i) = count_lone(1,i) + 1;
        end
    end
end

% Construct data structures for average over cycles
oh_zangle_chain_cycle = zeros(1,1000);
oh_zangle_lone_cycle = zeros(1,1000);
co_zangle_chain_cyle = zeros(1,1000);
co_zangle_lone_cycle = zeros(1,1000);
static_dipole_chain_cycle = zeros(1,1000);
static_dipole_lone_cycle = zeros(1,1000);
induced_dipole_chain_cycle = zeros(1,1000);
induced_dipole_lone_cycle = zeros(1,1000);
count_chain_cycle = zeros(1,1000);
count_lone_cycle = zeros(1,1000);
cycle_contributions = zeros(1,1000);

% Collect sums for average over cycles
ct = 0;
for i = 1:1:size(count_chain,2)
    ct = ct + 1;
    if ct > 1000
        ct = 1;
    end
    oh_zangle_chain_cycle(1,ct) = oh_zangle_chain_cycle(1,ct) + oh_zangle_chain(1,i);
    oh_zangle_lone_cycle(1,ct) = oh_zangle_lone_cycle(1,ct) + oh_zangle_lone(1,i);
    co_zangle_chain_cycle(1,ct) = co_zangle_chain_cycle(1,ct) + co_zangle_chain(1,i);
    co_zangle_lone_cycle(1,ct) = co_zangle_lone_cycle(1,ct) + co_zangle_lone(1,i);
    static_dipole_chain_cycle(1,ct) = static_dipole_chain_cycle(1,ct) + static_dipole_chain(1,i);
    static_dipole_lone_cycle(1,ct) = static_dipole_lone_cycle(1,ct) + static_dipole_lone(1,i);
    induced_dipole_chain_cycle(1,ct) = induced_dipole_chain_cycle(1,ct) + induced_dipole_chain(1,i);
    induced_dipole_lone_cycle(1,ct) = induced_dipole_lone_cycle(1,ct) + induced_dipole_lone(1,i);
    count_chain_cycle(1,ct) = count_chain_cycle(1,ct) + count_chain(1,i);
    count_lone_cycle(1,ct) = count_lone_cycle(1,ct) + count_lone(1,i);
    cycle_contributions(1,ct) = cycle_contributions(1,ct) + 1;
end

% Normalization -> per molecule
oh_zangle_chain_cycle = oh_zangle_chain_cycle./count_chain_cycle;
oh_zangle_lone_cycle = oh_zangle_lone_cycle./count_lone_cycle;
co_zangle_chain_cycle = co_zangle_chain_cycle./count_chain_cycle;
co_zangle_lone_cycle = co_zangle_lone_cycle./count_lone_cycle;
static_dipole_chain_cycle_permol = static_dipole_chain_cycle./count_chain_cycle;
static_dipole_lone_cycle_permol = static_dipole_lone_cycle./count_lone_cycle;
induced_dipole_chain_cycle_permol = induced_dipole_chain_cycle./count_chain_cycle;
induced_dipole_lone_cycle_permol = induced_dipole_lone_cycle./count_lone_cycle;

% Normalization -> Total
static_dipole_chain_cycle = static_dipole_chain_cycle./cycle_contributions;
static_dipole_lone_cycle = static_dipole_lone_cycle./cycle_contributions;
induced_dipole_chain_cycle = induced_dipole_chain_cycle./cycle_contributions;
induced_dipole_lone_cycle = induced_dipole_lone_cycle./cycle_contributions;

% Save values
save('./average_oh_chains.mat','oh_zangle_chain_cycle')
save('./average_co_chains.mat','co_zangle_chain_cycle')
save('./average_static_chains.mat','static_dipole_chain_cycle')
save('./average_induced_chains.mat','induced_dipole_chain_cycle')
save('./average_static_permol_chains.mat','static_dipole_chain_cycle_permol')
save('./average_induced_permol_chains.mat','induced_dipole_chain_cycle_permol')
save('./average_oh_lone.mat','oh_zangle_lone_cycle')
save('./average_co_lone.mat','co_zangle_lone_cycle')
save('./average_static_lone.mat','static_dipole_lone_cycle')
save('./average_induced_lone.mat','induced_dipole_lone_cycle')
save('./average_static_permol_lone.mat','static_dipole_lone_cycle_permol')
save('./average_induced_permol_lone.mat','induced_dipole_lone_cycle_permol')