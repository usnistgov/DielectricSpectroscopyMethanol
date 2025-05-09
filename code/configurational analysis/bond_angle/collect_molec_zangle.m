% Bond angle to z axis collection

% Inputs
tmax = 50000; % Number of time steps in simulation
nmol = 500; % Number of molecules in simulation
natom_permol = 6; % Number of atoms in one molecule
atom1 = 1; % Index in first molecule of first atom in bond of interest
atom2 = 2; % Index in first molecule of second atom in bond of interest
file = './methanol.arc'; % Output tinker configuration file

natom = nmol*natom_permol % Number of atoms in simulation

% Incidental data structure(s)
pos = zeros(natom,3); % Positions of atoms at current time step

% Output data structures
zangle_all = zeros(nmol,tmax); % bond angle to z axis of each molecule at each time step

fid = fopen(file,'r'); % Open file for reading

for i = 1:1:tmax % Loop through all time steps in simulation
    % Read in positions at current time
    a = fgetl(fid); % Read header line
    a = fgetl(fid); % Read header line
    for j = 1:1:natom % Loop through all atoms in simulation
        a = fgetl(fid); % Read in jth atom line
        b = strsplit(a,' '); % Parse line by spaces
        pos(j,1) = str2num(b{1,4}); % Parse and assign x position of jth atom
        pos(j,2) = str2num(b{1,5}); % Parse and assign y position of jth atom
        pos(j,3) = str2num(b{1,6}); % Parse and assign z position of jth atom
    end
    % Determine angles of bonds to z axis for each molecule
    for k = 1:1:nmol % Loop through all molecules in simulation
        index = (k-1)*natom_permol; % Index before start of current molecule in position array
        vec = pos(index+atom1,:) - pos(index+atom2,:); % Vector of atom1-atom2 in current molecule
        mag_vec = sqrt(sum(vec.^2)); % Magnitude of vector atom1-atom2
        angle_radian = acos(vec(3)/mag_vec); % Angle atom2-atom1-z-axis in radians
        angle_degree = angle_radian*(180/pi); % Angle in degrees
        zangle_all(k,i) = angle_degree; % Assign angle to array position
    end
end

fclose(fid); % Close tinker output configuration file

save('./zangle.mat','zangle_all') % Save