% Script to retrieve molecular static and induced dipole moments (z-component only)
natom_permol = 6; % Number of atoms in one molecule
nmol = 500; % Number of molecules in simulation
tmax = 50000; % Maximum number of simulation steps to read in

natom = nmol*natom_permol; % Number of atoms in simulation

static = zeros(nmol,tmax); % Structure to hold static dipole

fid = fopen('./methanol.ustc','r'); % Open raw tinker output file for reading

% Read in atomic dipoles and convert to molecular dipoles (z-component only)
for i = 1:1:tmax % Loop through all steps in simulation
    a = fgetl(fid); % Read in header line
    a = fgetl(fid); % Read in header line
    for j = 1:1:nmol % Loop through all molecules in simulation
        for k = 1:1:natom_permol % Loop through atoms in current molecule
            a = fgetl(fid); % Read in line for atom
            b = strsplit(a,' '); % Parse line by spaces
            static(j,i) = static(j,i) + str2num(b{1,6}); % Parse to number and assign
        end
    end
end

fclose(fid); % Close file

save('./static_molec_dipole.mat','static') % Save static molecular dipole (z-component)

induced = zeros(nmol,tmax); % Structure to hold induced dipole

fid = fopen('./methanol.uind','r'); % Open raw tinker output file for reading

% Read in atomic dipoles and convert to molecular dipoles (z-component
% only)
for i = 1:1:tmax % Loop through all steps in simulation
    a = fgetl(fid); % Read in header line
    a = fgetl(fid); % Read in header line
    for j = 1:1:nmol % Loop through all molecules in simulation
        for k = 1:1:natom_permol % Loop through atoms in current molecule
            a = fgetl(fid); % Read in line for atom
            b = strsplit(a,' '); % Parse line by spaces
            induced(j,i) = induced(j,i) + str2num(b{1,6}); % Parse to number and assign
        end
    end
end

fclose(fid); % Close file

save('./induced_molec_dipole.mat','induced'); % Save induced molecular dipole (z-component)

disp('Complete')
