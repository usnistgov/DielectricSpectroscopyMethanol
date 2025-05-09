% Lone molecule identification
nmol = 500; % Number of molecules in simulation
tmax = 50000; % Number of time steps in simulation
folder = './';

% Data structures
lone = zeros(nmol,tmax); % Flag array for if each molecule is alone or in a chain

% Identify lone molecules
for i = 1:1:tmax % Loop through time steps in simulation
    a = importdata([folder,'hbmol',num2str(i)]);
    T = zeros(nmol,nmol);
    for j = 1:1:size(a,1)
        T(a(j,1),a(j,2)) = 1;
        T(a(j,2),a(j,1)) = 1;
    end
    for k = 1:1:nmol
        if sum(T(:,k)) == 0
            lone(k,i) = 0;
        else
            lone(k,i) = 1;
        end
    end
end

% Save identities
save('./lone_ids.mat','lone')
