Collection of molecular statistics before and after chain lifetime events

The matlab code at issue are the following:
(1) collect_prepost_chainform.m
(2) collect_prepost_addition.m
(3) collect_prepost_removal.m
(4) collect_prepost_chaindeath.m
Each of these pieces of code determines the angle of the OH bond to the z-axis, the angle of the CO bond to the z-axis, the static molecular dipole, and the induced molecular dipole before and after (1) chain formation, (2) addition to an already existing chain, (3) removal from a chain not dying, and (4) chain death. This code requires the molecular static and induced dipole at each time step and the angle of the OH and CO bonds to the z-axis (the axis upon which the electric field is applied) for each molecule at each time step as well as the chain tracking data provided by other code. This code then provides averaged values just before and just after each chain event averaged over cycles as is done elsewhere in this work.