% static ADRF-CP matching condition between 4 protons and one 13C nucleus- reading from simpson spin file


%function ADRF_CPMAS_4protons()

% System specification
[sys,inter]=s2spinach('ADRF_CPMAS_5spins_strongerWd_CH0.2.in');

% Magnet
sys.magnet=18.8;

% Basis set (increase!)
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=5;

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Experiment parameters
parameters.rate=12000;
parameters.axis=[1 1 1];
parameters.max_rank=5;
parameters.spins={'1H','13C'};
parameters.offset=[0 0];
parameters.grid='rep_2ang_100pts_sph';
parameters.ramp_amp=2*pi*linspace(0,1.7e3,400);
parameters.ramp_dt=25e-6;
parameters.cw_amp=2*pi*1000*[0:0.1:1.5 2:0.5:9.5 10:0.1:14 14.5:0.5:21.5 22:0.1:26 26.5:0.5:30]; %linspace(0,50e3,201)
parameters.arrfsteps=600;
parameters.initial_amp=0.9;
%parameters.cw_dur=5e-3;
parameters.sweep=12e3;
parameters.nsteps=1;
parameters.zerofill=1024;
parameters.verbose=1;

% Initial and detection state
parameters.rho0=(state(spin_system,'L+','1H')+...
                 state(spin_system,'L-','1H'))/2;
parameters.coil=(state(spin_system,'L+','13C')+...
                 state(spin_system,'L-','13C'))/2;

% Simulation
mag=singlerot(spin_system,@adrf_cpmasrampup,parameters,'nmr');

% Plotting  
figure(); plot(parameters.cw_amp/(2*pi),real(mag));

%end

