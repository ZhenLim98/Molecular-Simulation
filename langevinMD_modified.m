function [ nSteps, t, E_t, T_t, R2_t, D_t Fmax_t Frms_t] = langevinMD( nAtoms, gamma, targetTemp, timestep, equilTime, runTime, mdmode, startconf)
%LANGEVINMD Perform Langevin Molecular Dynamics on a strand of polymer
%
%   Input variables:
%   nAtoms          Number of atoms
%   gamma           Damping rate in Langevin dynamics in fs^-1
%   targetTemp      Temperature of atoms in kelvin
%   timestep        Molecular dynamics time step in fs
%   equilTime       Time for equilibrating the system in fs
%   runTime         Time over which statistics are taken in fs
%
%   mdmode:
%    'constantE'  : after the equilibration, the friction and random forces are 
%                   turned off and so energy is conserved if the time step is small enough
%    'constantT'  : Langevin dynamics (default)
%    'relax'      : Steepest descent relaxation
%
%   startconf:
%    'new'        : Construct positions and velocities of the polymer from scratch
%    'restart'    : Read positions, velocities, and number of atoms in the polymer from file 'restart_file' 
%
%   Output variables:
%   nSteps          The number of MD steps used for the statistics
%   t               The time at each MD step
%   E_t             The energy at each MD step
%   T_t             The temperature at each MD step
%   P_t             The potential energy at each MD step
%   K_t             The kinetic energy at each MD step
%   R2_t            The square end-to-end length at each MD step
%   rms             The root-mean-squared displacement of the polymer's center of mass.
%   Fmax_t          The maximum force component
%   Frms_t          The root mean squared force component
%
%   This function performs Langevin Dynamics on a single strand of a model
%   polymer in a solvent. 
%
%   The units used in the program are:
%   Energy      eV
%   Force       eV/Angstrom
%   Distance    Angstrom
%   Time        fs
%   Mass        eV fs^2/A^2 = 0.00964855 AMU
%
intalgo = 1;   %Specifies the integration algorithm for the equations of motion
%
% The following potential parameters are taken from C. Liu and M.
% Muthukumar, "Langevin dynamics simulations of early-stage polymer
% nucleation and crystallization", J. Chem. Phys., 109, 2536 (1998)
sigma = 4.5;          % Angstroms
epsilon = 0.00485678; % eV; 0.112 kcal/mol
bondLength = 1.53;    % Angstroms
Kr = 15.18; % eV/A^2; 350 kcal/mol A^2
run_once = false; % the verlet algorithm needs a step back in time
run_count = 0;
previous_pos_2 = 0;
%
% When relaxing the structure the timestep is scaled by this factor: 
%
timestep_factor = 100.0;
%
% 1 AMU = 1.6605e-27 kg ; 1 eV fs^2/A^2 = 1.602e-29 kg
% 1 AMU = 103.65 eV fs^2/A^2 ; 14 AMU = 1451 eV fs^2/A^2
atomicmass = 1451.0;   % eV fs^2/A^2; mass of CH2 = 14 AMU
%
% collisiontime = atomicmass/gamma
% time between collisions in fs
%
% Define other useful quantities
kB = 8.617e-5;        % Boltzmann's constant
kT = kB*targetTemp;   % eV
beta = 1.0/kT;        % 1/eV
%
movieInterval = 100.0; % Time, in fs, between frames printed out to the xyz file. 
%
% Set some defaults in case mdmode and/or startconf aren't specified.
%
switch nargin
    case 6
      mdmode='constantE';
      startconf='new';
    case 7
      startconf='restart';
end
% Initialise the random number generator
rng('shuffle');
%
% Initialise the atomic positions
%
pos = zeros(nAtoms, 3);
vel = zeros(nAtoms, 3);
%
switch startconf
    case 'restart'
%
        fprintf ('\nReading input configuration from file.\n');
        inputFid = fopen('restart_file', 'r');
        nAtoms=fscanf(inputFid, '%i', 1);
        size = [nAtoms 3];
        pos=fscanf(inputFid, '%g', size);        
        vel=fscanf(inputFid, '%g', size);
        fclose(inputFid);
%
    otherwise
     fprintf ('\nConstructing new input configuration.\n');
     pos(:,3) = linspace(0, (nAtoms-1)*bondLength, nAtoms);
% Initialise the atomic velocities according to the Boltzmann distribution.
% The temperature has been increased above the target value to allow for
% the initial relaxation.
      vel = normrnd(0.0, sqrt(1.5*kT/atomicmass), nAtoms, 3);
end
fprintf ('\nConfiguration contains %i atoms.\n', nAtoms);
%
% Translate to the origin
%
cmpos=mean(pos);
pos(:,1) = pos(:,1) - cmpos(1);
pos(:,2) = pos(:,2) - cmpos(2);
pos(:,3) = pos(:,3) - cmpos(3);
%
% Calculate initial potential energy
%
pot1_i = getBondingEnergy()/nAtoms;
pot2_i = getNonBondingEnergy()/nAtoms;
pot_i = pot1_i + pot2_i;
%
% Remove centre of mass velocity
%
cmvel = mean(vel);
vel(:,1) = vel(:,1) - cmvel(1);
vel(:,2) = vel(:,2) - cmvel(2);
vel(:,3) = vel(:,3) - cmvel(3);
%
% Calculate initial temperature and kinetic energy
%
kinstart = getKinetic();
temperature = getTemperature();  
fprintf ('\nStarting configuration:');
fprintf ('\n-----------------------');
fprintf ('\nTemperature = %f K \n', temperature);
%
% Initialize the atomic forces
%
force = zeros(nAtoms, 3);
force = getForce();
%
fprintf ('Bonding energy = %g eV/atom \n', pot1_i);
fprintf ('Non-bonding energy = %g eV/atom \n', pot2_i);
fprintf ('Total potential energy = %g eV/atom \n', pot_i);
Fmax_i = max(max(abs(force)));
Frms_i = sqrt(sum(dot(force,force))/(3*nAtoms));
fprintf ('Rms force component = %g eV/A \n', Frms_i);
fprintf ('Max force component = %g eV/A \n', Fmax_i);
%
% Equilibrate the system
%
if equilTime > 0.1
  nSteps = round(equilTime/timestep)+1;
  for i=1:nSteps
    moveAtoms();
  end
  temperature = getTemperature();
  fprintf ('\nAfter thermalization:');
  fprintf ('\n---------------------');
  fprintf ('\nTemperature = %f K \n', temperature);
  pot1_i = getBondingEnergy()/nAtoms;
  pot2_i = getNonBondingEnergy()/nAtoms;
  pot_i = pot1_i + pot2_i;  
  fprintf ('Bonding energy = %g eV/atom \n', pot1_i);
  fprintf ('Non-bonding energy = %g eV/atom \n', pot2_i);
  fprintf ('Total potential energy = %g eV/atom \n', pot_i);
  Fmax_i = max(max(abs(force)));
  Frms_i = sqrt(sum(dot(force,force))/(3*nAtoms));
  fprintf ('Rms force component = %g eV/A \n', Frms_i);
  fprintf ('Max force component = %g eV/A \n', Fmax_i);     
else
  fprintf ('\nNo thermalization.\n');
end
%
switch mdmode
  case 'constantE'
%   turn off the thermostat
   gamma = 0.0;
   fprintf ('\nConstant energy mode\n');
  case 'relax'
   gamma = 0.0;
   fprintf('\nEnergy minimization mode\n');
   vel = zeros(nAtoms, 3);
  otherwise
   fprintf ('\nConstant temperature mode\n');
end
%   translate to the origin
   cmpos = mean(pos);
   pos(:,1) = pos(:,1) - cmpos(1);
   pos(:,2) = pos(:,2) - cmpos(2);
   pos(:,3) = pos(:,3) - cmpos(3);
%   Remove center of mass velocity
   cmvel = mean(vel);
   vel(:,1) = vel(:,1) - cmvel(1);
   vel(:,2) = vel(:,2) - cmvel(2);
   vel(:,3) = vel(:,3) - cmvel(3);
%
% Calculate how many MD/relaxation steps to do
%
nSteps = round(runTime/timestep)+1;
%
% Initialise some arrays for taking statistics
%
t = linspace(0.0, 0.0, nSteps);
E_t = linspace(0.0, 0.0, nSteps);
T_t = linspace(0.0, 0.0, nSteps);
R2_t = linspace(0.0, 0.0, nSteps);
P_t = linspace(0.0, 0.0, nSteps);
K_t = linspace(0.0, 0.0, nSteps);
D_t = linspace(0.0, 0.0, nSteps);
Fmax_t = linspace(0.0, 0.0, nSteps);
Frms_t = linspace(0.0, 0.0, nSteps);
%
%
tNextFrame = 0.0;
filename1 = strcat('movie_',int2str(nAtoms),'_',int2str(targetTemp),'.xyz');
movieFid = fopen(filename1, 'w');
%
%
%
for i=1:nSteps
    switch mdmode
        case 'relax'
         relaxAtoms();
        otherwise
          moveAtoms();
    end
    t(i) = (i-1)*timestep;
    pot1 = getBondingEnergy();
    pot2 = getNonBondingEnergy();
    P_t(i) = pot1 + pot2;
    K_t(i) = getKinetic();
    T_t(i) = getTemperature();
    E_t(i) = P_t(i) + K_t(i);
    D_t(i) = sqrt(dot(mean(pos),mean(pos)));
    R2_t(i) = (pos(nAtoms,1)-pos(1,1))*(pos(nAtoms,1)-pos(1,1))...
            + (pos(nAtoms,2)-pos(1,2))*(pos(nAtoms,2)-pos(1,2)) ...
            + (pos(nAtoms,3)-pos(1,3))*(pos(nAtoms,3)-pos(1,3));
    Fmax_t(i) = max(max(abs(force)));
    Frms_t(i) = sqrt(sum(dot(force,force))/(3*nAtoms));
%
% Write out the configuration to an xyz file:
%
    if t(i) >= tNextFrame
        fprintf(movieFid, '%i\n\n', nAtoms);
        fprintf(movieFid, 'C %14.7g %14.7g %14.7g\n', pos');
        tNextFrame = tNextFrame + movieInterval;
    end
end
fclose(movieFid);
%
% Write final configuration to 'restart_file'
%
outputFid = fopen('restart_file', 'w');
fprintf(outputFid, '%i', nAtoms);
fprintf(outputFid, '%20.10g', pos);
fprintf(outputFid, '%20.10g', vel);
fclose(outputFid);
%
% Write data to 'output_data' or 'relax_output_data'
%
switch mdmode
    case 'relax'
     filename2 = strcat('relax_output_data.csv');
     dataFid = fopen(filename2, 'a');   
     pot1_f = getBondingEnergy()/nAtoms;
     pot2_f = getNonBondingEnergy()/nAtoms;
     pot_f = pot1_f + pot2_f;   
     
     fprintf ('\nAfter relaxation:');
     fprintf ('\n-----------------');
     fprintf ('\nBonding energy = %g eV/atom \n', pot1_f);
     fprintf ('Non-bonding energy = %g eV/atom \n', pot2_f);
     fprintf ('Total potential energy = %g eV/atom \n', pot_f);
     force = getForce();
     Fmax_f = max(max(abs(force)));
     Frms_f = sqrt(sum(dot(force,force))/(3*nAtoms));
     fprintf ('Rms force component = %g eV/A \n', Frms_f);
     fprintf ('Max force component = %g eV/A \n', Fmax_f);
     fprintf(dataFid, ' \n %i, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g', ... 
       nAtoms, pot_i, pot_f, Fmax_i, Fmax_f, Frms_i, Frms_f, sqrt(R2_t(1)), sqrt(R2_t(nSteps)));
%
     fprintf ('\nChange in bonding energy = %g eV/atom \n', pot1_f-pot1_i);
     fprintf ('Change in non-bonding energy = %g eV/atom \n', pot2_f-pot2_i);
     fprintf ('Total change in potential energy = %g eV/atom \n', pot_f-pot_i);
     fprintf ('Change in rms force = %g eV/atom \n', Frms_f-Frms_i);
     fprintf ('Change in max force = %g eV/atom \n', Fmax_f-Fmax_i);
    otherwise
      filename2 = strcat('output_data.csv');
      dataFid = fopen(filename2, 'a');           
      fprintf(dataFid, ' \n %i, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g', ... 
            nAtoms, targetTemp, mean(T_t), mean(E_t), mean(P_t), mean(K_t) ...
            ,D_t(nSteps), mean(R2_t), sqrt(mean(R2_t)), mean(R2_t)/nAtoms...
            ,std(T_t), std(E_t), std(R2_t));
      fprintf ('\nAfter production run:'); 
      fprintf ('\n---------------------');
      fprintf ('\nAverage temperature = %f K \n', mean(T_t));
      fprintf ('Standard deviation of temperature = %f K \n', std(T_t));
      fprintf ('\nAverage energy = %g eV \n', mean(E_t));
      fprintf ('Standard deviation of energy = %g eV \n', std(E_t));
      fprintf ('\nAverage potential energy = %g eV \n', mean(P_t));
      fprintf ('Average kinetic energy = %g eV \n', mean(K_t));
      fprintf ('\nDistance travelled by center of mass = %g A \n', D_t(nSteps));
      fprintf ('Mean squared end-to-end distance = %g A^2 \n', mean(R2_t));
      fprintf ('RMS end-to-end distance = %g A \n', sqrt(mean(R2_t)));
      fprintf ('\nCenter of mass position = %f %f %f A \n', mean(pos));
      fprintf ('Center of mass velocity = %f %f %f A/fs \n\n', mean(vel));
end
fclose(dataFid);

% -------------------------------
% Nested functions start here ...
% -------------------------------
    function kinetic = getKinetic ()
% Computes the kinetic energy from the velocities and masses
      kinetic = 0.5*atomicmass*(dot(vel(:,1), vel(:,1)) ...
          + dot(vel(:,2), vel(:,2)) + dot(vel(:,3), vel(:,3)));
    end
%
    function potential1 = getBondingEnergy()  
  % Compute the bonding contribution to the potential energy which
  % is simply a harmonic bond-stretching term
        potential1 = 0.0;
        for i1=1:nAtoms-1
            j1=i1+1;
            r = sqrt((pos(i1,1)-pos(j1,1))*(pos(i1,1)-pos(j1,1)) ...
                + (pos(i1,2)-pos(j1,2))*(pos(i1,2)-pos(j1,2)) ...
                + (pos(i1,3)-pos(j1,3))*(pos(i1,3)-pos(j1,3)));
            potential1 = potential1 + Kr*(r-bondLength)*(r-bondLength);
        end
    end
%
    function potential2 = getNonBondingEnergy ()
        %   Computes the non-bonding Lennard-Jones energy term
        %   between remote atoms (third neighbour or further).
        potential2 = 0.0;
        for i1=1:nAtoms-3
            for j1=i1+3:nAtoms
                r2 = (pos(i1,1)-pos(j1,1))*(pos(i1,1)-pos(j1,1)) ...
                    + (pos(i1,2)-pos(j1,2))*(pos(i1,2)-pos(j1,2)) ...
                    + (pos(i1,3)-pos(j1,3))*(pos(i1,3)-pos(j1,3));
                x2 = sigma*sigma/r2;
                x6 = x2*x2*x2;
                x12 = x6*x6;
                potential2 = potential2 + epsilon*(x12-2.0*x6);
            end
        end
    end

   
    function temperature = getTemperature ()
        %GETTEMPERATURE computes the temperature
        vcm = mean(vel);
        temperature = atomicmass*(dot(vel(:,1)-vcm(1), ....
            vel(:,1)-vcm(1)) + dot(vel(:,2)-vcm(2), ...
            vel(:,2)-vcm(2)) + dot(vel(:,3)-vcm(3), ...
            vel(:,3)-vcm(3)))/(3.0*kB*(nAtoms-1));
    end

    function F = getForce ()
        %GETFORCE computes the forces on the atoms
        %   The forces are computed from the gradients of the potential
        %   energy. The potential energy consists of two components: a
        %   bonding term due to bond stretching, and a non-bonding
        %   Lennard-Jones term between remote atoms (third neighbour or
        %   further). Stochastic forces included, but NOT the friction
        %   forces
        
        % Compute the bonding contribution to the atomic forces
        F = zeros(nAtoms, 3);
        for i2=1:nAtoms-1
            j2=i2+1;
            dPos(1) = pos(j2,1)-pos(i2,1);
            dPos(2) = pos(j2,2)-pos(i2,2);
            dPos(3) = pos(j2,3)-pos(i2,3);
            r = sqrt(dPos(1)*dPos(1) + dPos(2)*dPos(2) + dPos(3)*dPos(3));
            fac = 2.0*Kr*(1.0-bondLength/r);
            F(i2,1) = F(i2,1) + fac*dPos(1);
            F(i2,2) = F(i2,2) + fac*dPos(2);
            F(i2,3) = F(i2,3) + fac*dPos(3);
            F(j2,1) = F(j2,1) - fac*dPos(1);
            F(j2,2) = F(j2,2) - fac*dPos(2);
            F(j2,3) = F(j2,3) - fac*dPos(3);
        end
        
        % Add the non-bonding contribution to the atomic forces
        fac1 = -12.0*epsilon/(sigma*sigma);
        for i2=1:nAtoms-3
            for j2=i2+3:nAtoms
                dPos(1) = pos(j2,1)-pos(i2,1);
                dPos(2) = pos(j2,2)-pos(i2,2);
                dPos(3) = pos(j2,3)-pos(i2,3);
                r2 = dPos(1)*dPos(1) + dPos(2)*dPos(2) + dPos(3)*dPos(3);
                x2 = sigma*sigma/r2;
                x6 = x2*x2*x2;
                x12 = x6*x6;
                fac = fac1*x2*(x12-x6);
                F(i2,1) = F(i2,1) + fac*dPos(1);
                F(i2,2) = F(i2,2) + fac*dPos(2);
                F(i2,3) = F(i2,3) + fac*dPos(3);
                F(j2,1) = F(j2,1) - fac*dPos(1);
                F(j2,2) = F(j2,2) - fac*dPos(2);
                F(j2,3) = F(j2,3) - fac*dPos(3);
            end
        end
       
    end
%
    function moveAtoms ()
        %MOVEATOMS updates the atomic positions,velocities and forces
        %   using the OVRVO integrator of Sivak, Chodera, and Crooks 
        %   ( J. Phys. Chem. B, 118, 6466, 2014). Reduces to 
        %   velocity Verlet algorithm if b=1 and gamma=0.
        %   the atomic positins and velocities. The friction forces
        %   included in the velocity update.
        
        % Update positions
        if intalgo == 1
            if run_once == false
                %velocity at next timestep uses 2 terms taylor expansion
                %the position at next timestep uses 3 terms taylor
                %expansion
                  
                  a = 1;
                %
                  roota = sqrt(a);
                %
                  b = 1.0;
                  %saves current position which becomes previous position
                  %in next timestep
                  previous_pos_1 = pos;
                  factor2 = b*timestep/2.0;
                % Find velocity at the quarter step:
                  v1q = roota*vel ;
                % Find velocity at the half step:
                  v2q = v1q + factor2*force/atomicmass;
                % Find position at the half step:
                  p2q = pos + factor2*v2q;
                % Update Hamiltonian:
                % not relevant in this case 
                % Find position at the next step:
                  pos = p2q + factor2*v2q;
                % Find velocity at the three quarter step:
                  force = getForce();
                  v3q = v2q + factor2*force/atomicmass;
                % Find velocity at the next step:
                  vel = roota*v3q ;
                  %runs this algorithm only once before moving on to the
                  %next algorithm
                  run_once = true;
                  run_count = 1;
            
            end
            
            if run_once == true
                
                %current force
                force_1 = getForce();
                
                %position at next timestep
                if run_count == 1
                    %save the current position for the iteration after the
                    %first
                    previous_pos_2 = pos;
                    pos = 2 * pos - previous_pos_1 + timestep ^ (2) * force_1/atomicmass ;
                
                else
                    a = pos;
                    pos = 2 * pos - previous_pos_2 + timestep ^ (2) * force_1/atomicmass ;
                    previous_pos_2 = a;
                    
                end
                
                %force at next timestep
                force_2 = getForce();
                %velocity at next timestep
                vel = vel + 0.5 * timestep * (force_1 + force_2) / atomicmass;
                %saves position at the next timestep as that 
                run_count = run_count + 1;
                                
            end
        
        elseif intalgo == 2
        %MOVEATOMS updates the atomic positions,velocities and forces
        %   using the BBK integrator. 
          factor1 = gamma*timestep/2.0 ;
          factor2 = timestep/(2.0*atomicmass);
          factor3 = sqrt(gamma*kT/factor2);
          normal1 = normrnd(0.0, 1.0, nAtoms, 3);
          normal2 = normrnd(0.0, 1.0, nAtoms, 3);
        % Find velocity at the half step:
          v2q = (1-factor1)*vel + factor2*(force + factor3*normal1);
        % Find position at the next step:
          pos = pos + timestep*v2q;
        % Find velocity at the next step:
          force = getForce();
          vel = (v2q + factor2*(force + factor3*normal2))/(1+factor1);
        end 
    end


     function relaxAtoms ()
         factor = timestep_factor*timestep/(2.0*atomicmass);
         force = getForce();
         pos = pos + factor*force;
     end

end
