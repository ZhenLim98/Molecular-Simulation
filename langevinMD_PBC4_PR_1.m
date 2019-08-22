function [ nSteps, t, E_t, T_t, R2_t, D_t Fmax_t Frms_t] = langevinMD_PBC4_PR_1( nAtoms, gamma, targetTemp,P_target, timestep, equilTime, runTime, nocells1D, h,V,mdmode, startconf)
%LANGEVINMD Perform Langevin Molecular Dynamics on a strand of polymer
%
%   Input variables:
%   nAtoms          Number of atoms
%   gamma           Damping rate in Langevin dynamics in fs^-1
%   targetTemp      Temperature of atoms in kelvin
%   timestep        Molecular dynamics time step in fs
%   equilTime       Time for equilibrating the system in fs
%   runTime         Time over which statistics are taken in fs
%   nocells1D       Number of additional cells on each side of the unit
%                   cell e.g nocells1D = 1 means 7 cells in total
%                   (Note: HAS TO BE EVEN)
%   h               Matrix of basis vectors (as column vectors) that span the unit cell 
%                   (Note: REMEMBER IT'S IN ANGSTROM)
%   V               In P-R (at least) its =  h_dot h^-1 where h_dot is the
%                   matrix of unit cell vector velocities (in columns) and
%                   h is given above
%   P_target        Target pressure (kPa)       
%                  
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
%   nocells
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
intalgo = 3;   %Specifies the integration algorithm for the equations of motion

% The following potential parameters are taken from C. Liu and M.
% Muthukumar, "Langevin dynamics simulations of early-stage polymer
% nucleation and crystallization", J. Chem. Phys., 109, 2536 (1998)
sigma = 4.5;          % Angstroms
epsilon = 0.00485678; % eV; 0.112 kcal/mol
bondLength = 10;    % Angstroms
Kr = 15.18; % eV/A^2; 350 kcal/mol A^2

run_once = false; 
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
W = 3 * atomicmass * nAtoms / 4 * pi ^2; %Calculated using Parrinello-Rahman
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
pos = zeros(3, nAtoms );
vel = zeros(3, nAtoms );
%
switch startconf
    case 'restart'
%
        fprintf ('\nReading input configuration from file.\n');
        inputFid = fopen('restart_file', 'r');
        nAtoms=fscanf(inputFid, '%i', 1);
        size1 = [ 3 nAtoms  ];
        pos=fscanf(inputFid, '%g', size1);        
        vel=fscanf(inputFid, '%g', size1);
        fclose(inputFid);
%
    otherwise
     fprintf ('\nConstructing new input configuration.\n');
     pos(3,:) = linspace(0, (nAtoms-1)*bondLength, nAtoms);
% Initialise the atomic velocities according to the Boltzmann distribution.
% The temperature has been increased above the target value to allow for
% the initial relaxation.
      vel = normrnd(0.0, sqrt(1.5*kT/atomicmass), 3, nAtoms );
end
fprintf ('\nConfiguration contains %i atoms.\n', nAtoms);
%
% Translate to the origin
%
pos_transpose = pos';
cmpos= mean(pos_transpose);
pos(1,:) = pos(1,:) - cmpos(1);
pos(2,:) = pos(2,:) - cmpos(2);
pos(3,:) = pos(3,:) - cmpos(3);
%
% Calculate initial potential energy
%
pot1_i = getBondingEnergy()/nAtoms;
pot2_i = getNonBondingEnergy()/nAtoms;
pot_i = pot1_i + pot2_i;
%
% Remove centre of mass velocity
%
vel_transpose = vel';
cmvel = mean(vel_transpose);
vel(1,:) = vel(1,:) - cmvel(1);
vel(2,:) = vel(2,:) - cmvel(2);
vel(3,:) = vel(3,:) - cmvel(3);
%
% Calculate initial temperature and kinetic energy
%
kinstart = getKinetic();
disp(kinstart)
temperature = getTemperature();  
fprintf ('\nStarting configuration:');
fprintf ('\n-----------------------');
fprintf ('\nTemperature = %f K \n', temperature);
%
% Initialize the atomic forces
%
force = zeros(3,nAtoms);
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
   vel = zeros(3, nAtoms );
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
    R2_t(i) = (pos(nAtoms,1)-pos(1,1))'*(pos(nAtoms,1)-pos(1,1))'...
            + (pos(nAtoms,2)-pos(1,2))'*(pos(nAtoms,2)-pos(1,2))' ...
            + (pos(nAtoms,3)-pos(1,3))'*(pos(nAtoms,3)-pos(1,3))';
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
% Computes the kinetic energy from the velocities and masses for ALL the
% particles in a single. what about moving basis vectors/KE of cell(PR)?

  kinetic = 0.5*atomicmass*(dot(vel(1,:), vel(1,:)) ...
  + dot(vel(2,:), vel(2,:)) + dot(vel(3,:), vel(3,:)));

    end
%
    function potential1 = getBondingEnergy()  
  % Compute the bonding contribution to the potential energy for the unit cell which
  % is simply a harmonic bond-stretching term between neighbouring molecules
  % only within the unit cell.
    potential1 = 0.0;
                for i1=1 : nAtoms -1
                    j1=i1+1;
                    r = sqrt((pos(1,i1)-pos(1,j1))^2 ...
                        + (pos(2,i1)-pos(2,j1))^2 ...
                        + (pos(3,i1)-pos(3,j1))^2);
                    potential1 = potential1 + Kr * (r-bondLength)*(r-bondLength);
                end
       
       
%       % additional potential terms for the start and the end of the chain
%        r_ends = sqrt((pos(1,1)-pos(nAtoms,1))^2 ...
%             + (pos(1,2)-pos(nAtoms,2))^2 ...
%             + (pos(1,3)-pos(nAtoms,3))^2);
%        potential1 = potential1 +  2 * Kr *(r_ends - bondLength)^2 ;
    end
%
    function potential2 = getNonBondingEnergy ()
        %   Computes the non-bonding Lennard-Jones energy term
        %   for the unit cell surrounded by extra nocells1D on each side in
        %   each dimension
        % molecule interacts with particles and its images in all other unit cells 
        % but not itself

       potential2 = 0.0;
       for p = -nocells1D:nocells1D
           for q = -nocells1D:nocells1D
               for r = -nocells1D:nocells1D
                    for l=1:nAtoms
                        for i1=1:nAtoms
                            if l == i1 && p == 0 && q == 0 && r == 0
                                %non self-interacting 
                            else
                                component_1 = p * h(1,1) + q * h(2,1) + r * h(3,1);
                                component_2 = p * h(1,2) + q * h(2,2) + r * h(3,2);
                                component_3 = p * h(1,3) + q * h(2,3) + r * h(3,3);
                                r2 = (pos(1,l) - pos(1,i1) -  component_1 )^2 ...
                                + (pos(2,l) - pos(2,i1)- component_2 )^2 ...
                                + (pos(3,l) - pos(3,i1)- component_3 )^2 ;
                                x2 = sigma*sigma/r2;
                                x6 = x2*x2*x2; 
                                x12 = x6*x6;
                                potential2 = potential2 + 0.5 * epsilon * (x12-2.0*x6);
                            end
                        end
                    end
               end
            end
        end
    end


    function temperature = getTemperature ()
        %GETTEMPERATURE computes the temperature
        %calculate for one unit cell only as temperature is intensive 
        vel_transpose = vel';
        vcm = mean(vel_transpose);
        temperature = atomicmass*(dot(vel(1,:)-vcm(1), ....
            vel(1,:)-vcm(1)) + dot(vel(2,:)-vcm(2), ...
            vel(2,:)-vcm(2)) + dot(vel(3,:)-vcm(3), ...
            vel(3,:)-vcm(3)))/(3.0*kB*(nAtoms-1));
    end

    function F = getForce ()
        %GETFORCE computes the forces on the atoms
        %   The forces are computed from the gradients of the potential
        %   energy. The potential energy consists of two components: a
        %   bonding term due to bond stretching, and a non-bonding
        %   Lennard-Jones term between remote atoms (third neighbour or
        %   further). Stochastic forces included, but NOT the friction
        %   forces
        %should the molecules get to the edge of the unit cell,
        %we let them go on (but compute interactions with the nearest images)

        % Compute the bonding contribution to the atomic forces
        %only between adjacent molecules
        F = zeros(3, nAtoms);
        for i2=1:nAtoms-1
            j2=i2+1;
            dPos(1) = pos(1,j2)-pos(1,i2);
            dPos(2) = pos(2,j2)-pos(2,i2);
            dPos(3) = pos(3,j2)-pos(3,i2);


            r = sqrt(dPos(1)*dPos(1) + dPos(2)*dPos(2) + dPos(3)*dPos(3));
            fac = 2.0*Kr*(1.0-bondLength/r);
            F(1,i2) = F(1,i2) + fac*dPos(1);
            F(2,i2) = F(2,i2) + fac*dPos(2);
            F(3,i2) = F(3,i2) + fac*dPos(3);
            F(1,j2) = F(1,j2) - fac*dPos(1);
            F(2,j2) = F(2,j2) - fac*dPos(2);
            F(3,j2) = F(3,j2) - fac*dPos(3);
% %             % forces on first and last molecule outside the unit cell
% %             if i2 == 1
% %                 F(1, 1) = F(1,1) - fac *( pos(1, 1) - pos(nAtoms,1));
% %                 F(1, 2) = F(1,2) - fac*( pos(1, 2) - pos(nAtoms,2));
% %                 F(1, 3) = F(1,3) - fac*( pos(1, 3) - pos(nAtoms,3));
% % 
% %             elseif i2 == nAtoms
% % 
% %                 F(nAtoms, 1) = F(nAtoms, 1) + fac*( pos(1, 1) - pos(nAtoms,1) );
% %                 F(nAtoms, 2) = F(nAtoms, 2) + fac*(pos(1, 2) - pos(nAtoms,2));
% %                 F(nAtoms, 3) = F(nAtoms, 3) + fac*(pos(1, 3) - pos(nAtoms,3));
% % 
%          end

       end

        % Add the non-bonding contribution to the atomic forces with
        % each image in every other unit cell.
        fac1 = -12.0*epsilon/(sigma*sigma);
        for i2=1:nAtoms
            for p = -nocells1D:nocells1D
               for q = -nocells1D:nocells1D
                   for r = -nocells1D:nocells1D
                        for j2=1:nAtoms
                            if  i2 == j2 && p == 0 && q == 0 && r == 0
                                %no self interaction
                            else
                                %interaction with all other images in all
                                %other cells
                                component_1 = p * h(1,1) + q * h(2,1) + r * h(3,1);
                                component_2 = p * h(1,2) + q * h(2,2) + r * h(3,2);
                                component_3 = p * h(1,3) + q * h(2,3) + r * h(3,3);
                                dPos(1) = pos(1,j2)- pos(1,i2) + component_1;
                                dPos(2) = pos(2,j2)- pos(2,i2) + component_2;
                                dPos(3) = pos(3,j2)- pos(3,i2) + component_3;
                                r2 = dPos(1)^2 + dPos(2)^2 + dPos(3)^2 ;
                                x2 = sigma*sigma/r2;
                                x6 = x2*x2*x2;
                                x12 = x6*x6;
                                fac = fac1*x2*(x12-x6);
                                F(1,i2) = F(1,i2) + fac * dPos(1);
                                F(2,i2) = F(2,i2) + fac * dPos(2);
                                F(3,i2) = F(3,i2) + fac * dPos(3);
                                
                            end
                        end
                    end
               end
           end
        end


    end
%
    function moveAtoms ()
        %MOVEATOMS updates the atomic positions,velocities and forces
        %   using the Verlet integrator a SINGLE timestep delta T.
        %No friction. Possibly map positions back if too far out?

        % Update positions
        if intalgo == 1
           disp(pos)
            if run_once == false
                %velocity at next timestep uses 2 terms in taylor expansion
                %position at next timestep uses 3 terms in taylor expansion
                  %saves current position which becomes previous position
                  previous_pos_1 = pos;
                  factor2 = timestep/2.0;
                % Find velocity at the quarter step:
                  v1q = vel ;
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
                  vel = v3q ;
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
                    %save the current position (becomes next previous position)
                    previous_pos_2 = pos;
                    %position at next timestep
                    pos = 2 * pos - previous_pos_1 + timestep ^ (2) * force_1/atomicmass ;

                else
                    % saves current position for next loop
                    a = pos;
                    %position at the next timestep
                    pos = 2 * pos - previous_pos_2 + timestep ^ (2) * force_1/atomicmass ;
                    % position before updating becomes previous position
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
          normal1 = normrnd(0.0, 1.0, 3, nAtoms);
          normal2 = normrnd(0.0, 1.0, 3, nAtoms);
        % Find velocity at the half step:
          v2q = (1-factor1)*vel + factor2*(force + factor3*normal1);
        % Find position at the next step:
          pos = pos + timestep*v2q;
        % Find velocity at the next step:
          force = getForce();
          vel = (v2q + factor2*(force + factor3*normal2))/(1+factor1);
         
        
        elseif intalgo == 3
            %updates the atomic positions, velocities, forces using NPH
            %integrator
%             halfstep = 1;
%             quarterstep = 0;

            %1.
            Molec_velocities_qq ()
            %CALL A_v_
            %2.

            halfstep = 0;
            quarterstep = 1;
            V_update(halfstep,quarterstep)
            %+q
            %3.
            Molec_positions_qqqq ()
            h_qqqq() 
            getForce()
            getStressTensor() %- updates S
            %A_r_qqqq using (H_qqqq, S_qqqq and Ptarget)

            %4.

            V_update(halfstep, quarterstep) 
            %5.
            Molec_velocities_qq()
            %A_v()
            %6.

            halfstep = 1;
            quarterstep = 0;
            V_update (halfstep, quarterstep)

        end
         
    end

 %%NPH integrator functions below
 %what is N_f?
 %each particle occup
 
    function Molec_velocities_qq()
      %updates velocity at the half step, using equation derived (in book). 
      %forces and velocities are used as row vectors in matrix, h is used
      %as a column matrix
%       acc = getForce()/atomicmass;
%       %constructing delta and delta_tilda
%       %transpose to undergo transformations on column vectors
%       vel = vel';
%       % (Eq. 50) in wrapper
%       [c,~] = size(V);
%       C = V + 1/N_f * trace(V) * eye(c);
%       [a, ~] = size(acc');
%       factor1 = expm(-timestep * C);
%       factor2 = C \(eye(a) - expm(-timestep * C));
%       vel = factor1 * vel +  factor2 * acc;
%       
%     %transpose back after linear transformations
%       vel = vel'; 
    I = eye(3);
    force = getForce();
    %using V, and h at the current time, we calculate these other matrices
    %from it
    h_dot = V * h;
    G = h' * h;
    G_dot = h_dot' * h + h' * h_dot;
    s = h \ pos;
    s_dot = h \ vel;
    S = getStressTensor();
    factor1 = 1/W * det(h) * (S - P_target * I) * ( h' \  s); 
    factor2 = 2 * h_dot * s_dot ;
    factor3 = 1/atomicmass * (h \ force) - G \  G_dot * s_dot;
    M = factor1 + factor2 + factor3; 
    trace_v = trace(V);
    %propagating using Louvillian expression.
    vel = vel + timestep * M + (timestep)^2 * trace_v * M + (timestep)^3 * 2 * 2 / factorial(3) * (trace_v ^ 2) * M  ;
    end

    function V_update(halfstep,quarterstep)
       % Eq 55. 
        %constructing A_r
%        [s,;] = size(S);
%         A_r = det(H)*( W \ (S - eye(s) * P_target)) ;
%         
%     if halfstep == false && quarterstep == true:
%         %quarter timestep
%         V = V + timestep/4 * (A_r + A_v);
%     elseif halfstep ==true && quarterstep ==false:
%         %half timestep
%         V = V + timestep/2 * (A_r + A_v);
%     end
    S = getStressTensor();
    identity = eye(3);
    omega = det(h);
    A = ( (S - P_target * identity)/ h' ) * omega/W;
    if halfstep == 0 && quarterstep == 1
        factor1 = timestep/4 * A  / h + identity;
    elseif halfstep == 1 && quarterstep == 0
        factor1 = timestep/2 * A  / h + identity;
    end
    
    V = V + factor1;    

    end

    function Molec_positions_qqqq ()
    %see above 48
    %updates molecular positions
    %constructing D and D_tilda
    d_factor = eig(V);
    d_factor = d_factor';
    d_factor_2 = d_factor;
    for d = [1, length(d_factor)]
    d_factor(d) = exp(d_factor(d) * timestep);
    factor_a = timestep * d_factor_2(d)/2;
    sinh_factor = sinh(factor_a)/factor_a;
    d_factor_2(d) = exp(factor_a) * sinh_factor;
    end
    D = diag(d_factor)  ;
    D_tilda = diag(d_factor_2);
    
    %constructing the matrix of normalised eigenvectors
    
    [~,eig_matrix] = eig(V);
    [~,no_of_eigenvectors] = size(eig_matrix);
    for d1 = [1, no_of_eigenvectors]
         eig_matrix(:,d1) = eig_matrix(:,d1)/norm(eig_matrix(:,d1));        
    end
    O = eig_matrix;
    
    factor1 = O' * D * O;
    factor2 = O' * D_tilda * O;
     
    pos = factor1 * pos + factor2 * vel  ;
    end

    function h_qqqq ()
        %Eq 54
        %constructing D and D_tilda
        d_factor = eig(V);
        d_factor = d_factor';
        for d = [1, length(d_factor)]
        d_factor(d) = exp(d_factor(d) * timestep);

        end
        D = diag(d_factor)  ;
        %constructing the matrix of normalised eigenvectors

        [~,eig_matrix] = eig(V);
        [~,no_of_eigenvectors] = size(eig_matrix);
        for d1 = [1, no_of_eigenvectors]
             eig_matrix(:,d1) = eig_matrix(:,d1)/norm(eig_matrix(:,d1));        
        end
        O = eig_matrix;

        h = O' * D * O * h;
    end
 
%function getforce not updated

    function S = getStressTensor()
        %using the StressTensor from Parrinello-Rahman, Eq. 7
        %stress tensor is 'S', scaled size is 's'
    [~ ,b] = size(vel); %size is being used elsewhere  
    S = zeros(3);
    s =  h \ pos;
    for i1 = [1,b]
       factor1 = atomicmass * vel(:,b) .* vel(:,b);
       factor2 = ( force * s' ) * h';
       S = S + 1/det(h) * ( factor1 + factor2 ) ;
    end
    
    end
   
    function relaxAtoms ()
         factor = timestep_factor*timestep/(2.0*atomicmass);
         force = getForce();
         pos = pos + factor*force;
      end

 
 
end