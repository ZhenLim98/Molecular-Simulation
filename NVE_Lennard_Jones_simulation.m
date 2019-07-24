
function begin = NVE_Lennard_Jones_simulation(position_matrix,momentum_matrix, mass_vector, tstep, epsilon, sigma)
% implement energy?
% implement periodic BCs (how? - use pseudocode? - see plan for this week)
% print out code and annotate it and try to find bugs; if not they will
% persist. then execute the rest of the code. go steady and slow - too fast
% will mean it will take too long to fix it. Especially that thing
%initial conditions
position_matrix = [1,0,0 ; 0,0,0 ; -1,0,0]' ; % a,b,c creates row vector so need to transpose it
momentum_matrix = [0,0,0 ; 0,0,0 ; 0,0,0]' ;
mass_vector = [1,2,3];
tstep = 1;
epsilon = 1.65 * 10^(-21) ;
sigma = 3.4 * 10^(-10) ;
t = 0;
    
    for time = 1 :tstep: 100 
    t = t + tstep;
    position_matrix_saved = position_matrix;
    momentum_matrix_saved = momentum_matrix ;
    % saving the updated positions
    a = update_position(position_matrix, tstep, momentum_matrix, mass_vector, epsilon, sigma);
    position_matrix = a;
    %saving the updated momenta
    b = update_momentum( momentum_matrix_saved, position_matrix_saved, tstep, mass_vector, epsilon, sigma );
    %( momentum_matrix, saved_position_matrix, tstep, mass_vector, epsilon, sigma)
    momentum_matrix = b;
    disp(b)
    string = 'yolo';
    %x = hi(string);
    x = total_energy(position_matrix, momentum_matrix, epsilon, sigma, mass_vector);
    TotalE = x;
    
    end

   % plot(t , TotalE)


%%%%%%%%%%%%%%%%%%%%%%%%nested functions%%%%%%%%%%%%%%%%%%%%%%
%     function [xd] = hi(string)
%         xd = string;
%         disp(xd)
%         
%     end

    function [position_matrix] = update_position(position_matrix, tstep, momentum_matrix, mass_vector, epsilon, sigma)
    %this function updates the atomic position (Verlet algorithm) by a SINGLE timestep delta
    %the _initial subscript refers to the parameters at time t, updated is for
    % t + tstep ()#
    
    [~, columns] = size (position_matrix);
    for i = 1:columns 
        
        % column vector is position vector
        r_initial = position_matrix( : , i );
        
        %mass is mass vector
        mass_i = mass_vector( i );
        %column vector is momentum_matrix
        p_initial = momentum_matrix(: , i);
        % force becomes not a number
        r_updated = r_initial + tstep * p_initial + ( tstep ^ 2 )/(2 * mass_i )* ( force_n(position_matrix,i,  epsilon, sigma) );

        position_matrix(:, i) = ( r_updated );
        % stores the new position vector in the position matrix as a column
        % vector
    end    

    end

    function [ momentum_matrix ] = update_momentum( momentum_matrix, saved_position_matrix, tstep, mass_vector, epsilon, sigma)
    % this function updates momentum for every single  with input VECTORS
    [ ~ , columns ] = size(momentum_matrix);
    for i = 1:columns
        p_initial = momentum_matrix ( : , i) ;
        updated_position_matrix = update_position(saved_position_matrix, tstep, momentum_matrix, mass_vector, epsilon, sigma);

        p_updated = p_initial + tstep / 2 * ( force_n(saved_position_matrix, i, epsilon, sigma) + force_n(updated_position_matrix, i, epsilon, sigma));    
        momentum_matrix(: , i) = p_updated;
        
    end
    end


    function [ F_total ] = force_n(position_matrix, n, epsilon, sigma)
    %calculates the force on a SINGLE (nth) particle due to all others

    [~ , columns] = size(position_matrix) ;
    r_i = position_matrix(:,n);
    F_total = ( [0,0,0] )';

    for j = 1:columns  
        if ( j ~= n)
        % for all the particles in the position matrix except itself
        r_j = position_matrix( : , j ) ;
        sep_ij = r_i - r_j;
        F_addition = epsilon * ( 12 * ( sigma / norm ( sep_ij ) ^11 * (-sigma) * sep_ij / ( norm(sep_ij)^3 )) - 6 * ((sigma/norm(sep_ij))^5 *(- sigma ) * sep_ij / ( norm (sep_ij )^3 ) ) ) ;

        F_total = F_total + F_addition ;
        end
    end


    end

    function [totalenergy] = total_energy(position_matrix, momentum_matrix, epsilon, sigma, mass_vector)% check which variables i need
        %see langevin MD to see if there's anything in there i can use e.g the way he used potentials to find forces 
        
        [~ , columns] = size(position_matrix) ;
        
        KE = 0;
        PE = 0;
        
       
        for s = 1:columns
            %calculating the potential between neighbouring molecules and
            %the j'th molecule
            for v= 1:columns
                if ( v ~= s)
                    r_1 = position_matrix(:,v);
                    r_2 = position_matrix(:,s);
                    sep_12 = (r_2 - r_1);
                    norm_sep_12 = norm(sep_12);
                   
                    %muliply each potential term by 0.5 as each molecular
                    %interaction is counted twice.
                    
                    PE_single = 0.5 * epsilon * ( ( sigma / norm_sep_12 ) ^(12) - ( sigma / norm_sep_12 )^(6) );
                    %calculate PE using Lennard_jones due to the j'th molecule 

                    PE = PE + PE_single;
                end

            end
            
            mass_s = mass_vector(s);
            %calculating the kinetic energy
            p = momentum_matrix (: , s);
            KE_single = norm(p)^2 / (2 * mass_s) ;
            KE = KE + KE_single ; 
            
        end
        totalenergy = KE + PE;
    end
    
end


