% The MIT License (MIT)
% 
% Copyright (c) 2018 Sebastian Savidan
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software")
% , to deal in the Software without restriction, including without 
% limitation the rights to use, copy, modify, merge, publish, distribute, 
% sublicense, and/or sell copies of the Software, and to permit persons to 
% whom the Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [general_best_position,general_best_score] = pso_simple(n_param,eval_func,lower_bound,upper_bound,n_particle,max_iter)
%PSO_simple
%   pso_simple does an optimization of parameters for a certain evaluation function:
%
% Parameters:
% -n_param: 	number of parameters to optimize
% -eval_func: 	evaluation function to optimize
% -lower_bound: lower boundary of search space, should be a vector of length n_param
% -upper_bound: upper boundary of search space, should be a vector of length n_param
% -n_particle: 	number of particle for PSO algorithm
% -max_iter: 	max number of iteration
%
% Outputs:
% -general_best_position: 	best parameters set found by the optimization for the particular evaluation function
% -general_best_score: 		best score obtained with the evaluation function

%Verifications TODO:
%-nargin
%-upperBound>lowerBound
%-length(lowerBound) = length(upperBound)
%-nParam >=1
%-evalFunc return something corectly nReturned = 1
%-initialize lB and uB and nParticle


%Initialisation of parameters
c1 = 2;
c2 = 2;
inertia_max = 0.98;
inertia_min = 0.4;

tolerance = 10^-10;
flag_tolerance = 1;
iteration = 1;

%Initialisation of particles
particle_position       = lower_bound + (upper_bound-lower_bound).*rand(n_particle,n_param);
particle_velocity       = zeros(n_particle,n_param);
particle_best_position  = particle_position;

particle_score = zeros(1,n_particle);
particle_best_score = zeros(1,n_particle);

for i=1:n_particle
    particle_best_score(i) = eval_func(particle_position(i,:)');
end

[general_best_score,index] = min(particle_best_score);
general_best_position = particle_position(index,:);


% Main loop
while iteration <max_iter && flag_tolerance

% Update inertia
    inertia = inertia_max - (inertia_max - inertia_min)*iteration/max_iter;
    
% Evaluate partcles
    for i=1:n_particle
        particle_score(i) = eval_func(particle_position(i,:));
        if particle_score(i) < particle_best_score(i)
            particle_best_position(i,:) = particle_position(i,:);
            particle_best_score(i) = particle_score(i);
        end
    end
    [g_best_score(iteration),index] = min(particle_best_score);
    general_best_position = particle_best_position(index,:);

% Update position and velocitiy of particles
    for i = 1:n_particle
        particle_velocity(i,:) = inertia*particle_velocity(i,:) + c1*rand*(particle_best_position(i,:) - particle_position(i,:)) + c2*rand*(general_best_position - particle_position(i,:));
        particle_position(i,:) = particle_position(i,:) + particle_velocity(i,:);
        for j = 1:length(lower_bound)
            if particle_position(i,j) < lower_bound(j)
                particle_position(i,j) = lower_bound(j);
            end
            if particle_position(i,j) > upper_bound(j)
                particle_position(i,j) = upper_bound(j);
            end
        end
    end
    
% Look for stop criterion
    if iteration > 10
        if (g_best_score(iteration)-g_best_score(iteration - 10)) < tolerance
            flag_tolerance = 0;
            disp("Pso ended because of relative change")
        end
    end
    iteration = iteration + 1;
end
end

