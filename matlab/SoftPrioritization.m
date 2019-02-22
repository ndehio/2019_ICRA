% Author: Niels Dehio
%
%
% Writing code takes time. Polishing it and making it available to others takes longer! 
% If some parts of the code were useful for your research or for a better understanding 
% of the algorithms, please reward the authors by citing the related publications, 
% and consider making your own research available in this way.
%
% @inproceedings{Dehio2015,
% title={{Multiple Task Optimization with a Mixture of Controllers for Motion Generation}},
% author={Dehio, Niels and Reinhart, R. Felix and Steil, Jochen J.},
% booktitle={IEEE/RSJ Int. Conf. on Intelligent Robots and Systems},
% pages={6416--6421}, 
% year={2015}
% }
% 
% @inproceedings{Dehio2016,
% title={{Continuous Task-Priority Rearrangement during Motion Execution with a Mixture of Torque Controllers}},
% author={Dehio, Niels and Reinhart, R. Felix and Steil, Jochen J.},
% booktitle={IEEE/RAS Int. Conf. on Humanoid Robots},
% pages={264--270}, 
% year={2016}
% }
% 
% @article{Dehio2018c,
% title={{Prioritized Multi-Objective Robot Control}},
% author={Dehio, Niels},
% journal={PhD dissertation, Technical University of Braunschweig, Germany},
% year={2018},
% url={https://publikationsserver.tu-braunschweig.de/receive/dbbs_mods_00066108}
% }
% 
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 3 as
% published by the Free Software Foundation.
% 
% This code is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License.
% If not, see <http://www.gnu.org/licenses/>.
%
% #########################################################################

classdef SoftPrioritization < PrioritizationScheme
% Description: Algorithm according to the "soft weighted superposition of multiple tasks" approach

    properties
%         regCoeff = 1e-10;
        regCoeff = 0;
    end

    methods
        function p = SoftPrioritization(numTasks, tasksize, DOFsize, doVelocityControl)
            p = p@PrioritizationScheme(numTasks, tasksize, DOFsize, doVelocityControl);
        end
        
        function jointCommands = computePrioritizedJointCommands(p)
            %approach: scalar weight for each task 
            jointCommands = computePrioritizedJointCommands@PrioritizationScheme(p);
            for taskNr=1:p.numTasks
                if (p.doVelocityControl)
                    J=p.jacobians{taskNr};
                    pinvJ = J' * inv(J * J' + eye(size(J,1))*p.regCoeff); %TODO do we need to add the weightingMat here?
                    jointCommands = jointCommands + p.priorities(taskNr) * pinvJ * p.cartVelocities{taskNr};
                else
                    jointCommands = jointCommands + p.priorities(taskNr) * p.jacobians{taskNr}' * p.wrenches{taskNr};
                end
            end
        end
        
        function ok=setPriorities(p, priorities)
            assert(size(priorities,1) == p.numTasks);
            ok=true;
            for i=1:p.numTasks
                if(1.0 < priorities(i) || 0.0 > priorities(i))
                    ok = false;
                    break;
                end
            end
            if ok
                p.priorities = priorities;
            end
        end
    end
end