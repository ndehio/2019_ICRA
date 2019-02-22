% Author: Niels Dehio
%
%
% Writing code takes time. Polishing it and making it available to others takes longer! 
% If some parts of the code were useful for your research or for a better understanding 
% of the algorithms, please reward the authors by citing the related publications, 
% and consider making your own research available in this way.
%
% @article{Dehio2018c,
% title={{Prioritized Multi-Objective Robot Control}},
% author={Dehio, Niels},
% journal={PhD dissertation, Technical University of Braunschweig, Germany},
% year={2018},
% url={https://publikationsserver.tu-braunschweig.de/receive/dbbs_mods_00066108}
% }
% 
% @article{Silverio2019,
% author={Silv{\'{e}}rio, Jo{\~{a}}o and Calinon, Sylvain and Rozo, Leonel and Caldwell, Darwin G.},
% journal={IEEE Transactions on Robotics},
% title={Learning Task Priorities from Demonstrations},
% year={2019},
% volume={35},
% number={1},
% pages={78-94},
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

classdef WeightedSumOfSOT < PrioritizationScheme
% Description: Algorithm according to the "weighted sum of multiple stack-of-tasks hierarchy candidates" approach

    properties
        numCandidates = 0;
        candidates = [];
        weightingMat = [];
%         regCoeff = 1e-10;
        regCoeff = 0.0;
    end

    methods
        function p = WeightedSumOfSOT(numTasks, tasksize, DOFsize, doVelocityControl)
            p = p@PrioritizationScheme(numTasks, tasksize, DOFsize, doVelocityControl);
            p.setWeightingMatrix(eye(DOFsize,DOFsize));
        end

        function setWeightingMatrix(p, weightingMat)
            assert(size(weightingMat,1) == p.DOFsize);
            assert(size(weightingMat,2) == p.DOFsize);
            p.weightingMat = weightingMat;
        end
        
        function setCandidates(p, candidates)
            p.numCandidates = length(candidates);
            assert(p.numCandidates>1);
            for i=1:p.numCandidates
                assert(length(candidates{i})<=p.numTasks);
                assert(length(candidates{i})>=1);
            end
            p.candidates = candidates;
        end

        function jointCommands = computePrioritizedJointCommands(p)
            %approach: project each task onto the nullspace of all tasks, but shape projection according to priorities
            jointCommands = computePrioritizedJointCommands@PrioritizationScheme(p);
            candidateJointCommands = p.getAllCandidateJointCommands();
            for taskNr=1:p.numTasks
                jointCommands = jointCommands + p.priorities(taskNr) * candidateJointCommands{taskNr};
            end
        end
        
        function candidateJointCommands = getAllCandidateJointCommands(p)
            candidateJointCommands = {};
            for taskNr=1:p.numTasks
                candidateJointCommands{taskNr} = p.getCandidateJointCommands(taskNr);
            end
        end

        function jointCommands = getCandidateJointCommands(p, i)
            assert(i >= 1 && i <= p.numCandidates);

            taskPriorityOrder = p.candidates{i};
            
            jointCommands = zeros(p.DOFsize,1);
            JacAugmented = [];
            for i=1:length(taskPriorityOrder)
                if i>1
                    if (p.doVelocityControl)
                        proj = p.identityDOFsizeDOFsize - p.weightingMat * JacAugmented' * inv(JacAugmented * p.weightingMat * JacAugmented' + p.regCoeff*eye(size(JacAugmented,1))) * JacAugmented;
                    else
                        proj = p.identityDOFsizeDOFsize - JacAugmented' * inv(JacAugmented * p.weightingMat * JacAugmented' + p.regCoeff*eye(size(JacAugmented,1))) * JacAugmented * p.weightingMat;
                    end
                else
                    proj = p.identityDOFsizeDOFsize;
                end
                
                if (p.doVelocityControl)
                    J=p.jacobians{taskPriorityOrder(i)};
                    pinvJ = J' * inv(J * J' + eye(size(J,1))*p.regCoeff); %TODO do we need to add the weightingMat here?
                    jointCommands = jointCommands + proj * pinvJ * p.cartVelocities{taskPriorityOrder(i)};
                else
                    jointCommands = jointCommands + proj * p.jacobians{taskPriorityOrder(i)}' * p.wrenches{taskPriorityOrder(i)};
                end
                JacAugmented = [JacAugmented; p.jacobians{taskPriorityOrder(i)}];
            end
        end
        
        function ok=setPriorities(p, priorities)
            assert(size(priorities,1) == p.numCandidates);
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