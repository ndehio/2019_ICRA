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

classdef PrioritizationScheme < handle
% Description: Generic super class fpr prioritization schemes

    properties
        priorities = [];
        gravityCompensation = [];
        wrenches = {};
        cartVelocities = {};
        jacobians = {};
        Tasksize = [];
        overallTasksize = 0;
        identityDOFsizeDOFsize = [];
        numTasks = 0;
        DOFsize = 0;
        doVelocityControl = NaN;
    end

    methods
        function p = PrioritizationScheme(numTasks, tasksize, DOFsize, doVelocityControl)
            p.setNumTasks(numTasks);
            p.setTasksize(tasksize);
            p.setDOFsize(DOFsize);
            p.doVelocityControl = doVelocityControl;
            p.identityDOFsizeDOFsize = eye(DOFsize,DOFsize);
        end

        function setJacobianMatrices(p, jacobians)
            for taskNr=1:p.numTasks
                assert(size(jacobians{taskNr},1) == p.Tasksize(taskNr));
                assert(size(jacobians{taskNr},2) == p.DOFsize);
            end
            p.jacobians=jacobians;
        end
        
        function ok = setGravityCompensation(p, gravityCompensation)
            if (p.doVelocityControl==true)
                ok = false;
            else
                assert(size(gravityCompensation,1)==p.DOFsize);
                p.gravityCompensation = gravityCompensation;
                ok = true;
            end
        end
        
        function ok = setWrenches(p, wrenches)
            if (p.doVelocityControl==true)
                ok = false;
            else
                for taskNr=1:p.numTasks
                    assert(size(wrenches{taskNr},1)==p.Tasksize(taskNr));
                end
                p.wrenches = wrenches;
                ok = true;
            end
        end
        
        function ok = setCartVelocities(p, cartVelocities)
            if (p.doVelocityControl==false)
                ok = false;
            else
                for taskNr=1:p.numTasks
                    assert(size(cartVelocities{taskNr},1)==p.Tasksize(taskNr));
                end
                p.cartVelocities = cartVelocities;
                ok = true;
            end
        end
        
        function setNumTasks(p, numTasks)
            assert(numTasks > 0);
            p.numTasks = numTasks;
        end

        function setTasksize(p, tasksize)
            assert(size(tasksize,1) == p.numTasks);
            assert(size(tasksize,2) == 1);
            p.Tasksize=tasksize;
            p.overallTasksize = sum(p.Tasksize);
        end

        function setDOFsize(p, DOFsize)
            assert(DOFsize > 0);
            assert(size(p.Tasksize,1) == p.numTasks);
            p.DOFsize = DOFsize;
        end
        
        function jointCommands = computePrioritizedJointCommands(p)
            if (p.doVelocityControl)
                jointCommands = zeros(p.DOFsize,1);
            else
                jointCommands = p.gravityCompensation;
            end
            %TODO: implement specific function in individual subclass
        end
        
        function ok=setPriorities(p, priorities)
            %TODO: implement specific function in individual subclass
        end
    end
end