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
% @article{Dehio2018c,
% title={{Prioritized Multi-Objective Robot Control}},
% author={Dehio, Niels},
% journal={PhD dissertation, Technical University of Braunschweig, Germany},
% year={2018},
% url={https://publikationsserver.tu-braunschweig.de/receive/dbbs_mods_00066108}
% }
% 
% @inproceedings{Dehio2019,
% title={{Dynamically-consistent Generalized Hierarchical Control}},
% author={Dehio, Niels and Steil, Jochen J.},
% booktitle={IEEE/RSJ Int. Conf. on Robotics and Automation},
% year={2019}
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

clc;
close all;
clear all;

% A demo program that computes projections for a set of tasks according to the "Generalized Hierarchical Control" (GHC) Approach
DOFsize = 4;
numTasks = 3; %in this example "task0" and "task1" already consume all degress of freedom

%choose dimensionality of each task
tasksize = zeros(numTasks,1);
tasksize(1) = 2;
tasksize(2) = 2;
tasksize(3) = 3;

%initialize prioritization schemes
DG = DynGHC(numTasks, tasksize, DOFsize, false);
WS = WeightedSumOfSOT(numTasks, tasksize, DOFsize, false);
SP = SoftPrioritization(numTasks, tasksize, DOFsize, false);

%choose random Jacobians for each task
allJacobians{1} =rand(tasksize(1),DOFsize);
allJacobians{2} =rand(tasksize(2),DOFsize);
allJacobians{3} =rand(tasksize(3),DOFsize);
DG.setJacobianMatrices(allJacobians);
WS.setJacobianMatrices(allJacobians);
SP.setJacobianMatrices(allJacobians);

%choose random Wrenches for each task
allWrenches{1} = ones(tasksize(1),1);
allWrenches{2} = ones(tasksize(2),1);
allWrenches{3} = ones(tasksize(3),1);
DG.setWrenches(allWrenches);
WS.setWrenches(allWrenches);
SP.setWrenches(allWrenches);

%choose gravity compensation
gravityCompensation = rand(DOFsize,1);
DG.setGravityCompensation(gravityCompensation);
WS.setGravityCompensation(gravityCompensation);
SP.setGravityCompensation(gravityCompensation);

%choose weighting matrix
weightingMat = eye(DOFsize);
DG.setWeightingMatrix(weightingMat);
WS.setWeightingMatrix(weightingMat);

%define candidate SoT hierarchies
candidates{1} = [1,2,3];
candidates{2} = [1,3,2];
candidates{3} = [2,1,3];
candidates{4} = [2,3,1];
candidates{5} = [3,1,2];
candidates{6} = [3,2,1];
candidates{7} = [1,2];
candidates{8} = [1,3];
candidates{9} = [2,1];
candidates{10}= [2,3];
candidates{11}= [3,1];
candidates{12}= [3,2];
candidates{13}= [1];
candidates{14}= [2];
candidates{15}= [3];
WS.setCandidates(candidates);

%choose scalar GHC priorities (between one and zero) for each pair of tasks -> there exists 0.5*(numTasks*numTasks+numTasks) pairs!
prioritiesVector = zeros(0.5*(numTasks*numTasks+numTasks),1);
%exampleA: strict hierachy with "task0" strict more important that "task1" and "task1" strict more important that "task2"
prioritiesVector(1) = 0.0;
prioritiesVector(2) = 0.0;
prioritiesVector(3) = 0.0;
prioritiesVector(4) = 0.0;
prioritiesVector(5) = 0.0;
prioritiesVector(6) = 0.0;
%exampleB: strict hierachy with "task2" strict more important that "task1" and "task1" strict more important that "task0"
% prioritiesVector(1) = 0.0;
% prioritiesVector(2) = 1.0;
% prioritiesVector(3) = 1.0;
% prioritiesVector(4) = 0.0;
% prioritiesVector(5) = 1.0;
% prioritiesVector(6) = 0.0;
%exampleC: all tasks with equal (soft) priority
% prioritiesVector(1) = 0.0;
% prioritiesVector(2) = 0.5;
% prioritiesVector(3) = 0.5;
% prioritiesVector(4) = 0.0;
% prioritiesVector(5) = 0.5;
% prioritiesVector(6) = 0.0;

counter=1;
for i=1:numTasks
 for j=i:numTasks
  if(prioritiesVector(counter) < 0)
   disp([ " a( " num2str( i ) "," num2str( j ) ") = " num2str( prioritiesVector(counter) ) " < 0" ]);
   assert(false);
  end
  if(prioritiesVector(counter) > 1)
   disp([ " a( " num2str( i ) "," num2str( j ) ") = " num2str( prioritiesVector(counter) ) " > 1" ]);
   assert(false);
  end
  DG.setAlphaIJ(i,j,prioritiesVector(counter));
  counter = counter+1;
 end
end
assert(counter-1==length(prioritiesVector));

%choose scalar weighted-SOT priorities 
prioritiesWS = rand(length(candidates),1);
prioritiesSP = rand(numTasks,1);
WS.setPriorities(prioritiesWS);
SP.setPriorities(prioritiesSP);

%compute final torque commands
torqueCommandDG = DG.computePrioritizedJointCommands();
torqueCommandWS = WS.computePrioritizedJointCommands();
torqueCommandSP = SP.computePrioritizedJointCommands();
display(torqueCommandDG);
display(torqueCommandWS);
display(torqueCommandSP);

disp("done!");