/* 
Author: Niels Dehio


Writing code takes time. Polishing it and making it available to others takes longer! 
If some parts of the code were useful for your research or for a better understanding 
of the algorithms, please reward the authors by citing the related publications, 
and consider making your own research available in this way.

@inproceedings{Dehio2018b,
title={{Continuously Shaping Projections and Operational Space Tasks}},
author={Dehio, Niels and Kubus, Daniel and Steil, Jochen J.},
booktitle={IEEE/RSJ Int. Conf. on Intelligent Robots and Systems},
pages={5995--6002}, 
year={2018}
}

@article{Dehio2018c,
title={{Prioritized Multi-Objective Robot Control}},
author={Dehio, Niels},
journal={PhD dissertation, Technical University of Braunschweig, Germany},
year={2018},
url={https://publikationsserver.tu-braunschweig.de/receive/dbbs_mods_00066108}
}

@inproceedings{Dehio2019,
title={{Dynamically-consistent Generalized Hierarchical Control}},
author={Dehio, Niels and Steil, Jochen J.},
booktitle={IEEE/RSJ Int. Conf. on Robotics and Automation},
year={2019}
}

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License.
If not, see <http://www.gnu.org/licenses/>.
 */
// #########################################################################


#include "GHCProjections.hpp"

int main (int argc, char *argv[]) {
srand((unsigned int) time(0)); //initialize the random number sequence 
int DOFsize = 4;
int numTasks = 3;

//choose dimensionality of each task
Eigen::VectorXd tasksize;
tasksize = Eigen::VectorXd::Zero(numTasks);
//in this example "task0" and "task1" already consume all degress of freedom
tasksize[0] = 2;
tasksize[1] = 2;
tasksize[2] = 3;

//initialize GHC class
GHCProjections GP;
GP = GHCProjections();
GP.init(numTasks, tasksize, DOFsize);

//choose random Jacobians and wrenches for each task
Eigen::MatrixXd jacobian0;
jacobian0 = Eigen::MatrixXd::Random(tasksize[0],DOFsize);
Eigen::MatrixXd jacobian1;
jacobian1 = Eigen::MatrixXd::Random(tasksize[1],DOFsize);
Eigen::MatrixXd jacobian2;
jacobian2 = Eigen::MatrixXd::Random(tasksize[2],DOFsize);

std::vector<Eigen::MatrixXd> allJacobians;
allJacobians.push_back(jacobian0);
allJacobians.push_back(jacobian1);
allJacobians.push_back(jacobian2);
GP.setJacobianMatrices(allJacobians);

Eigen::VectorXd wrench0;
wrench0 = Eigen::VectorXd::Ones(tasksize[0]);
Eigen::VectorXd wrench1;
wrench1 = Eigen::VectorXd::Ones(tasksize[1]);
Eigen::VectorXd wrench2;
wrench2 = Eigen::VectorXd::Ones(tasksize[2]);

//choose random inertia
Eigen::MatrixXd inertia;
inertia = Eigen::MatrixXd::Random(DOFsize,DOFsize);
inertia = inertia * inertia.transpose();
GP.setInertiaMatrix(inertia);

//choose scalar priorities (between one and zero) for each pair of tasks -> there exists 0.5*(numTasks*numTasks+numTasks) pairs!
Eigen::VectorXd prioritiesVector;
prioritiesVector = Eigen::VectorXd::Zero(0.5*(numTasks*numTasks+numTasks));
//exampleA: strict hierachy with "task0" strict more important that "task1" and "task1" strict more important that "task2"
prioritiesVector[0] = 0.0;
prioritiesVector[1] = 0.0;
prioritiesVector[2] = 0.0;
prioritiesVector[3] = 0.0;
prioritiesVector[4] = 0.0;
prioritiesVector[5] = 0.0;
//exampleB: strict hierachy with "task2" strict more important that "task1" and "task1" strict more important that "task0"
prioritiesVector[0] = 0.0;
prioritiesVector[1] = 1.0;
prioritiesVector[2] = 1.0;
prioritiesVector[3] = 0.0;
prioritiesVector[4] = 1.0;
prioritiesVector[5] = 0.0;
//exampleC: all tasks with equal (soft) priority
prioritiesVector[0] = 0.0;
prioritiesVector[1] = 0.5;
prioritiesVector[2] = 0.5;
prioritiesVector[3] = 0.0;
prioritiesVector[4] = 0.5;
prioritiesVector[5] = 0.0;

int counter=0;
for(unsigned int i=0; i<numTasks; i++){
 for(unsigned int j=i; j<numTasks; j++){
  if(prioritiesVector(counter) < 0){
   std::cerr << " a( " << i << "," << j << ") = " << prioritiesVector(counter) << " < 0" << std::endl;
   return false;
  }
  if(prioritiesVector(counter) > 1){
   std::cerr << " a( " << i << "," << j << ") = " << prioritiesVector(counter) << " > 1" << std::endl;
   return false;
  }
  GP.setAlphaIJ(i,j,prioritiesVector(counter));
  counter++;
 }
}
assert(counter==prioritiesVector.size());


//approach: project each task onto the nullspace of all tasks, but shape projection according to priorities
std::vector<Eigen::MatrixXd> allProjections;
for(unsigned int i=0; i<numTasks; i++){
 allProjections.push_back(Eigen::MatrixXd::Zero(DOFsize,DOFsize));
}
Eigen::VectorXd ranks;
ranks = Eigen::VectorXd::Zero(numTasks);
bool ok = GP.getAllGeneralizedProjectors(allProjections, ranks);

//compute final torque commands
Eigen::VectorXd gravityCompensation;
gravityCompensation = Eigen::VectorXd::Random(DOFsize);
Eigen::VectorXd torqueCommand;
torqueCommand = Eigen::VectorXd::Zero(DOFsize);
torqueCommand = gravityCompensation
+ allProjections[0] * jacobian0.transpose() * wrench0
+ allProjections[1] * jacobian1.transpose() * wrench1
+ allProjections[2] * jacobian2.transpose() * wrench2;

//print some results
std::cout << "projection for task0 \n" << allProjections[0] << std::endl;
std::cout << "projection for task1 \n" << allProjections[1] << std::endl;
std::cout << "projection for task2 \n" << allProjections[2] << std::endl;

std::cout << "torques for task0 after being projected \n" << allProjections[0] * jacobian0.transpose() * wrench0<< std::endl;
std::cout << "torques for task1 after being projected \n" << allProjections[1] * jacobian1.transpose() * wrench1<< std::endl;
std::cout << "torques for task2 after being projected \n" << allProjections[2] * jacobian2.transpose() * wrench2<< std::endl;

Eigen::MatrixXd priorityMatrix;
priorityMatrix = Eigen::MatrixXd::Zero(DOFsize,DOFsize);
priorityMatrix = GP.getAlphas();
std::cout << "priorityMatrix \n" << priorityMatrix << std::endl;

//GP.showHierarchy();

std::cout << "done =)" << std::endl;
return 0;
}

