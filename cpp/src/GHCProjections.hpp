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


#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Cholesky>

typedef std::pair<double,int> mypair;
bool comparator (const mypair& l, const mypair& r);

class GHCProjections {
public:
    GHCProjections();
    void init(unsigned int numTasks, Eigen::VectorXd tasksize, unsigned int DOFsize);
    void setJacobianMatrices(std::vector<Eigen::MatrixXd> jacobian);
    void setInertiaMatrix(Eigen::MatrixXd inertia);
    bool getAllGeneralizedProjectors(std::vector<Eigen::MatrixXd> & Pvector, Eigen::VectorXd & ranks);
    unsigned int getGeneralizedProjector(unsigned int i, Eigen::MatrixXd & PiAis);
    bool setAlphas(Eigen::MatrixXd a);
    Eigen::MatrixXd getAlphas();
    bool setAlphaIJ(unsigned int i, unsigned int j, double a);
    double getAlphaIJ(unsigned int i, unsigned int j);
    void showHierarchy();
   

private:
    void setNumTasks(unsigned int numTasks);
    void setTasksize(Eigen::VectorXd tasksize);
    void setDOFsize(unsigned int DOFsize);
    void getJaugmented(Eigen::MatrixXd & Jaugmented);
    void getTaskIDXaugmented(Eigen::VectorXd & TaskIDXaugmented);
    void getAiVec(unsigned int i, Eigen::VectorXd & Ai);
    void getAiMat(unsigned int i, Eigen::MatrixXd & Ai);
    void computeQRD(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const & Mat,
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & Q,
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & R,
               int & rank,
               Eigen::VectorXi & Pvec);
    bool checkOrthogonality(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const & Mat);
    bool checkOrthogonality(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const & Mat, const double & threshold);

    std::vector<Eigen::MatrixXd> jacobian;
    std::vector<unsigned int> Tasksize;
    Eigen::MatrixXd alphas;
    Eigen::MatrixXd Jt;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd R;
    Eigen::VectorXi origin;
    Eigen::MatrixXd identityDOFsizeDOFsize;
    Eigen::MatrixXd inertiaInv;
    Eigen::MatrixXd L;
    Eigen::VectorXd TaskIDXaugmented;
    Eigen::VectorXd AiVec;
    Eigen::MatrixXd Jaugmented;
    Eigen::MatrixXd JaugmentedSi;
    Eigen::VectorXd AiVecS;
    Eigen::VectorXd AiVecS2;
    Eigen::MatrixXd AiMatS;
    std::vector<std::pair<double,int> > AiSorted;
    std::pair<double,int> mypair;
    int rankQRD;
    unsigned int overallTasksize;
    unsigned int numTasks;
    unsigned int DOFsize;
    double epsilon;
    bool changeRows;
    Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > chol;
};

