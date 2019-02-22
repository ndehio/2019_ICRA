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

typedef std::pair<double,int> mypair;
bool comparator (const mypair& l, const mypair& r) {
    return l.first > r.first;
}

GHCProjections::GHCProjections() {
    this->mypair = std::pair<double,int>();
    this->rankQRD = -1;

    this->epsilon = 1.0e-06; //small value greater zero
    this->changeRows = false;
    this->changeRows = true; //selecting only the most important rows is necessary when overallTasksize > DOFsize
    this->chol = Eigen::LLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >();
}

void GHCProjections::init(unsigned int numTasks, Eigen::VectorXd tasksize, unsigned int DOFsize) {
    this->setNumTasks(numTasks);
    this->setTasksize(tasksize);
    this->setDOFsize(DOFsize);

    Jt = Eigen::MatrixXd::Zero(DOFsize,overallTasksize);
    Q = Eigen::MatrixXd::Zero(DOFsize,overallTasksize);
    R = Eigen::MatrixXd::Zero(overallTasksize,overallTasksize);
    origin = Eigen::VectorXi::Zero(overallTasksize);
    identityDOFsizeDOFsize = Eigen::MatrixXd::Identity(DOFsize,DOFsize);
    inertiaInv = Eigen::MatrixXd::Identity(DOFsize,DOFsize);
    L = Eigen::MatrixXd::Identity(DOFsize,DOFsize);
    TaskIDXaugmented = Eigen::VectorXd::Zero(overallTasksize);
    AiVec = Eigen::VectorXd::Zero(overallTasksize);
    Jaugmented = Eigen::MatrixXd::Zero(overallTasksize,DOFsize);
    JaugmentedSi = Eigen::MatrixXd::Zero(overallTasksize,DOFsize);
    AiVecS = Eigen::VectorXd::Zero(overallTasksize);
    AiVecS2 = Eigen::VectorXd::Zero(overallTasksize);
    AiMatS = Eigen::MatrixXd::Zero(overallTasksize,overallTasksize);
    AiSorted = std::vector<std::pair<double,int> >();
}

void GHCProjections::setJacobianMatrices(std::vector<Eigen::MatrixXd> jacobian) {
    this->jacobian = jacobian;
}

void GHCProjections::setInertiaMatrix(Eigen::MatrixXd inertia) {
    this->inertiaInv = inertia.inverse();
    this->L = chol.compute(this->inertiaInv).matrixL(); //lower triangular matrix
}

bool GHCProjections::getAllGeneralizedProjectors(std::vector<Eigen::MatrixXd> & ProjArray, Eigen::VectorXd & ranks) {
    assert(ProjArray.size() == numTasks);

    for(unsigned int taskNr=0; taskNr<numTasks; taskNr++){
        ranks(taskNr) = this->getGeneralizedProjector(taskNr, ProjArray[taskNr]);
    }
    return true;
}

unsigned int GHCProjections::getGeneralizedProjector(unsigned int i, Eigen::MatrixXd & proj){
    assert(i >= 0 && i < numTasks);
    assert(proj.rows() == DOFsize);
    assert(proj.cols() == DOFsize);

    // ########################################
    //step 1: sort Ai and Jaugmented
    // ########################################
    //augmented vector concatenating the coefficients a_ij which indicate the priorities of all tasks j with respect to task i
    this->getAiVec(i, AiVec);

    //augmented Jacobian concatenating the Jacobian matrices of all tasks
    this->getJaugmented(Jaugmented); //TODO no necessity to compute for each iteration i ...

    //augmented vector concatenating the task-indexes
    this->getTaskIDXaugmented(TaskIDXaugmented); //TODO no necessity to compute for each iteration i ...

    //assign pairs of value and indexes for all a_ij
    AiSorted.clear();
    for(unsigned int idx=0; idx<overallTasksize; idx++){
        mypair.first = AiVec(idx);
        mypair.second = idx;
        AiSorted.push_back(mypair);
    }
    assert(AiSorted.size() == overallTasksize);

    //sort in descending order with saving an index vector
    std::sort(AiSorted.begin(), AiSorted.end(), comparator);

    //sort AIvector and Jaugmented
    //JaugmentedSi is thus constructed so that tasks which should be
    //the least influenced by task i appear in its first rows, while
    //tasks which can be the most influenced by task i appear
    //in its last rows.
    //The values in Ai are sorted accordingly.
    for(unsigned int idx=0; idx<overallTasksize; idx++){
        JaugmentedSi.row(idx) = Jaugmented.row(AiSorted[idx].second);
        AiVecS(idx) = AiVec(AiSorted[idx].second);
    }

    // ########################################
    //step 2: create orthonormal basis (orthonormalization process)
    // ########################################
    Jt = (JaugmentedSi*L).transpose().cast<double>();
    this->computeQRD(Jt, Q, R, rankQRD, origin);

    if(rankQRD <= 0){
        //this case should not happen!
        proj = identityDOFsizeDOFsize;
    }
    else{
        assert((unsigned int)rankQRD <= overallTasksize);

        if((unsigned int)rankQRD > DOFsize){
            rankQRD = DOFsize;
        }

        // ########################################
        //step 3: compute generalized projector PiAis
        // ########################################
        AiVecS2.setZero();
        for(unsigned int ii=0; ii<(unsigned int)rankQRD; ii++){
            assert(origin(ii) >= 0);
            AiVecS2(ii) = AiVecS(origin(ii));
        }
        AiMatS = AiVecS2.asDiagonal();
        proj = L.transpose().inverse() * ( identityDOFsizeDOFsize - Q.leftCols(rankQRD) * AiMatS.topLeftCorner(rankQRD,rankQRD) * Q.leftCols(rankQRD).transpose() ) * L.transpose();
    }
    return rankQRD;
}

void GHCProjections::setNumTasks(unsigned int numTasks){
    assert(numTasks > 0);

    this->numTasks = numTasks;
    this->alphas = Eigen::MatrixXd::Zero(numTasks,numTasks);

    //create standard alphas
    Eigen::MatrixXd tmp = 0.5 * Eigen::MatrixXd::Ones(numTasks,numTasks);
    for(unsigned int i=0; i<numTasks; i++){
        tmp(i,i)=0.0;
    }
    for(unsigned int i=0; i<numTasks; i++){
        for(unsigned int j=0; j<numTasks; j++){
            if (j>i){
                tmp(i,j) = 0.0;
                tmp(j,i) = 1.0;
            }
        }
    }
    this->setAlphas(tmp);
}

bool GHCProjections::setAlphas(Eigen::MatrixXd a){
    // a(i,i) = 0 => task fully activated
    // 0 < a(i,i) < 1 => gradually task (de-)activation
    // a(i,i) = 1 => task fully deactivated (projected in its own nullspace)

    // a(i,j) = 0 => task j has strict lower priority with respect to task i
    // 0 < a(i,j) < 1 => non strict prioritization:
    // a(i,j) = 1 => task j has strict higher priority with respect to task i

    assert(numTasks > 0);
    assert(a.rows() == numTasks);
    assert(a.cols() == numTasks);

    for(unsigned int i=0; i<numTasks; i++){
        for(unsigned int j=0; j<numTasks; j++){
            if(a(i,j) < 0){
                std::cerr << " a( " << i << "," << j << ") = " << a(i,j) << " < 0" << std::endl;
                return false;
            }
            if(a(i,j) > 1){
                std::cerr << " a( " << i << "," << j << ") = " << a(i,j) << " > 1" << std::endl;
                return false;
            }
            if(i!=j){
                if(a(i,j) + a(j,i) != 1.0){
                    std::cerr << " a( " << i << "," << j << ") = " << a(i,j) << " + " << " a( " << j << "," << i << ") = " << a(j,i) << " != 1.0 " << std::endl;
                    //return false;
                }
            }
        }
    }
    alphas = a;
    return true;
}

Eigen::MatrixXd GHCProjections::getAlphas(){
    return alphas;
}

bool GHCProjections::setAlphaIJ(unsigned int i, unsigned int j, double a){
    assert(i >= 0 && i <= numTasks);
    assert(j >= 0 && j <= numTasks);
    assert(a >= 0 && a <= 1.0);
    if(i==j){
        alphas(i,j) = a;
    }
    else{
        alphas(i,j) = a;
        alphas(j,i) = 1.0-a;//TODO makes sense in most scenarios! (but maybe not always?)
    }
    return true;
}

double GHCProjections::getAlphaIJ(unsigned int i, unsigned int j){
    return alphas(i,j);
}

void GHCProjections::setTasksize(Eigen::VectorXd tasksize){
    assert(tasksize.size() == numTasks);

    Tasksize.clear();
    overallTasksize = 0;
    for(unsigned int taskNr=0; taskNr<numTasks; taskNr++){
        Tasksize.push_back(tasksize(taskNr));
        overallTasksize += tasksize(taskNr);
    }
}

void GHCProjections::setDOFsize(unsigned int DOFsize){
    assert(DOFsize > 0);
    assert(Tasksize.size() == numTasks);

    this->DOFsize = DOFsize;
    jacobian.clear();
    for(unsigned int taskNr=0; taskNr<numTasks; taskNr++){
        jacobian.push_back(Eigen::MatrixXd::Zero(Tasksize[taskNr],DOFsize));
    }
}

void GHCProjections::getJaugmented(Eigen::MatrixXd & Jaugmented){
    assert(Jaugmented.rows() == overallTasksize);
    assert(Jaugmented.cols() == DOFsize);

    Jaugmented.setZero();
    unsigned int pos = 0;
    for(unsigned int taskNr=0; taskNr<numTasks; taskNr++){
        Jaugmented.block(pos,0,Tasksize[taskNr],DOFsize) = jacobian[taskNr];
        pos += Tasksize[taskNr];
    }
    assert(pos == overallTasksize);
}

void GHCProjections::getTaskIDXaugmented(Eigen::VectorXd & TaskIDXaugmented){
    assert(TaskIDXaugmented.rows() == overallTasksize);
    assert(TaskIDXaugmented.cols() == 1);

    TaskIDXaugmented.setZero();
    unsigned int pos = 0;
    for(unsigned int taskNr=0; taskNr<numTasks; taskNr++){
        for(unsigned int taskDim=0; taskDim<Tasksize[taskNr]; taskDim++){
            TaskIDXaugmented(pos) = taskNr;
            pos += 1;
        }
    }
    assert(pos == overallTasksize);
}

void GHCProjections::getAiVec(unsigned int i, Eigen::VectorXd & Ai){
    assert(Ai.rows() == overallTasksize);
    assert(Ai.cols() == 1);
    Ai.setZero();
    unsigned int pos = 0;
    for(unsigned int taskNr=0; taskNr<numTasks; taskNr++){
        for(unsigned int taskDim=0; taskDim<Tasksize[taskNr]; taskDim++){
            Ai(pos) = alphas(i,taskNr);
            pos += 1;
        }
    }
    assert(pos == overallTasksize);
}

void GHCProjections::getAiMat(unsigned int i, Eigen::MatrixXd & Ai){
    assert(Ai.rows() == overallTasksize);
    assert(Ai.cols() == overallTasksize);

    Ai.setZero();
    unsigned int pos = 0;
    for(unsigned int taskNr=0; taskNr<numTasks; taskNr++){
        Ai.block(pos,pos,Tasksize[taskNr],Tasksize[taskNr]) = alphas(i,taskNr) * Eigen::MatrixXd::Identity(Tasksize[taskNr],Tasksize[taskNr]);
        pos += Tasksize[taskNr];
    }
    assert(pos == overallTasksize);
}

void GHCProjections::showHierarchy(){
    // a(i,i) = 0 => task fully activated
    // 0 < a(i,i) < 1 => gradually task (de-)activation
    // a(i,i) = 1 => task fully deactivated (projected in its own nullspace)

    // a(i,j) = 0 => task j has strict lower priority with respect to task i
    // 0 < a(i,j) < 1 => non strict prioritization:
    // a(i,j) = 1 => task j has strict higher priority with respect to task i

    assert(numTasks > 0);
    assert(alphas.rows() == numTasks);
    assert(alphas.cols() == numTasks);
    for(unsigned int i=0; i<numTasks; i++){
        for(unsigned int j=0; j<numTasks; j++){
            if(i==j){
                if(alphas(i,j) == 0.0){
                    std::cout << i << " fully activated" << std::endl;
                }
                else if(alphas(i,j) > 0.0 && alphas(i,j) < 1.0){
                    std::cout << i << " gradually de-/activated" << std::endl;
                }
                else if(alphas(i,j) == 1.0){
                    std::cout << i << " fully deactivated" << std::endl;
                }
                else{
                    std::cout << i << " & " << j << " is wrong = " << alphas(i,j) << std::endl;
                }
            }
            else{
                if(alphas(i,j) == 0.0){
                    std::cout << i << " >> " << j << std::endl;
                }
                else if(alphas(i,j) > 0.0 && alphas(i,j) < 0.5){
                    std::cout << i << " >  " << j << std::endl;
                }
                else if(alphas(i,j) == 0.5){
                    std::cout << i << " =  " << j << std::endl;
                }
                else if(alphas(i,j) > 0.5 && alphas(i,j) < 1.0){
                    std::cout << i << " <  " << j << std::endl;
                }
                else if(alphas(i,j) == 1.0){
                    std::cout << i << " << " << j << std::endl;
                }
                else{
                    std::cout << i << " & " << j << " is wrong = " << alphas(i,j) << std::endl;
                }
            }
        }
    }
}

void GHCProjections::computeQRD(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const & Mat,
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & Q,
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> & R,
               int & rank,
               Eigen::VectorXi & Pvec){
    //expect Mat=Jacobian.transpose()

    unsigned int mRows = Mat.rows();
    unsigned int nCols = Mat.cols();
    assert(mRows > 0);
    assert(nCols > 0);
    assert(Q.rows() == mRows);
    assert(Q.cols() == nCols);
    assert(R.rows() == nCols);
    assert(R.cols() == nCols);
    assert(Pvec.rows() == nCols);
    assert(Pvec.cols() == 1);

    Q.setZero();
    R.setZero();
    Pvec.setConstant(-1);

    if(Mat.isZero()){
        rank = 0;
        return;
    }

    unsigned int i = 0;
    unsigned int l;
    for(unsigned int k=0; k<=nCols-1; k++){
        if (changeRows){
            l=i;
        }
        else{
            l=k;
        }

        if (l >= mRows+1){
            break;
        }
        Q.col(l) = Mat.col(k);

        for(int j=0; j<=(int(l)-1); j++){
            R(j,l) = Q.col(j).transpose() * Mat.col(k);
            //R(j,l) = Q.col(j).transpose() * Q.col(l);
            Q.col(l) = Q.col(l) - Q.col(j) * (Q.col(l).transpose() * Q.col(j));
        }
        R(l,l) = Q.col(l).norm();
        if (R(l,l) > epsilon){
            Q.col(l) = Q.col(l) / R(l,l);
            Pvec(l) = k;
            i = i+1;
        }
    }
    rank = i;
}

bool GHCProjections::checkOrthogonality(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const & Mat){
    return this->checkOrthogonality(Mat, 0.00001);
}

bool GHCProjections::checkOrthogonality(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const & Mat, const double & threshold){
    //check if each combination of two task directions are orthogonal (Mat(k) * Mat(j)^T == 0)
    assert(threshold>0 && threshold<1);
    bool ok = true;
    double orthogonality;
    for(unsigned int a=0; a<Mat.cols(); a++){
        for(unsigned int b=a+1; b<Mat.cols(); b++){
            orthogonality = Mat.col(a).transpose() * Mat.col(b);
            if(std::isinf(orthogonality)){
                ok = false;
                std::cout << "orthogonality inf-issue for a=" << a << " and b=" << b << " => Mat(a)*Mat(b)^T = " << orthogonality << std::endl;
            }
            if(orthogonality!=orthogonality){
                ok = false;
                std::cout << "orthogonality nan-issue for a=" << a << " and b=" << b << " => Mat(a)*Mat(b)^T = " << orthogonality << std::endl;
            }
            if(orthogonality > threshold || orthogonality < -threshold){
                ok = false;
                std::cout << "orthogonality issue for a=" << a << " and b=" << b << " => Mat(a)*Mat(b)^T = " << orthogonality << std::endl;
            }
        }
    }
    return ok;
}
