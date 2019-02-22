% Author: Niels Dehio
%
%
% Writing code takes time. Polishing it and making it available to others takes longer! 
% If some parts of the code were useful for your research or for a better understanding 
% of the algorithms, please reward the authors by citing the related publications, 
% and consider making your own research available in this way.
%
% @inproceedings{Dehio2018b,
% title={{Continuously Shaping Projections and Operational Space Tasks}},
% author={Dehio, Niels and Kubus, Daniel and Steil, Jochen J.},
% booktitle={IEEE/RSJ Int. Conf. on Intelligent Robots and Systems},
% pages={5995--6002}, 
% year={2018}
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

classdef DynGHC < PrioritizationScheme
% Description: Algorithm according to the "Dynamically-consistent Generalized Hierarchical Control" (DynGHC) approach

    properties
        Jaugmented = [];
        alphas = [];
        weightingMat = [];
%         regCoeff = 1e-10;
        regCoeff = 0.0;
        %epsilon = 1.0e-06; 
        epsilon = realmin(); %small value greater zero
        %changeRows = false;
        changeRows = true; %selecting only the most important rows is necessary when overallTasksize > DOFsize
    end

    methods
        function p = DynGHC(numTasks, tasksize, DOFsize, doVelocityControl)
            p = p@PrioritizationScheme(numTasks, tasksize, DOFsize, doVelocityControl);
            p.setWeightingMatrix(eye(DOFsize,DOFsize));
        end

        function setJacobianMatrices(p, jacobians)
            setJacobianMatrices@PrioritizationScheme(p, jacobians);

            p.Jaugmented = zeros(p.overallTasksize, p.DOFsize);
            pos = 1;
            for taskNr=1:p.numTasks
                p.Jaugmented(pos:pos-1+p.Tasksize(taskNr),:) = jacobians{taskNr};
                pos = pos + p.Tasksize(taskNr);
            end
            assert(pos-1 == p.overallTasksize);
            assert(size(p.Jaugmented,1) == p.overallTasksize);
            assert(size(p.Jaugmented,2) == p.DOFsize);
        end

        function setWeightingMatrix(p, weightingMat)
            assert(size(weightingMat,1) == p.DOFsize);
            assert(size(weightingMat,2) == p.DOFsize);
            p.weightingMat = weightingMat;
        end
        
        function jointCommands = computePrioritizedJointCommands(p)
            %approach: project each task onto the nullspace of all tasks, but shape projection according to priorities
            jointCommands = computePrioritizedJointCommands@PrioritizationScheme(p);
            [allProjections, ~] = p.getAllGeneralizedProjectors();
            for taskNr=1:p.numTasks
                if (p.doVelocityControl)
                    J=p.jacobians{taskNr}; 
                    pinvJ = J' * inv(J * J' + eye(size(J,1))*p.regCoeff); 
                    jointCommands = jointCommands + allProjections{taskNr} * pinvJ * p.cartVelocities{taskNr}; 
                else
                    jointCommands = jointCommands + allProjections{taskNr} * p.jacobians{taskNr}' * p.wrenches{taskNr};
                end
            end
        end
        
        function [ProjArray, ranks] = getAllGeneralizedProjectors(p)
            ProjArray = {};
            ranks = zeros(p.numTasks,1);
            for taskNr=1:p.numTasks
                [ProjArray{taskNr}, ranks(taskNr)] = p.getGeneralizedProjector(taskNr);
            end
        end

        function [proj, rankQRD] = getGeneralizedProjector(p, i)
            assert(i >= 1 && i <= p.numTasks);

            % ########################################
            %step 1: sort Ai and Jaugmented
            % ########################################
            %augmented vector concatenating the coefficients a_ij which indicate the priorities of all tasks j with respect to task i
            AiVec = p.getAiVec(i);

            %sort in descending order with saving an index vector
            %sort AIvector and Jaugmented
            [AiVecS,idx]=sort(AiVec, 'descend');

            %JaugmentedSi is thus constructed so that tasks which should be
            %the least influenced by task i appear in its first rows, while
            %tasks which can be the most influenced by task i appear
            %in its last rows.
            %The values in Ai are sorted accordingly.
            for ii=1:p.overallTasksize
                JaugmentedSi(ii,:) = p.Jaugmented(idx(ii),:);
            end

            % ########################################
            %step 2: create orthonormal basis (orthonormalization process)
            % ########################################
            %extension for dynamic consistency:
            L = chol(p.weightingMat, 'lower'); % weightingMat=L*L'
            Jt = (JaugmentedSi*L)';
            [Q, R, rankQRD, origin] = p.computeQRD(Jt);

            assert(rankQRD > 0);
            assert(rankQRD <= p.overallTasksize);

            if(rankQRD > p.DOFsize)
                rankQRD = p.DOFsize;
            end
            %p.checkOrthogonality(Q(:,1:rankQRD));

            % ########################################
            %step 3: compute generalized projector PiAis
            % ########################################
            AiVecS2 = zeros(size(AiVecS));
            for ii=1:rankQRD
                assert(origin(ii) >= 0);
                AiVecS2(ii) = AiVecS(origin(ii));
            end
            if (p.doVelocityControl)
                proj = p.identityDOFsizeDOFsize - L * Q(:,1:rankQRD) * diag(AiVecS2(1:rankQRD)) * Q(:,1:rankQRD)' * inv(L);
            else
                proj = p.identityDOFsizeDOFsize - inv(L') * Q(:,1:rankQRD) * diag(AiVecS2(1:rankQRD)) * Q(:,1:rankQRD)' * L';
            end
        end

        function setNumTasks(p, numTasks)
            assert(numTasks > 0);
            p.numTasks = numTasks;

            %create standard alphas
            tmp = 0.5 * ones(p.numTasks,p.numTasks);
            for i=1:p.numTasks
                tmp(i,i)=0.0;
            end
            for i=1:p.numTasks
                for j=1:p.numTasks
                    if (j>i)
                        tmp(i,j) = 0.0;
                        tmp(j,i) = 1.0;
                    end
                end
            end
            p.setAlphasMat(tmp);
        end

        function ok = setAlphasMat(p, a)
            % a(i,i) = 0 => task fully activated
            % 0 < a(i,i) < 1 => gradually task (de-)activation
            % a(i,i) = 1 => task fully deactivated (projected in its own nullspace)

            % a(i,j) = 0 => task j has strict lower priority with respect to task i
            % 0 < a(i,j) < 1 => non strict prioritization:
            % a(i,j) = 1 => task j has strict higher priority with respect to task i
            assert(p.numTasks > 0);
            assert(size(a,1) == p.numTasks);
            assert(size(a,2) == p.numTasks);
            ok = true;
            for i=1:p.numTasks
                for j=1:p.numTasks
                    if(a(i,j) < 0)
                        %disp([ " a( " num2str( i ) "," num2str( j ) ") = " num2str( a(i,j) ) " < 0" ]);
                        ok = false;
                        break;
                    end
                    if(a(i,j) > 1)
                        %disp([ " a( " num2str( i ) "," num2str( j ) ") = " num2str( a(i,j) ) " > 1" ]);
                        ok = false;
                        break;
                    end
%                     if(i~=j)
%                         if(a(i,j) + a(j,i) ~= 1.0)
%                             disp([ " a( " num2str( i ) "," num2str( j ) ") = " num2str( a(i,j) ) " +  a( " num2str( j ) "," num2str( i ) ") = " num2str( a(j,i) ) " != 1.0 " ]);
%                             %return false; %TODO
%                         end
%                     end
                end
            end
            if ok
                p.alphas = a;
            end
        end

        function a = getAlphasMat(p)
            a = p.alphas;
        end
        
        function ok = setAlphasVec(p, GHC_prioVec)
            GHC_prioMat = p.GHC_Vec2Mat(GHC_prioVec);
            ok = p.setAlphasMat(GHC_prioMat);
        end
        
        function GHC_prioVec = getAlphasVec(p)
            GHC_prioMat = p.getAlphasMat();
            GHC_prioVec = p.GHC_Mat2Vec(GHC_prioMat);
        end

        function ok = setAlphaIJ(p, i,  j, a)
            assert(i >= 0 && i <= p.numTasks);
            assert(j >= 0 && j <= p.numTasks);
            assert(a >= 0 && a <= 1.0);
            if(i==j)
                p.alphas(i,j) = a;
            else
                p.alphas(i,j) = a;
                p.alphas(j,i) = 1.0-a;%TODO makes sense in most scenarios! (but maybe not always?)
            end
            ok = true;
        end

        function a = getAlphaIJ(p, i,  j)
            a = p.alphas(i,j);
        end

        function Ai = getAiVec(p, i)
            Ai = zeros(p.overallTasksize, 1);
            pos = 1;
            for taskNr=1:p.numTasks
                for taskDim=1:p.Tasksize(taskNr)
                    Ai(pos) = p.alphas(i,taskNr);
                    pos = pos + 1;
                end
            end
            assert(pos-1 == p.overallTasksize);
        end

        function Ai = getAiMat(p, i)
            pos = 1;
            for taskNr=1:p.numTasks
                Ai(pos:pos-1+p.Tasksize(taskNr),pos:pos-1+p.Tasksize(taskNr)) = p.alphas(i,taskNr) * eye(p.Tasksize(taskNr));
                pos = pos + p.Tasksize(taskNr);
            end
            assert(pos-1 == p.overallTasksize);
            assert(size(Ai,1) == p.overallTasksize);
            assert(size(Ai,2) == p.overallTasksize);
        end
        
        function GHC_Mat = GHC_Vec2Mat(p,GHC_Vec)
            GHC_Mat = zeros(p.numTasks,p.numTasks);
            row = 1;
            col = 1;
            for i=1:length(GHC_Vec)
%                 assert(1.0 >= GHC_Vec(i));
%                 assert(0.0 <= GHC_Vec(i));
                if row == col
                    GHC_Mat(row,col) = GHC_Vec(i);
                else
                    GHC_Mat(row,col) = GHC_Vec(i);
                    GHC_Mat(col,row) = 1.0-GHC_Vec(i);
                end
                col = col+1;
                if (col>p.numTasks)
                    row=row+1;
                    col=row;
                end
            end
        end

        function GHC_Vec = GHC_Mat2Vec(p,GHC_Mat)
            assert(size(GHC_Mat,1) == size(GHC_Mat,2));
            p.numTasks = size(GHC_Mat,1);
            nbVecElements = p.numTasks + 0.5*(p.numTasks*p.numTasks-p.numTasks);
            GHC_Vec = zeros(nbVecElements,1);
            counter = 1;
            for row=1:p.numTasks
                for col=row:p.numTasks
                    assert(1.0 >= GHC_Mat(row,col));
                    assert(0.0 <= GHC_Mat(row,col));
                    GHC_Vec(counter) = GHC_Mat(row,col);
                    counter = counter +1;
                end
            end
            assert(nbVecElements == length(GHC_Vec));
        end

        function showPrioritization(p)
            % a(i,i) = 0 => task fully activated
            % 0 < a(i,i) < 1 => gradually task (de-)activation
            % a(i,i) = 1 => task fully deactivated (projected in its own nullspace)

            % a(i,j) = 0 => task j has strict lower priority with respect to task i
            % 0 < a(i,j) < 1 => non strict prioritization:
            % a(i,j) = 1 => task j has strict higher priority with respect to task i

            assert(p.numTasks > 0);
            assert(size(p.alphas,1) == p.numTasks);
            assert(size(p.alphas,2) == p.numTasks);
            for i=1:p.numTasks
                for j=1:p.numTasks
                    if(i==j)
                        if(p.alphas(i,j) == 0.0)
                            disp([num2str( i ) " fully activated" ]);
                        elseif(p.alphas(i,j) > 0.0 && p.alphas(i,j) < 1.0)
                            disp([num2str( i ) " gradually de-/activated" ]);
                        elseif(p.alphas(i,j) == 1.0)
                            disp([num2str( i ) " fully deactivated" ]);
                        else
                            disp([num2str( i ) " & " num2str( j ) " is wrong = " num2str( p.alphas(i,j) )]);
                        end
                    else
                        if(p.alphas(i,j) == 0.0)
                            disp([num2str( i ) " >> " num2str( j )]);
                        elseif(p.alphas(i,j) > 0.0 && p.alphas(i,j) < 0.5)
                            disp([num2str( i ) " >  " num2str( j )]);
                        elseif(p.alphas(i,j) == 0.5)
                            disp([num2str( i ) " =  " num2str( j )]);
                        elseif(p.alphas(i,j) > 0.5 && p.alphas(i,j) < 1.0)
                            disp([num2str( i ) " <  " num2str( j )]);
                        elseif(p.alphas(i,j) == 1.0)
                            disp([num2str( i ) " << " num2str( j )]);
                        else
                            disp([num2str( i ) " & " num2str( j ) " is wrong = " num2str( p.alphas(i,j) )]);
                        end
                    end
                end
            end
        end

        function [Q,R,rank,origin]=computeQRD(p, A,varargin)
            if nargin==2
                p.epsilon = realmin('double');
            elseif nargin==3
                p.epsilon = varargin1end;
            end

            p.changeRows = true;
            p.changeRows = false;
            [m,n] = size(A);
            assert(m>0);
            assert(n>0);    
            Q = zeros(m,n);
            R = zeros (n,n);
            origin = zeros(n,1);
            i = 1;
            for k=1:n
                if p.changeRows
                    l=i;
                else
                    l=k;
                end

                if (l >= m+1)
                    break;
                end
                Q(:,l) = A(:,k);
                for j=1:l-1 %TODO -1 correct?
                    R(j,l) = Q(:,j)'*A(:,k);
                    %R(j,l) = Q(:,j)'*Q(:,l);
                    Q(:,l) = Q(:,l) - Q(:,j) * (Q(:,l)' * Q(:,j));
                end
                R(l,l) = norm(Q(:,l));
                if (R(l,l) > p.epsilon)
                    Q(:,l) = Q(:,l) ./ R(l,l);
                    origin(l) = k;
                    i = i+1;
                end
            end
            rank = i-1;
        end

        function ok = checkOrthogonality(p, Mat, varargin)
            if nargin==2
                threshold = realmin('double');
            elseif nargin==3
                threshold = varargin1end;
            end

            %check if each combination of two task directions are orthogonal (Mat(k) * Mat(j)^T == 0)
            assert(threshold>0 && threshold<1);
            ok = true;
            for a=1:size(Mat,2)
                for b=a+1:size(Mat,2)
                    orthogonality = Mat(:,a)' * Mat(:,b);
                    if(isinf(orthogonality))
                        ok = false;
                        disp([ "orthogonality inf-issue for a=" num2str( a ) " and b=" num2str( b ) " => Mat(a)*Mat(b)^T = " num2str( orthogonality )]);
                    end
                    if(orthogonality~=orthogonality)
                        ok = false;
                        disp([ "orthogonality nan-issue for a=" num2str( a ) " and b=" num2str( b ) " => Mat(a)*Mat(b)^T = " num2str( orthogonality )]);
                    end
                    if(orthogonality > threshold || orthogonality < -threshold)
                        ok = false;
                        disp([ "orthogonality issue for a=" num2str( a ) " and b=" num2str( b ) " => Mat(a)*Mat(b)^T = " num2str( orthogonality )]);
                    end
                end
            end
        end
        
        function ok=setPriorities(p, priorities)
            ok=p.setAlphasVec(priorities);
        end
    end
end