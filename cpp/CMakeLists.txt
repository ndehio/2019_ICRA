# Author: Niels Dehio
#
#
# Writing code takes time. Polishing it and making it available to others takes longer! 
# If some parts of the code were useful for your research or for a better understanding 
# of the algorithms, please reward the authors by citing the related publications, 
# and consider making your own research available in this way.
#
# @inproceedings{Dehio2018b,
# title={{Continuously Shaping Projections and Operational Space Tasks}},
# author={Dehio, Niels and Kubus, Daniel and Steil, Jochen J.},
# booktitle={IEEE/RSJ Int. Conf. on Intelligent Robots and Systems},
# pages={5995--6002}, 
# year={2018}
# }
# 
# @article{Dehio2018c,
# title={{Prioritized Multi-Objective Robot Control}},
# author={Dehio, Niels},
# journal={PhD dissertation, Technical University of Braunschweig, Germany},
# year={2018},
# url={https://publikationsserver.tu-braunschweig.de/receive/dbbs_mods_00066108}
# }
# 
# @inproceedings{Dehio2019,
# title={{Dynamically-consistent Generalized Hierarchical Control}},
# author={Dehio, Niels and Steil, Jochen J.},
# booktitle={IEEE/RSJ Int. Conf. on Robotics and Automation},
# year={2019}
# }
# 
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as
# published by the Free Software Foundation.
# 
# This code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#
# #########################################################################

project(DynGHC)
cmake_minimum_required (VERSION 2.8)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

add_executable(Demo src/demo.cpp src/GHCProjections.hpp src/GHCProjections.cpp)
