/*!
 * @file Physics.h
 *
 * Handles the physics of the Shallow water equations.
 *
 * @authors Gabriel Suau, Remi Pegouret, Lucas Trautmann
 *
 * @version 0.1.0
 *
 * @copyright © 2021 Gabriel Suau
 * @copyright © 2021 Remi Pegouret
 * @copyright © 2021 Lucas Trautmann
 * 
 * @copyright This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * @copyright This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * @copyright You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef PHYSICS_H
#define PHYSICS_H

#include "DataFile.h"
#include "Mesh.h"
#include "termcolor.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

class Physics
{
private:
  // Pointer to useful objects
  DataFile* _DF;
  Mesh* _mesh;

  // Useful variables
  double _g;
  int _nCells;
  Eigen::Matrix<double, Eigen::Dynamic, 2> _cellCenters;

  // Initial condition
  Eigen::Matrix<double, Eigen::Dynamic, 3> _Sol0;
  
  // Topography and source term
  Eigen::VectorXd _topography;
  Eigen::Matrix<double, Eigen::Dynamic, 3> _source;
  
public:
  // Constructeur
  Physics();
  Physics(DataFile* DF, Mesh* mesh);

  // Initialisation
  void Initialize();
  void Initialize(DataFile* DF, Mesh* mesh);

  // Getters
  const Eigen::Matrix<double, Eigen::Dynamic, 3>& getInitialCondition() const {return _Sol0;};
  const Eigen::VectorXd& getTopography() const {return _topography;};
  const Eigen::Matrix<double, Eigen::Dynamic, 3>& getSourceTerm() const {return _source;};
  
  // Construit le terme source
  void buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol);

  // Conditions aux limites
  Eigen::Vector3d dirichletFunction(double x, double y, double t);
  Eigen::Vector3d neumannFunction(double x, double y, double t);

  // Compute the physical flux
  Eigen::Matrix<double, 3, 2> physicalFlux(const Eigen::Vector3d& Sol) const;

  // Compute the eigenvalues of the flux jacobian
  void computeWaveSpeed(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal, double& lambda1, double& lambda2) const;
};

#endif // PHYSICS_H
