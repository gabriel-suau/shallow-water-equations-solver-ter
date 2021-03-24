/*!
 * @file FiniteVolume.h
 *
 * Defines classes for the discrete fluxes.
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

#ifndef FINITE_VOLUME_H
#define FINITE_VOLUME_H

#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

//--------------------------------------------------//
//--------------------Base Class--------------------//
//--------------------------------------------------//
class FiniteVolume
{
protected:
  // Pointeurs vers les trucs importants
  DataFile* _DF;
  Mesh* _mesh;
  Physics* _physics;

  // Nom du flux numérique
  std::string _fluxName;

  // Vecteur des flux
  Eigen::Matrix<double, Eigen::Dynamic, 3> _fluxVector;
  
public:
  // Constructeurs
  FiniteVolume();
  FiniteVolume(DataFile* DF, Mesh* mesh, Physics* physics);

  // Destructeur
  virtual ~FiniteVolume() = default;
  
  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Getters
  const std::string& getFluxName() const {return _fluxName;};
  const Eigen::Matrix<double, Eigen::Dynamic, 3>& getFluxVector() const {return _fluxVector;};
  
  // Fluxes
  virtual Eigen::Vector3d numFlux1D(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal) const = 0;
  virtual void buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol) = 0;
};


//--------------------------------------------------//
//-------------------Rusanov flux-------------------//
//--------------------------------------------------//
class Rusanov: public FiniteVolume
{
public:
  // Constructeur
  Rusanov();
  Rusanov(DataFile* DF, Mesh* mesh, Physics* physics);

  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Build flux vector
  void buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol);
  Eigen::Vector3d numFlux1D(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal) const;
};


//--------------------------------------------------//
//---------------------HLL flux---------------------//
//--------------------------------------------------//
class HLL: public FiniteVolume
{
public:
  // Constructeur
  HLL();
  HLL(DataFile* DF, Mesh* mesh, Physics* physics);

  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Build flux vector
  void buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol);
  Eigen::Vector3d numFlux1D(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal) const;
};

#endif //FINITE_VOLUME_H
