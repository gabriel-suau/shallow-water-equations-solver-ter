/*!
 * @file TimeScheme.h
 *
 * Defines classes for time integration.
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

#ifndef TIME_SCHEME_H
#define TIME_SCHEME_H

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"
#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"
#include "FiniteVolume.h"

class TimeScheme
{
protected:
  // Pointeur vers les trucs importants
  DataFile* _DF;
  Mesh* _mesh;
  Physics* _physics;
  FiniteVolume* _finVol;

  // Solution
  Eigen::Matrix<double, Eigen::Dynamic, 3> _Sol;
  
  // Paramètres de temps
  double _timeStep;
  double _initialTime;
  double _finalTime;
  double _currentTime;
  
public:
  // Constructeurs
  TimeScheme();
  TimeScheme(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);

  // Initialiseur
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);
  // Destructeur
  virtual ~TimeScheme() = default;

  // Getters
  const Eigen::Matrix<double, Eigen::Dynamic, 3>& getSolution() const {return _Sol;};
  double getTimeStep() const {return _timeStep;};
  double getInitialTime() const {return _initialTime;};
  double getFinalTime() const {return _finalTime;};
  double getCurrentTime() const {return _currentTime;};
  
  // Solve and save solution
  virtual void oneStep() = 0;
  void saveCurrentSolution(std::string& fileName) const;
  void solve();
};

class ExplicitEuler: public TimeScheme
{
public:
  // Constructeurs
  ExplicitEuler();
  ExplicitEuler(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);

  // Initialiseur
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol);

  // One time step
  void oneStep();
};

#endif // TIME_SCHEME_H
