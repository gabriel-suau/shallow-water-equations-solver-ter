/*!
 * @file TimeScheme.cpp
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

#include "TimeScheme.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"
#include "FiniteVolume.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


//--------------------------------------------------//
//--------------------Base Class--------------------//
//--------------------------------------------------//
TimeScheme::TimeScheme()
{
}

TimeScheme::TimeScheme(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol):
  _DF(DF), _mesh(mesh), _physics(physics), _finVol(finVol), _Sol(_physics->getInitialCondition()), _timeStep(DF->getTimeStep()), _initialTime(DF->getInitialTime()), _finalTime(DF->getFinalTime()), _currentTime(_initialTime)
{
}

void TimeScheme::Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _finVol = finVol;
  _Sol = _physics->getInitialCondition();
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime;
}

void TimeScheme::saveCurrentSolution(std::string& fileName) const
{
  std::ofstream outputFile(fileName, std::ios::out);
  outputFile.precision(7);
  Eigen::Matrix<double, Eigen::Dynamic, 2> cellCenters (_mesh->getCellsCenter());

  // Vérifications
  if (_Sol.rows() != _mesh->getNumberOfCells())
    {
      std::cout << termcolor::red << "ERROR::TIMESCHEME : The size of the solution is not the same that the number of cells !" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }

  // Informations générales
  outputFile << "# vtk DataFile Version 3.0 " << std::endl;
  outputFile << "2D Unstructured Grid" << std::endl;
  outputFile << "ASCII" << std::endl;
  outputFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // Sauvegarde des sommets
  int nbVertices(_mesh->getNumberOfVertices());
  outputFile << "POINTS " << nbVertices << " float " << std::endl;
  for (int i(0) ; i < nbVertices ; ++i)
    {
      outputFile << _mesh->getVertices()[i].getCoordinates()[0] << " " << _mesh->getVertices()[i].getCoordinates()[1] << " 0." << std::endl;
    }
  outputFile << std::endl;

  // Sauvegarde des cellules
  int nbCells(_mesh->getNumberOfCells());
  int nbVerticesPCell(_mesh->getNumberOfVerticesPerCell());
  outputFile << "CELLS " << nbCells << " " << nbCells * (nbVerticesPCell + 1) << std::endl;
  for (int i(0) ; i < nbCells ; ++i)
    {
      outputFile << nbVerticesPCell;
      for (int j(0) ; j < nbVerticesPCell ; ++j)
        {
          outputFile << " " << _mesh->getCells()[i].getVerticesIndex()[j];
        }
      outputFile << std::endl;
    }
  outputFile << std::endl;

  // Sauvegarde du type de cellules
  outputFile << "CELL_TYPES " << nbCells << std::endl;
  for (int i(0) ; i < nbCells ; ++i)
    {
      outputFile << 5 << std::endl;
    }
  outputFile << std::endl;

  outputFile << "CELL_DATA " << nbCells << std::endl;

  // Sauvegarde de la hauteur
  outputFile << "SCALARS h float 1" << std::endl;
  outputFile << "LOOKUP_TABLE default" << std::endl;
  for (int i(0) ; i < nbCells ; ++i)
    {
      outputFile << _Sol(i,0) << std::endl;
    }
  outputFile << std::endl;

  // Sauvegarde de la vitesse
  outputFile << "VECTORS vel float" << std::endl;
  for (int i(0) ; i < _mesh->getNumberOfCells() ; ++i)
    {
      outputFile << _Sol(i,1)/_Sol(i,0) << " " << _Sol(i,2)/_Sol(i,0) << " 0" << std::endl;
    }
  outputFile << std::endl;
}

void TimeScheme::solve()
{
  // Logs de début
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Time loop..." << std::endl;

  // Variables pratiques
  int n(0);
  std::string resultsDir(_DF->getResultsDirectory());
  std::string fluxName(_finVol->getFluxName());

  // Sauvegarde la condition initiale
  std::string fileName(resultsDir + "/solution_" + fluxName + "_" + std::to_string(n) + ".vtk");
  saveCurrentSolution(fileName);

  // Boucle en temps
  while (_currentTime < _finalTime)
    {
      _physics->buildSourceTerm(_Sol);
      _finVol->buildFluxVector(_Sol);
      oneStep();
      ++n;
      _currentTime += _timeStep;
      if (n % _DF->getSaveFrequency() == 0)
        {
          std::cout << "Saving solution at t = " << _currentTime << std::endl;
          std::string fileName(resultsDir + "/solution_" + fluxName + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".vtk");
          saveCurrentSolution(fileName);
        }
    }

  // Logs de fin
  std::cout << termcolor::green << "SUCCESS::TIMESCHEME : Solved 2D St-Venant equations successfully !" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}


//-------------------------------------------------------------//
//--------------------Explicit Euler scheme--------------------//
//-------------------------------------------------------------//
ExplicitEuler::ExplicitEuler():
  TimeScheme()
{
}

ExplicitEuler::ExplicitEuler(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol):
  TimeScheme(DF, mesh, physics, finVol)
{
}

void ExplicitEuler::Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _finVol = finVol;
  _Sol = _physics->getInitialCondition();
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime;
}

void ExplicitEuler::oneStep()
{
  // Récupération des trucs importants
  double dt(_timeStep);
  const Eigen::VectorXd& cellsArea(_mesh->getCellsArea());
  const Eigen::Matrix<double, Eigen::Dynamic, 3>& fluxVector(_finVol->getFluxVector());
  // const Eigen::Matrix<double, Eigen::Dynamic, 3>& sourceTerm(_physics->getSourceTerm());
  
  // Mise à jour de la solution
  for (int i(0) ; i < _Sol.rows() ; ++i)
    {
      double cellArea(cellsArea(i));
      _Sol.row(i) += - dt / cellArea * fluxVector.row(i);
    }
}
