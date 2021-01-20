#include "TimeScheme.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Function.h"
#include "FiniteVolume.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

TimeScheme::TimeScheme()
{
}

TimeScheme::TimeScheme(DataFile* DF, Mesh* mesh, Function* function, FiniteVolume* finVol):
  _DF(DF), _mesh(mesh), _function(function), _finVol(finVol), _Sol(_function->getInitialCondition()), _timeStep(DF->getTimeStep()), _initialTime(DF->getInitialTime()), _finalTime(DF->getFinalTime()), _currentTime(_initialTime)
{
}

void TimeScheme::Initialize(DataFile* DF, Mesh* mesh, Function* function, FiniteVolume* finVol)
{
  _DF = DF;
  _mesh = mesh;
  _function = function;
  _finVol = finVol;
  _Sol = _function->getInitialCondition();
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime;
}

void TimeScheme::saveCurrentSolution(std::string& fileName) const
{
  std::ofstream outputFile(fileName, std::ios::out);
  outputFile.precision(7);
  Eigen::Matrix<double, Eigen::Dynamic, 2> cellCenters (_mesh->getTrianglesCenter());

  // Vérifications
  if (_Sol.size() != _mesh->getNumberOfTriangles())
    {
      std::cout << termcolor::red << "ERROR::TIMESCHEME : The size of the solution is not the same that the number of triangles !" << std::endl;
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
  int nbTriangles(_mesh->getNumberOfTriangles());
  outputFile << "CELLS " << nbTriangles << " " << nbTriangles * 4 << std::endl;
  for (int i(0) ; i < nbTriangles ; ++i)
    {
      outputFile << 3 << _mesh->getTriangles()[i].getVerticesReference()[0]
                 << " " << _mesh->getTriangles()[i].getVerticesReference()[1]
                 << " " << _mesh->getTriangles()[i].getVerticesReference()[2] << std::endl;
    }
  outputFile << std::endl;

  // Sauvegarde du type de cellules
  outputFile << "CELL_TYPES " << nbTriangles << std::endl;
  for (int i(0) ; i < nbTriangles ; ++i)
    {
      outputFile << 5 << std::endl;
    }
  outputFile << std::endl;

  outputFile << "CELL_DATA " << nbTriangles << std::endl;

  // Sauvegarde de la hauteur
  outputFile << "SCALARS h float 1" << std::endl;
  outputFile << "LOOKUP_TABLE default" << std::endl;
  for (int i(0) ; i < nbTriangles ; ++i)
    {
      outputFile << _Sol(i,0) << std::endl;
    }
  outputFile << std::endl;

  // Sauvegarde de u
  outputFile << "SCALARS u float 1" << std::endl;
  outputFile << "LOOKUP_TABLE default" << std::endl;
  for (int i(0) ; i < nbTriangles ; ++i)
    {
      if (_Sol(i,0) > 1e-10)
        {
          outputFile << _Sol(i,1)/_Sol(i,0) << std::endl; 
        }
      else
        {
          outputFile << 0. << std::endl;
        }
    }
  outputFile << std::endl;

  // Sauvegarde de v
  outputFile << "SCALARS v float 1" << std::endl;
  outputFile << "LOOKUP_TABLE default" << std::endl;
  for (int i(0) ; i < nbTriangles ; ++i)
    {
      if (_Sol(i,0) > 1e-10)
        {
          outputFile << _Sol(i,2)/_Sol(i,0) << std::endl; 
        }
      else
        {
          outputFile << 0. << std::endl;
        }
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
      _function->buildSourceTerm(_Sol);
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

ExplicitEuler::ExplicitEuler():
  TimeScheme()
{
}

ExplicitEuler::ExplicitEuler(DataFile* DF, Mesh* mesh, Function* function, FiniteVolume* finVol):
  TimeScheme(DF, mesh, function, finVol)
{
}

void ExplicitEuler::Initialize(DataFile* DF, Mesh* mesh, Function* function, FiniteVolume* finVol)
{
  _DF = DF;
  _mesh = mesh;
  _function = function;
  _finVol = finVol;
  _Sol.resize(mesh->getNumberOfTriangles(), 2);
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime; 
}

void ExplicitEuler::oneStep()
{
  // Récupération des trucs importants
  double dt(_timeStep);

  // TODO
}
