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
  Eigen::Matrix<double, Eigen::Dynamic, 2> cellCenters (_mesh->getTrianglesCenter());

  // TODO
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

  // Sauvegarde la topographie
  std::string topoFileName(resultsDir + "/topography.txt");
  std::ofstream topoFile(topoFileName, std::ios::out);
  for (int i(0) ; i < _Sol.rows() ; ++i)
    {
      topoFile << _function->getTopography()(i,0) << " " << _function->getTopography()(i,1) << _function->getTopography()(i,2) << std::endl;
    }
  
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
  std::cout << termcolor::green << "TIMESCHEME::SUCCESS : Solved 1D St-Venant equations successfully !" << std::endl;
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
