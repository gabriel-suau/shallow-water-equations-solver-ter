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


//----------------------------------------------------------//
//------------------Time Scheme base class------------------//
//----------------------------------------------------------//
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
  std::cout << "Saving solution at t = " << _currentTime << std::endl;
  std::ofstream outputFile(fileName, std::ios::out);
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
  double g(_DF->getGravityAcceleration());
  // Gnuplot comments for the user
  outputFile << "# x  H=h+z   h       u       q       Fr=|u|/sqrt(gh)" << std::endl;
  for (int i(0) ; i < _Sol.rows() ; ++i)
    {
      outputFile << cellCenters(i) << " " <<
        _Sol(i,0) + _function->getTopography()(i,1) << " " <<
        _Sol(i,0) << " " <<
        _Sol(i,1)/_Sol(i,0) << " " <<
        _Sol(i,1) << " " <<
        abs(_Sol(i,1)/_Sol(i,0))/sqrt(g * _Sol(i,0)) << std::endl;
    }
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
  std::string fileName(resultsDir + "/solution_" + fluxName + "_" + std::to_string(n) + ".txt");
  saveCurrentSolution(fileName);

  // Sauvegarde la topographie
  std::string topoFileName(resultsDir + "/topography.txt");
  std::ofstream topoFile(topoFileName, std::ios::out);
  for (int i(0) ; i < _Sol.rows() ; ++i)
    {
      topoFile << _function->getTopography()(i,0) << " " << _function->getTopography()(i,1) << std::endl;
    }
  
  // Boucle en temps
  while (_currentTime < _finalTime)
    {
      oneStep();
      ++n;
      _currentTime += _timeStep;
      if (!_DF->isSaveFinalTimeOnly() &&  n % _DF->getSaveFrequency() == 0)
        {
          std::string fileName(resultsDir + "/solution_" + fluxName + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".txt");
          saveCurrentSolution(fileName);
        }
    }
  if (_DF->isSaveFinalTimeOnly())
    {
      std::string fileName(resultsDir + "/solution_" + fluxName + "_" + std::to_string(n/_DF->getSaveFrequency()) + ".txt");
      saveCurrentSolution(fileName);
    }
  
  // Logs de fin
  std::cout << termcolor::green << "TIMESCHEME::SUCCESS : Solved 1D St-Venant equations successfully !" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}


//--------------------------------------------------//
//------------------Explicit Euler------------------//
//--------------------------------------------------//
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
  _Sol.resize(mesh->getNumberOfCells(), 2);
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime; 
}



void ExplicitEuler::oneStep()
{
  // Récupération des trucs importants
  double dt(_timeStep);
  double dx(_mesh->getSpaceStep());
  int nCells(_mesh->getNumberOfCells());

  // Construction du terme source et du flux numérique
  _function->buildSourceTerm(_Sol);
  _finVol->buildFluxVector(_currentTime, _Sol);
  // Recuperation du terme source et du flux numerique
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& source(_function->getSourceTerm());
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& fluxVector(_finVol->getFluxVector());

  // Mise à jour de la solution sur chaque cellules
  _Sol += dt * (fluxVector / dx + source);
}


//-------------------------------------------------//
//------------------Runge Kutta 2------------------//
//-------------------------------------------------//
RK2::RK2():
  TimeScheme()
{
}



RK2::RK2(DataFile* DF, Mesh* mesh, Function* function, FiniteVolume* finVol):
  TimeScheme(DF, mesh, function, finVol)
{
}



void RK2::Initialize(DataFile* DF, Mesh* mesh, Function* function, FiniteVolume* finVol)
{
  _DF = DF;
  _mesh = mesh;
  _function = function;
  _finVol = finVol;
  _Sol.resize(mesh->getNumberOfCells(), 2);
  _timeStep = DF->getTimeStep();
  _initialTime = DF->getInitialTime();
  _finalTime = DF->getFinalTime();
  _currentTime = _initialTime; 
}



void RK2::oneStep()
{
  // Récupération des trucs importants
  double dt(_timeStep);
  double dx(_mesh->getSpaceStep());

  Eigen::Matrix<double, Eigen::Dynamic, 2> k1, k2;

  // Calcul de k1
  _function->buildSourceTerm(_Sol);
  _finVol->buildFluxVector(_currentTime, _Sol);
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& source1(_function->getSourceTerm());
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& fluxVector1(_finVol->getFluxVector());
  k1 = fluxVector1 / dx + source1;
  // Calcul de k2
  _function->buildSourceTerm(_Sol + dt * k1);
  _finVol->buildFluxVector(_currentTime + dt, _Sol + dt * k1);
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& source2(_function->getSourceTerm());
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& fluxVector2(_finVol->getFluxVector());
  k2 = fluxVector2 / dx + source2;
  // Mise a jour de la solution
  _Sol += 0.5 * dt * (k1 + k2);
}
