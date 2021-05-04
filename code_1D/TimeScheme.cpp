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


//----------------------------------------------------------//
//------------------Time Scheme base class------------------//
//----------------------------------------------------------//
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
#if VERBOSITY>0
  std::cout << "Saving solution at t = " << _currentTime << std::endl;
#endif
  std::ofstream outputFile(fileName, std::ios::out);
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
  double g(_DF->getGravityAcceleration());
  // Gnuplot comments for the user
  outputFile << "# x  H=h+z   h       u       q       Fr=|u|/sqrt(gh)" << std::endl;
  for (int i(0) ; i < _Sol.rows() ; ++i)
    {
      outputFile << cellCenters(i) << " " <<
        _Sol(i,0) + _physics->getTopography()(i) << " " <<
        _Sol(i,0) << " " <<
        _Sol(i,1)/_Sol(i,0) << " " <<
        _Sol(i,1) << " " <<
        abs(_Sol(i,1)/_Sol(i,0))/sqrt(g * _Sol(i,0)) << std::endl;
    }
}



void TimeScheme::solve()
{
  // Logs de début
#if VERBOSITY>0
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Time loop..." << std::endl;
#endif
  
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
      topoFile << _mesh->getCellCenters()(i) << " " << _physics->getTopography()(i) << std::endl;
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
  if (_DF->isTestCase())
    {
      _physics->buildExactSolution(_currentTime);
      std::string fileName(resultsDir + "/solution_exacte.txt");
      _physics->saveExactSolution(fileName);
      double L2errorHeight(computeL2Error()(0));
      double L2errorDischarge(computeL2Error()(1));
      std::cout << "Error h  L2 = " << L2errorHeight << " and error q L2 = " << L2errorDischarge << " for dx = " << _DF->getDx() << std::endl;
      double L1errorHeight(computeL1Error()(0));
      double L1errorDischarge(computeL1Error()(1));
      std::cout << "Error h  L1 = " << L1errorHeight << " and error q L1 = " << L1errorDischarge << " for dx = " << _DF->getDx() << std::endl;
    }
  // Logs de fin
#if VERBOSITY>0
  std::cout << termcolor::green << "TIMESCHEME::SUCCESS : Solved 1D St-Venant equations successfully !" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
#endif
}


Eigen::Vector2d TimeScheme::computeL2Error() const
{
  Eigen::Vector2d error(0., 0.);
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& exactSol(_physics->getExactSolution());
  error(0) = (_Sol.col(0) - exactSol.col(0)).norm();
  error(1) = (_Sol.col(1) - exactSol.col(1)).norm();
  error *= sqrt(_DF->getDx());
  return error;
}


Eigen::Vector2d TimeScheme::computeL1Error() const
{
  Eigen::Vector2d error(0., 0.);
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& exactSol(_physics->getExactSolution());
  for (int i(0) ; i < _Sol.rows() ; ++i)
    {
      error(0) += abs(_Sol(i,0) - exactSol(i,0));
      error(1) += abs(_Sol(i,1) - exactSol(i,1));
    }
  error *= _DF->getDx();
  return error;
}


//--------------------------------------------------//
//------------------Explicit Euler------------------//
//--------------------------------------------------//
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

  // Construction du terme source et du flux numérique
  _finVol->buildFluxVector(_currentTime, _Sol);
  _physics->buildSourceTerm(_Sol);
  // Recuperation du terme source et du flux numerique
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& source(_physics->getSourceTerm());
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



RK2::RK2(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol):
  TimeScheme(DF, mesh, physics, finVol)
{
}



void RK2::Initialize(DataFile* DF, Mesh* mesh, Physics* physics, FiniteVolume* finVol)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
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
  _physics->buildSourceTerm(_Sol);
  _finVol->buildFluxVector(_currentTime, _Sol);
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& source1(_physics->getSourceTerm());
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& fluxVector1(_finVol->getFluxVector());
  k1 = fluxVector1 / dx + source1;
  // Calcul de k2
  _physics->buildSourceTerm(_Sol + dt * k1);
  _finVol->buildFluxVector(_currentTime + dt, _Sol + dt * k1);
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& source2(_physics->getSourceTerm());
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& fluxVector2(_finVol->getFluxVector());
  k2 = fluxVector2 / dx + source2;
  // Mise a jour de la solution
  _Sol += 0.5 * dt * (k1 + k2);
}
