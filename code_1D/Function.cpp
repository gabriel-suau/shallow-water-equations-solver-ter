#include "Function.h"
#include "DataFile.h"
#include "termcolor.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <cmath>

Function::Function()
{
}

Function::Function(DataFile* DF, Mesh* mesh):
  _DF(DF), _mesh(mesh), _xmin(mesh->getxMin()), _xmax(mesh->getxMax()), _nCells(mesh->getNumberOfCells()), _cellCenters(mesh->getCellCenters())
{
}

void Function::Initialize(DataFile* DF, Mesh* mesh)
{
  // Logs de début
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Building topography and initial condition..." << std::endl;

  // Initialisation
  _DF = DF;
  _mesh = mesh;
  _xmin = mesh->getxMin();
  _xmax = mesh->getxMax();
  _nCells = mesh->getNumberOfCells();
  _cellCenters = mesh->getCellCenters();
  
  // Resize la condition initiale, la topographie et le terme source
  _Sol0.resize(_nCells, 2);
  _topography.resize(_nCells, 2);
  _source.resize(_nCells, 2);

  // Initialise la topographie
  _topography.col(0) = _cellCenters;
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _topography.col(1).setZero();
    }
  else
    {
      std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Case not implemented" << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }

  // Initialise la condition initiale
  if (_DF->getScenario() == "ConstantWaterHeight")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0.row(i)(0) = 3.;
        }
    }
  else if (_DF->getScenario() == "RestingLake")
    {
      _Sol0.col(1).setZero();
      double H(3.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0.row(i)(0) = std::max(H - _topography.row(i)(1), 0.);
        }
    }
  else if (_DF->getScenario() == "DamBreak")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (_cellCenters(i) < 0.5*(_xmax + _xmin))
            {
              _Sol0.row(i)(0) = 1.0;
            }
          else
            {
              _Sol0.row(i)(0) = 0.;
            }
        }
    }
  else if (_DF->getScenario() == "SineWave")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0.row(i)(0) = 2. + 0.2 * cos(M_PI * _cellCenters(i));
        }
    }
  else if (_DF->getScenario() == "SinePerturbation")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (-1 < _cellCenters(i) && _cellCenters(i) < 1)
            {
              _Sol0.row(i)(0) = 2. + 0.2 * cos(M_PI * _cellCenters(i)); 
            }
          else
            {
              _Sol0.row(i)(0) = 1.8;
            }
        }
    }
  else
    {
      std::cout << termcolor::red << "ERROR::SCENARIO : Case not implemented" << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }

  // Logs de fin
  std::cout << termcolor::green << "SUCCESS : Topography and Initial Conditions were successfully built." << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}

void Function::Initialize()
{
  // Logs de début
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Building topography and initial condition..." << std::endl;

  // Initialisation
  _xmin = _mesh->getxMin();
  _xmax = _mesh->getxMax();
  _nCells = _mesh->getNumberOfCells();
  _cellCenters = _mesh->getCellCenters();
  
  // Resize la condition initiale, la topographie et le terme source
  _Sol0.resize(_nCells, 2);
  _topography.resize(_nCells, 2);
  _source.resize(_nCells, 2);

  // Initialise la topographie
  _topography.col(0) = _cellCenters;
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _topography.col(1).setZero();
    }
  else
    {
      std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Case not implemented" << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }

  // Initialise la condition initiale
  if (_DF->getScenario() == "ConstantWaterHeight")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0.row(i)(0) = 3.;
        }
    }
  else if (_DF->getScenario() == "RestingLake")
    {
      _Sol0.col(1).setZero();
      double H(3.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0.row(i)(0) = std::max(H - _topography.row(i)(1), 0.);
        }
    }
  else if (_DF->getScenario() == "DamBreak")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (_cellCenters(i) < 0.5*(_xmax + _xmin))
            {
              _Sol0.row(i)(0) = 1.0;
            }
          else
            {
              _Sol0.row(i)(0) = 0.;
            }
        }
    }
  else if (_DF->getScenario() == "SineWave")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0.row(i)(0) = 2. + 0.2 * cos(M_PI * _cellCenters(i));
        }
    }
  else if (_DF->getScenario() == "SinePerturbation")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (-1 < _cellCenters(i) && _cellCenters(i) < 1)
            {
              _Sol0.row(i)(0) = 2. + 0.2 * cos(M_PI * _cellCenters(i)); 
            }
          else
            {
              _Sol0.row(i)(0) = 1.8;
            }
        }
    }
  else
    {
      std::cout << termcolor::red << "ERROR::SCENARIO : Case not implemented" << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }

  // Logs de fin
  std::cout << termcolor::green << "SUCCESS : Topography and Initial Conditions were successfully built." << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}

void Function::buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  // Récupère la valeur de la gravité
  double g(_DF->getGravityAcceleration());

  // Construit le terme source en fonction de la topographie.
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _source.setZero();
    }
  else
    {
      std::cout << termcolor::red << "ERROR::SOURCETERM : Case not implemented." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
}
