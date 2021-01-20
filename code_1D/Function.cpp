#include "Function.h"
#include "DataFile.h"
#include "termcolor.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>

Function::Function()
{
}

Function::Function(DataFile* DF, Mesh* mesh):
  _DF(DF), _mesh(mesh), _xmin(mesh->getxMin()), _xmax(mesh->getxMax()), _g(_DF->getGravityAcceleration()), _nCells(mesh->getNumberOfCells()), _cellCenters(mesh->getCellCenters())
{
}

void Function::Initialize(DataFile* DF, Mesh* mesh)
{
  // Initialisation
  _DF = DF;
  _mesh = mesh;
  _xmin = mesh->getxMin();
  _xmax = mesh->getxMax();
  _g = DF->getGravityAcceleration();
  _nCells = mesh->getNumberOfCells();
  _cellCenters = mesh->getCellCenters();  
  this->Initialize();
}

void Function::Initialize()
{
  // Logs de début
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Building topography and initial condition..." << std::endl;
  
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
  else if (_DF->getTopographyType() == "LinearUp")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i,1) = 0.05 * (_cellCenters(i) - _xmin);
        }
    }
  else if (_DF->getTopographyType() == "LinearDown")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i,1) = 0.05 * (_xmax - _cellCenters(i));
        }
    }
  else if (_DF->getTopographyType() == "SineLinearUp")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i,1) = 0.05 * (_cellCenters(i) - _xmin) + 0.05*sin(20*M_PI*_cellCenters(i)/(_xmax - _xmin));
        }
    }
  else if (_DF->getTopographyType() == "SineLinearDown")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i,1) = 0.05 * (_xmax - _cellCenters(i)) + 0.05*sin(20*M_PI*_cellCenters(i)/(_xmax - _xmin));
        }
    }
  else if (_DF->getTopographyType() == "EllipticBump")
    {
      _topography.col(1).setZero();
      double bumpCenter(_xmin + 0.75 * (_xmax - _xmin));
      double bumpHeight(0.8);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double bumpFunction(bumpHeight - pow(_cellCenters(i) - bumpCenter,2));
          if (bumpFunction > 0)
            {
              _topography(i,1) = bumpFunction;
            }
        }
    }
  else if (_DF->getTopographyType() == "File")
    {
      const std::string topoFile(_DF->getTopographyFile());
      std::ifstream topoStream(topoFile);
      std::string line, properLine;
      double dummy;
      int i(0);
      if (!topoStream.is_open())
        {
          std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Unable to open the topography file : " << topoFile << std::endl;
          std::cout << termcolor::reset << "====================================================================================================" << std::endl;
          exit(-1);
        }
      else
        {
          std::cout << "Building the topography from file : " << topoFile << std::endl;
        }
      while(getline(topoStream, line))
        {
          properLine = regex_replace(line, std::regex(",") , std::string(" "));
          std::stringstream ss(properLine);
          ss >> dummy >> _topography(i,1);
          ++i;
        }
    }
  else
    {
      std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Case not implemented" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
  std::cout << termcolor::green << "SUCCESS::TOPOGRAPHY : Topography was successfully built." << std::endl;
  std::cout << termcolor::reset;

  // Initialise la condition initiale
  if (_DF->getScenario() == "ConstantWaterHeight")
    {
      _Sol0.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = 3.;
        }
    }
  else if (_DF->getScenario() == "RestingLake")
    {
      _Sol0.col(1).setZero();
      double H(3.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = std::max(H - _topography(i,1), 0.);
        }
    }
  else if (_DF->getScenario() == "DamBreak")
    {
      _Sol0.col(1).setZero();
      double Hg(2.), Hd(1.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (_cellCenters(i) < 0.5*(_xmax + _xmin))
            {
              _Sol0(i,0) = std::max(Hg - _topography(i,1), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(Hd - _topography(i,1), 0.);
            }
        }
    }
  else if (_DF->getScenario() == "SineWave")
    {
      _Sol0.col(1).setZero();
      double H(2.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = std::max(H + 0.2 * cos(M_PI * _cellCenters(i)) - _topography(i,1), 0.);
        }
    }
  else if (_DF->getScenario() == "SinePerturbation")
    {
      _Sol0.col(1).setZero();
      double H(2.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (-1 < _cellCenters(i) && _cellCenters(i) < 1)
            {
              _Sol0(i,0) = std::max(H + 0.2 * cos(M_PI * _cellCenters(i)) - _topography(i,1), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(1.8 - _topography(i,1),0.);
            }
        }
    }
  else
    {
      std::cout << termcolor::red << "ERROR::SCENARIO : Case not implemented" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
  std::cout << termcolor::green << "SUCCESS::SCENARIO : Initial Conditions was successfully built." << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}

void Function::buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  // Construit le terme source en fonction de la topographie.
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _source.setZero();
    }
  else if (_DF->getTopographyType() == "LinearUp")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * 0.05;
        }
    }
  else if (_DF->getTopographyType() == "LinearDown")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _source(i,1) = _g * Sol(i,0) * 0.05;
        }
    }
  else if (_DF->getTopographyType() == "SineLinearUp")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * (0.05 + M_PI/(_xmax - _xmin)*cos(20*M_PI*_cellCenters(i)/(_xmax - _xmin)));
        }
    }
  else if (_DF->getTopographyType() == "SineLinearDown")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * (-0.05 + M_PI/(_xmax - _xmin)*cos(20*M_PI*_cellCenters(i)/(_xmax - _xmin)));
        }
    }
  else if (_DF->getTopographyType() == "EllipticBump")
    {
      _source.setZero();
      double bumpCenter(_xmin + 0.75 * (_xmax - _xmin));
      double bumpHeight(0.8);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double bumpFunction(bumpHeight - pow(_cellCenters(i) - bumpCenter,2));
          if (bumpFunction > 0)
            {
              _source(i,1) = _g * Sol(i,0) * 2. * (_cellCenters(i) - bumpCenter);
            }
        }
    }

  // Pour un fichier de topographie, la dérivée est approchée par une formule de
  // différence finie centrée d'ordre deux à l'intérieur, et par une formule
  // décentrée d'ordre deux sur les bords
  else if (_DF->getTopographyType() == "File")
    {
      double dx(_mesh->getSpaceStep());
      _source.col(0).setZero();
      _source(0,1) = - _g * Sol(0,0) * (-_topography(2) + 4.*_topography(1) - 3.*_topography(0))/(2.*dx);
      for (int i(1) ; i < _nCells - 1 ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * (_topography(i+1) - _topography(i-1))/(2. * dx);
        }
      _source(_nCells - 1, 1) = - _g * Sol(_nCells - 1,0) * (3.*_topography(_nCells - 1) - 4.*_topography(_nCells - 2) + _topography(_nCells - 3))/(2.*dx);
    }
  else
    {
      std::cout << termcolor::red << "ERROR::SOURCETERM : Case not implemented." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
}

Eigen::Vector2d Function::dirichletFunction(double x, double t)
{
  Eigen::Vector2d g(0.,0.);
  // Condition d'entrée
  if (x == _mesh->getxMin())
    {
      g << 2., 10.;
    }
  // Condition de sortie
  if (x == _mesh->getxMax())
    {
      
    }
  return g;
}

Eigen::Vector2d Function::neumannFunction(double x, double t)
{
  Eigen::Vector2d h(0.,0.);
  // Condition d'entrée
  
  // Condition de sortie

  return h;
}
