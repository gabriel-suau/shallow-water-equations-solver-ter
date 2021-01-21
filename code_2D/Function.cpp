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
  _DF(DF), _mesh(mesh), _g(_DF->getGravityAcceleration()), _nCells(_mesh->getNumberOfTriangles()), _cellCenters(_mesh->getTrianglesCenter())
{
}

void Function::Initialize(DataFile* DF, Mesh* mesh)
{
  // Initialisation
  _DF = DF;
  _mesh = mesh;
  _g = DF->getGravityAcceleration();
  _nCells = _mesh->getNumberOfTriangles();
  _cellCenters = _mesh->getTrianglesCenter();
  this->Initialize();
}

void Function::Initialize()
{
  // Logs de début
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Building topography and initial condition..." << std::endl;
  
  // Resize la condition initiale, la topographie et le terme source
  _Sol0.resize(_nCells, 3);
  _topography.resize(_nCells);
  _source.resize(_nCells, 3);

  // Initialise la topographie
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _topography.setZero();
    }
  else if (_DF->getTopographyType() == "LinearUp")
    {
      // TODO
    }
  else if (_DF->getTopographyType() == "LinearDown")
    {
      // TODO
    }
  else if (_DF->getTopographyType() == "SineLinearDown")
    {
      // TODO
    }
  else if (_DF->getTopographyType() == "EllipticBump")
    {
      // TODO
    }
  else if (_DF->getTopographyType() == "File")
    {
      // TODO
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
      _Sol0.leftCols(2).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = 3.;
        }
    }
  else if (_DF->getScenario() == "RestingLake")
    {
      _Sol0.leftCols(2).setZero();
      double H(3.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = std::max(H - _topography(i,1), 0.);
        }
    }
  else if (_DF->getScenario() == "DamBreak")
    {
      for (int i(0); i < _nCells; i++)
        {
          if (_mesh->getTrianglesCenter()(i,1) < 0.5)
            {
              _Sol0(i,0) = 2;
            }
          else
            {
              _Sol0(i,0) = 1;
            }
        }
    }
  else if (_DF->getScenario() == "SineWave")
    {
      // TODO
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

void Function::buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol)
{
  // Construit le terme source en fonction de la topographie.
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _source.setZero();
    }
  else if (_DF->getTopographyType() == "LinearUp")
    {
      // TODO
    }
  else if (_DF->getTopographyType() == "LinearDown")
    {
      // TODO
    }
  else if (_DF->getTopographyType() == "SineLinearUp")
    {
      // TODO
    }
  else if (_DF->getTopographyType() == "SineLinearDown")
    {
      // TODO
    }
  else if (_DF->getTopographyType() == "EllipticBump")
    {
      // TODO
    }

  // Pour un fichier de topographie, la dérivée est approchée par une formule de
  // différence finie ?
  else if (_DF->getTopographyType() == "File")
    {
      // TODO
    }
  else
    {
      std::cout << termcolor::red << "ERROR::SOURCETERM : Case not implemented." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
}

Eigen::Vector3d Function::dirichletFunction(double x, double y, double t)
{
  Eigen::Vector3d g(0.,0.,0.);

  // TODO

  return g;
}

Eigen::Vector3d Function::neumannFunction(double x, double y, double t)
{
  Eigen::Vector3d h(0.,0.,0.);

  // TODO
  
  return h;
}
