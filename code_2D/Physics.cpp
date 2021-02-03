#include "Physics.h"
#include "DataFile.h"
#include "termcolor.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>

Physics::Physics()
{
}

Physics::Physics(DataFile* DF, Mesh* mesh):
  _DF(DF), _mesh(mesh), _g(_DF->getGravityAcceleration()), _nCells(_mesh->getNumberOfCells()), _cellCenters(_mesh->getCellsCenter())
{
}

void Physics::Initialize(DataFile* DF, Mesh* mesh)
{
  // Initialisation
  _DF = DF;
  _mesh = mesh;
  _g = DF->getGravityAcceleration();
  _nCells = _mesh->getNumberOfCells();
  _cellCenters = _mesh->getCellsCenter();
  this->Initialize();
}

void Physics::Initialize()
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
      _Sol0.rightCols(2).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i, 0) = 3.;
        }
    }
  else if (_DF->getScenario() == "RestingLake")
    {
      _Sol0.rightCols(2).setZero();
      double H(3.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i, 0) = std::max(H - _topography(i, 1), 0.);
        }
    }
  else if (_DF->getScenario() == "DamBreak")
    {
      _Sol0.rightCols(2).setZero();
      double Hg(2.), Hd(1.);
      for (int i(0); i < _nCells; i++)
        {
          if (_cellCenters(i, 0) < 0.)
            {
              _Sol0(i, 0) = std::max(Hg, 0.);
            }
          else
            {
              _Sol0(i, 0) = std::max(Hd, 0.);
            }
        }
    }
  else if (_DF->getScenario() == "SinePerturbation")
    {
      _Sol0.rightCols(2).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          // TODO
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

void Physics::buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol)
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

Eigen::Vector3d Physics::dirichletFunction(double x, double y, double t)
{
  Eigen::Vector3d g(0., 0., 0.);
  // TODO
  return g;
}

Eigen::Vector3d Physics::neumannFunction(double x, double y, double t)
{
  Eigen::Vector3d h(0., 0., 0.);
  // TODO
  return h;
}

Eigen::Matrix<double, 3, 2> Physics::physicalFlux(const Eigen::Vector3d& Sol) const
{
  Eigen::Matrix<double, 3, 2> flux;
  double h(Sol(0)), qx(Sol(1)), qy(Sol(2));
  flux(0,0) = qx;
  flux(0,1) = qy;
  flux(1,0) = qx*qx/h + 0.5*_g*h*h;
  flux(1,1) = qx*qy/h;
  flux(2,0) = qx*qy/h;
  flux(2,1) = qy*qy/h + 0.5*_g*h*h;
  return flux;
}

void Physics::computeWaveSpeed(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal, double& lambda1, double& lambda2) const
{
  double hG(SolG(0)), hD(SolD(0));
  Eigen::Vector2d velocityG(SolG(1)/hG, SolG(2)/hG);
  Eigen::Vector2d velocityD(SolD(1)/hD, SolD(2)/hD);
  double normalVelocityG(velocityG.dot(normal));
  double normalVelocityD(velocityD.dot(normal));
  lambda1 = std::min(normalVelocityG - sqrt(_g*hG), normalVelocityD - sqrt(_g*hD));
  lambda2 = std::max(normalVelocityG + sqrt(_g*hG), normalVelocityD + sqrt(_g*hD));
}
