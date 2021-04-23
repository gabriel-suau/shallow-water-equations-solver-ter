#include "FiniteVolume.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <cmath>



//--------------------------------------------------------//
//---------------Classe mère flux numérique---------------//
//--------------------------------------------------------//
FiniteVolume::FiniteVolume()
{
}



FiniteVolume::FiniteVolume(DataFile* DF, Mesh* mesh, Physics* physics):
  _DF(DF), _mesh(mesh), _physics(physics), _fluxVector(_mesh->getNumberOfCells(), 2)
{
}



void FiniteVolume::Initialize(DataFile* DF, Mesh* mesh, Physics* physics)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _fluxVector.resize(_mesh->getNumberOfCells(), 2);
}


void FiniteVolume::buildFluxVector(const double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  // Reset the flux
  _fluxVector.setZero();

  // Get mesh parameters
  int nCells(_mesh->getNumberOfCells());
  // Loop on the edges
  for (int i(0) ; i <= nCells ; ++i)
    {
      Eigen::Vector2d SolG, SolD;
      // Left Boundary
      if (i == 0)
        {
          SolG = _physics->leftBoundaryFunction(t + _DF->getTimeStep(), Sol);
          SolD = Sol.row(0);
          _fluxVector.row(0) += numFlux(SolG, SolD);
        }
      // Right Boundary
      else if (i == nCells)
        {
          SolG = Sol.row(nCells - 1);
          SolD = _physics->rightBoundaryFunction(t + _DF->getTimeStep(), Sol);
          _fluxVector.row(nCells - 1) -= numFlux(SolG, SolD);
        }
      // Interior
      else
        {
          SolG = Sol.row(i-1);
          SolD = Sol.row(i);
          _fluxVector.row(i-1) -= numFlux(SolG, SolD);
          _fluxVector.row(i) += numFlux(SolG, SolD);
        }
    }
}


//---------------------------------------------------//
//---------------Flux de LaxFriedrichs---------------//
//---------------------------------------------------//
LaxFriedrichs::LaxFriedrichs():
  FiniteVolume()
{
}



LaxFriedrichs::LaxFriedrichs(DataFile* DF, Mesh* mesh, Physics* function):
  FiniteVolume(DF, mesh, function)
{
  _fluxName = "LF";
}



void LaxFriedrichs::Initialize(DataFile* DF, Mesh* mesh, Physics* physics)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _fluxName = "LF";
  _fluxVector.resize(_mesh->getNumberOfCells(), 2);
}


Eigen::Vector2d LaxFriedrichs::numFlux(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD) const
{
  // Vecteur flux au travers d'une arete
  Eigen::Vector2d flux;
  
  // Recupere dt et dx
  double dt(_DF->getTimeStep()), dx(_DF->getDx());
  double b(dx/dt);

  // Calcul du flux
  flux = 0.5 * ((_physics->physicalFlux(SolD) + _physics->physicalFlux(SolG)) - b * (SolD - SolG));
  
  return flux;
}


//---------------------------------------------//
//---------------Flux de Rusanov---------------//
//---------------------------------------------//
Rusanov::Rusanov():
  FiniteVolume()
{
}



Rusanov::Rusanov(DataFile* DF, Mesh* mesh, Physics* physics):
  FiniteVolume(DF, mesh, physics)
{
  _fluxName = "Rusanov";
}



void Rusanov::Initialize(DataFile* DF, Mesh* mesh, Physics* physics)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _fluxName = "Rusanov";
  _fluxVector.resize(_mesh->getNumberOfCells(), 2);
}


Eigen::Vector2d Rusanov::numFlux(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD) const
{
  // Vecteur flux au travers d'une arete
  Eigen::Vector2d flux;
  
  // Calcul de b
  double lambda1, lambda2;
  _physics->computeWaveSpeed(SolG, SolD, &lambda1, &lambda2);
  double b(std::max(abs(lambda1),abs(lambda2)));

  // Calcul du flux
  flux = 0.5 * ((_physics->physicalFlux(SolD) + _physics->physicalFlux(SolG)) - b * (SolD - SolG));
  
  return flux;
}


//--------------------------------------//
//---------------Flux HLL---------------//
//--------------------------------------//
HLL::HLL():
  FiniteVolume()
{
}



HLL::HLL(DataFile* DF, Mesh* mesh, Physics* physics):
  FiniteVolume(DF, mesh, physics)
{
  _fluxName = "HLL";
}



void HLL::Initialize(DataFile* DF, Mesh* mesh, Physics* physics)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _fluxName = "HLL";
  _fluxVector.resize(_mesh->getNumberOfCells(), 2);
}


Eigen::Vector2d HLL::numFlux(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD) const
{
  // Vecteur flux au travers d'une arete
  Eigen::Vector2d flux;
  
  // Calcul de b
  double lambda1, lambda2;
  _physics->computeWaveSpeed(SolG, SolD, &lambda1, &lambda2);

  // Calcul du flux
  if (0 <= lambda1)
    {
      flux = _physics->physicalFlux(SolG);
    }
  else if (lambda1 < 0 && 0 < lambda2)
    {
      flux = (lambda2 * _physics->physicalFlux(SolG) - lambda1 * _physics->physicalFlux(SolD) + lambda2 * lambda1 * (SolD - SolG))/(lambda2 - lambda1);
    }
  else if (lambda2 <= 0)
    {
      flux = _physics->physicalFlux(SolD);
    }
  
  return flux;
}
