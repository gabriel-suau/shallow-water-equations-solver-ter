#include "FiniteVolume.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <cmath>
#include <algorithm>


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
  double dx(_mesh->getSpaceStep());

  // Get gravity
  double g(_DF->getGravityAcceleration());

  // Vectors to store the reconstruted values at the left and right of each interface
  Eigen::Matrix<double, Eigen::Dynamic, 2> SolD, SolG;
  SolD.resize(nCells + 1, 2);
  SolG.resize(nCells + 1, 2);

  // Select order of the scheme
  switch(_DF->getSchemeOrder())
    {
      // First order, the reconstructed values are the cell-centered approximations
    case 1:
      // Left boundary
      SolG.row(0) = _physics->leftBoundaryFunction(t + _DF->getTimeStep(), Sol);
      SolD.row(0) = Sol.row(0);
      // Right boundary
      SolG.row(nCells) = Sol.row(nCells - 1);
      SolD.row(nCells) = _physics->rightBoundaryFunction(t + _DF->getTimeStep(), Sol);
      // Interior edges
      for (int i(1) ; i < nCells ; ++i)
        {
          SolG.row(i) = Sol.row(i-1);
          SolD.row(i) = Sol.row(i);
        }
      break;
      
      // Second Order MUSCL, the reconstructed values are obtained via linear interpolation
      // + slope limitation (minmod limiter) to get a TVD scheme.
    case 2:
      // Vector to store the slopes and the limited slopes for the piecewise linear reconstruction
      Eigen::Matrix<double, Eigen::Dynamic, 2> slopes, limSlopes;
      slopes.resize(nCells + 1, 2);
      limSlopes.resize(nCells, 2);
      
      // Compute the slopes
      // Left boundary
      Eigen::Vector2d leftBoundarySol(_physics->leftBoundaryFunction(t + _DF->getTimeStep(), Sol));
      slopes(0,0) = (Sol(0,0) - leftBoundarySol(0)) / dx;
      slopes(0,1) = (Sol(0,1) - leftBoundarySol(1)) / dx;
      // Right boundary
      Eigen::Vector2d rightBoundarySol(_physics->rightBoundaryFunction(t + _DF->getTimeStep(), Sol));
      slopes(nCells, 0) = (rightBoundarySol(0) - Sol(nCells - 1, 0)) / dx;
      slopes(nCells, 1) = (rightBoundarySol(1) - Sol(nCells - 1, 1)) / dx;
      // Interior edges
      for (int i(1) ; i < nCells ; ++i)
        {
          slopes.row(i) = (Sol.row(i) - Sol.row(i-1)) / dx;
        }

      // Limit the slopes
      for (int i(0) ; i < nCells - 1 ; ++i)
        {
          limSlopes(i,0) = minmod(slopes(i,0), slopes(i+1,0));
          limSlopes(i,1) = minmod(slopes(i,1), slopes(i+1,1));
        }

      // Reconstruct the values at each edge
      // Left boundary
      SolG.row(0) = leftBoundarySol;
      SolD.row(0) = Sol.row(0) - 0.5 * dx * limSlopes.row(0);
      // Right boundary
      SolG.row(nCells) = Sol.row(nCells - 1) + 0.5 * dx * limSlopes.row(nCells - 1);
      SolD.row(nCells) = rightBoundarySol;
      // Interior edges
      for (int i(1) ; i < nCells ; ++i)
        {
          SolG.row(i) = Sol.row(i-1) + 0.5 * dx * limSlopes.row(i-1);
          SolD.row(i) = Sol.row(i) - 0.5 * dx * limSlopes.row(i);
        }
      break;
    }

  
  // // Select the well-balanced reconstruction
  // // Hydrostatic
  // Eigen::Matrix<double, Eigen::Dynamic, 2> tempG, tempD;
  // tempG.resize(nCells + 1, 2);
  // tempD.resize(nCells + 1, 2);
  // const Eigen::VectorXd& topography(_physics->getTopography());
  // // Left boundary
  // tempG(0,0) = std::max(0., SolG(0,0));
  // tempD(0,0) = std::max(0., SolD(0,0));
  // tempG(0,1) = tempG(0,0) * SolG(0,1) / SolG(0,0);
  // tempD(0,1) = tempD(0,0) * SolD(0,1) / SolD(0,0);
  // // Interior edges
  // for (int i(1) ; i < nCells ; ++i)
  //   {
  //     double zi(topography(i-1)), zip(topography(i));
  //     tempG(i,0) = std::max(0., SolG(i,0) + zi - std::max(zi, zip));
  //     tempD(i,0) = std::max(0., SolD(i,0) + zip - std::max(zi, zip));
  //     tempG(i,1) = tempG(i,0) * SolG(i,1)/SolG(i,0);
  //     tempD(i,1) = tempD(i,0) * SolD(i,1)/SolD(i,0);
  //   }
  // // Right boundary
  // tempG(nCells,0) = std::max(0., SolG(nCells,0));
  // tempD(nCells,0) = std::max(0., SolD(nCells,0));
  // tempG(nCells,1) = tempG(nCells,0) * SolG(nCells,1) / SolG(nCells,0);
  // tempD(nCells,1) = tempD(nCells,0) * SolD(nCells,1) / SolD(nCells,0);

  // // Source term reconstruction
  // Eigen::Matrix<double, Eigen::Dynamic, 2> sourceG, sourceD;
  // sourceG.resize(nCells + 1, 2);
  // sourceD.resize(nCells + 1, 2);
  // for (int i(0) ; i <= nCells ; ++i)
  //   {
  //     sourceG.row(i) << 0. , 0.5 * g * (pow(SolG(i,0), 2) - pow(tempG(i,0), 2));
  //     sourceD.row(i) << 0. , 0.5 * g * (pow(SolD(i,0), 2) - pow(tempD(i,0), 2));
  //   }

  // // oui
  // SolG = tempG;
  // SolD = tempD;
  
  // Build the flux vector using the reconstructed values at each edge
  // Left boundary contribution
  _fluxVector.row(0) += numFlux(SolG.row(0), SolD.row(0));
  // Interior fluxes contribution
  for (int i(1) ; i < nCells; ++i)
    {
      Eigen::Vector2d flux(numFlux(SolG.row(i), SolD.row(i)));
      _fluxVector.row(i-1) -= flux;
      _fluxVector.row(i) += flux;
    }
  // Right boundary contribution
  _fluxVector.row(nCells - 1) -= numFlux(SolG.row(nCells), SolD.row(nCells));
}



// Minmod slope limiter
double FiniteVolume::minmod(double a, double b) const
{
  if (a * b < 0)
    return 0.;
  else if (abs(a) < abs(b))
    return a;
  else
    return b;
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
  double hg(SolG(0));
  double hd(SolD(0));
  if (hg > 1e-6 && hd > 1e-6)
    {
      flux = 0.5 * ((_physics->physicalFlux(SolD) + _physics->physicalFlux(SolG)) - b * (SolD - SolG)); 
    }
  else if (hg < 1e-6 && hd > 1e-6)
    {
      flux = 0.5 * (_physics->physicalFlux(SolD) - b * SolD);
    }
  else if (hd < 1e-6 && hg > 1e-6)
    {
      flux = 0.5 * (_physics->physicalFlux(SolG) + b * SolG);
    }
  else
    {
      flux << 0. , 0.; 
    }
  
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
  double hg(SolG(0));
  double hd(SolD(0));
  if (0 <= lambda1)
    {
      if (hg < 1e-6)
        flux << 0. , 0.;
      else
        flux = _physics->physicalFlux(SolG);
    }
  else if (lambda1 < 0 && 0 < lambda2)
    {
      if (hg > 1e-6 && hd > 1e-6)
        flux = (lambda2 * _physics->physicalFlux(SolG) - lambda1 * _physics->physicalFlux(SolD) + lambda2 * lambda1 * (SolD - SolG))/(lambda2 - lambda1);
      else if (hg < 1e-6 && hd > 1e-6)
        flux = (- lambda1 * _physics->physicalFlux(SolD) + lambda2 * lambda1 * (SolD))/(lambda2 - lambda1);
      else if (hd < 1e-6 && hg > 1e-6)
        flux = (lambda2 * _physics->physicalFlux(SolG) + lambda2 * lambda1 * (-SolG))/(lambda2 - lambda1);
      else
        flux << 0. , 0.;
    }
  else if (lambda2 <= 0)
    {
      if (hd < 1e-6)
        flux << 0. , 0.;
      else
        flux = _physics->physicalFlux(SolD);
    }
  
  return flux;
}
