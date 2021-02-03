#include "FiniteVolume.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <cmath>

//--------------------------------------------------//
//--------------------Base Class--------------------//
//--------------------------------------------------//
FiniteVolume::FiniteVolume()
{
}

FiniteVolume::FiniteVolume(DataFile* DF, Mesh* mesh, Physics* physics):
  _DF(DF), _mesh(mesh), _physics(physics), _fluxVector(_mesh->getNumberOfCells(), 3)
{
}

void FiniteVolume::Initialize(DataFile* DF, Mesh* mesh, Physics* physics)
{
  _DF = DF;
  _mesh = mesh;
  _physics = physics;
  _fluxVector.resize(_mesh->getNumberOfCells(), 3);
}


//--------------------------------------------------//
//-------------------Rusanov flux-------------------//
//--------------------------------------------------//
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
  _fluxVector.resize(_mesh->getNumberOfCells(), 3);
}

// Compute the numerical flux across an edge
Eigen::Vector3d Rusanov::numFlux1D(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal) const
{
  // Vecteur flux au travers de l'arête.
  Eigen::Vector3d flux;

  // Calcul de b
  double lambda1, lambda2;
  _physics->computeWaveSpeed(SolG, SolD, normal, lambda1, lambda2);
  double b(std::max(lambda1,lambda2));

  // Calcul du flux
  flux = 0.5 * ((_physics->physicalFlux(SolD) + _physics->physicalFlux(SolG))*normal - b * (SolD - SolG));
  return flux;
}


void Rusanov::buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol)
{
  // Reset the flux 
  _fluxVector.setZero();

  // Get mesh parameters
  // Edges
  int nbEdges(_mesh->getNumberOfEdges());
  const std::vector<Edge>& edges(_mesh->getEdges());
  const Eigen::VectorXd& edgesLength(_mesh->getEdgesLength());
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& edgesNormal(_mesh->getEdgesNormal());

  // Boucle sur les arêtes
  for (int i(0) ; i < nbEdges ; ++i)
    {
      int c1(edges[i].getC1()), c2(edges[i].getC2());
      double edgeLength(edgesLength(i));
      Eigen::Vector2d edgeNormal(edgesNormal.row(i));
      // Boundary edges
      if (c2 == -1)
        {
          Eigen::Vector3d flux1D(numFlux1D(Sol.row(c1), Sol.row(c1), edgeNormal));
          _fluxVector.row(c1) += edgeLength * flux1D;
        }
      // Interior edges
      else
        {
          Eigen::Vector3d flux1D(numFlux1D(Sol.row(c1), Sol.row(c2), edgeNormal));
          _fluxVector.row(c1) += edgeLength * flux1D;
          _fluxVector.row(c2) -= edgeLength * flux1D;
        }
    }
}

//--------------------------------------------------//
//---------------------HLL flux---------------------//
//--------------------------------------------------//
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
  _fluxVector.resize(_mesh->getNumberOfCells(), 3);
}

// Compute the numerical flux across an edge
Eigen::Vector3d HLL::numFlux1D(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal) const
{
  Eigen::Vector3d flux;

  // TODO
}

void HLL::buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol)
{
  // TODO
}
