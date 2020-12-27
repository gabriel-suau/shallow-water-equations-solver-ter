#include "FiniteVolume.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Function.h"

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

FiniteVolume::FiniteVolume(DataFile* DF, Mesh* mesh, Function* function):
  _DF(DF), _mesh(mesh), _function(function), _fluxVector(_mesh->getNumberOfCells() + 1, 2)
{
}

void FiniteVolume::Initialize(DataFile* DF, Mesh* mesh, Function* function)
{
  _DF = DF;
  _mesh = mesh;
  _function = function;
  _fluxVector.resize(_mesh->getNumberOfCells() + 1, 2);
}

//---------------------------------------------------//
//---------------Flux de LaxFriedrichs---------------//
//---------------------------------------------------//
LaxFriedrichs::LaxFriedrichs():
  FiniteVolume()
{
}

LaxFriedrichs::LaxFriedrichs(DataFile* DF, Mesh* mesh, Function* function):
  FiniteVolume(DF, mesh, function)
{
  _fluxName = "LF";
}

void LaxFriedrichs::Initialize(DataFile* DF, Mesh* mesh, Function* function)
{
  _DF = DF;
  _mesh = mesh;
  _function = function;
  _fluxName = "LF";
  _fluxVector.resize(_mesh->getNumberOfCells() + 1, 2);
}

void LaxFriedrichs::buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  // Récupère les données importantes
  double g(_DF->getGravityAcceleration());
  double dx(_DF->getDx()), dt(_DF->getTimeStep());
  double b(dx/dt);
  int nCells(Sol.rows());

  for (int i(1) ; i < nCells ; ++i)
    {
      // Récupère les données importantes
      double hg(Sol.row(i-1)(0));
      double qg(Sol.row(i-1)(1));
      double hd(Sol.row(i)(0));
      double qd(Sol.row(i)(1));
      // Construit le flux
      if (hg != 0. && hd != 0.)
        {
          _fluxVector.row(i)(0) = 0.5 * (qd + qg - b * (hd - hg));
          _fluxVector.row(i)(1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) + pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (qd - qg));
        }
      else if (hg == 0. && hd != 0.)
        {
          _fluxVector.row(i)(0) = 0.5 * (qd - b * hd);
          _fluxVector.row(i)(1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) - b * qd);
        }
      else if (hd == 0 && hg != 0)
        {
          _fluxVector.row(i)(0) = 0.5 * (qg + b * hg);
          _fluxVector.row(i)(1) = 0.5 * (pow(qg,2)/hg + 0.5*g*pow(hg,2) + b * qg);
        }
      else
        {
          _fluxVector.row(i) << 0., 0.;
        }
    }
  _fluxVector.row(0) << 0., 0.;
  _fluxVector.row(nCells) << 0., 0.;
}

//---------------------------------------------//
//---------------Flux de Rusanov---------------//
//---------------------------------------------//
Rusanov::Rusanov():
  FiniteVolume()
{
}

Rusanov::Rusanov(DataFile* DF, Mesh* mesh, Function* function):
  FiniteVolume(DF, mesh, function)
{
  _fluxName = "Rusanov";
}

void Rusanov::Initialize(DataFile* DF, Mesh* mesh, Function* function)
{
  _DF = DF;
  _mesh = mesh;
  _function = function;
  _fluxName = "Rusanov";
  _fluxVector.resize(_mesh->getNumberOfCells() + 1, 2);
}

void Rusanov::buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  // Récupère les données importantes
  double g(_DF->getGravityAcceleration());
  int nCells(Sol.rows());

  for (int i(1) ; i < nCells ; ++i)
    {
      // Récupère les données importantes
      double hg(Sol.row(i-1)(0));
      double qg(Sol.row(i-1)(1));
      double hd(Sol.row(i)(0));
      double qd(Sol.row(i)(1));
      // Valeurs propres
      double lambda1(abs(qg/hg) + sqrt(g * hg));
      double lambda2(abs(qd/hd) + sqrt(g * hd));
      double b;
      // Construit le flux
      if (hg != 0. && hd != 0.)
        {
          b = std::max(lambda1,lambda2);
          _fluxVector.row(i)(0) = 0.5 * (qd + qg - b * (hd - hg));
          _fluxVector.row(i)(1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) + pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (qd - qg));
        }
      else if (hg == 0. && hd != 0.)
        {
          b = lambda2;
          _fluxVector.row(i)(0) = 0.5 * (qd - b * hd);
          _fluxVector.row(i)(1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) - b * qd);
        }
      else if (hd == 0 && hg != 0)
        {
          b = lambda1;
          _fluxVector.row(i)(0) = 0.5 * (qg + b * hg);
          _fluxVector.row(i)(1) = 0.5 * (pow(qg,2)/hg + 0.5*g*pow(hg,2) + b * qg);
        }
      else
        {
          _fluxVector.row(i) << 0., 0.;
        }
    }
  _fluxVector.row(0) << 0., 0.;
  _fluxVector.row(nCells) << 0., 0.;
}
