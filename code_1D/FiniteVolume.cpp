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

  // Boucle sur les arêtes intérieures
  for (int i(1) ; i < nCells ; ++i)
    {
      // Récupère les données importantes
      double hg(Sol(i-1,0));
      double qg(Sol(i-1,1));
      double hd(Sol(i,0));
      double qd(Sol(i,1));
      // Construit le flux
      if (hg != 0. && hd != 0.)
        {
          _fluxVector(i,0) = 0.5 * (qd + qg - b * (hd - hg));
          _fluxVector(i,1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) + pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (qd - qg));
        }
      else if (hg == 0. && hd != 0.)
        {
          _fluxVector(i,0) = 0.5 * (qd + qg - b * hd);
          _fluxVector(i,1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) - b * (qd - qg));
        }
      else if (hd == 0 && hg != 0)
        {
          _fluxVector(i,0) = 0.5 * (qd + qg + b * hg);
          _fluxVector(i,1) = 0.5 * (pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (qd - qg));
        }
      else
        {
          _fluxVector.row(i) << 0., 0.;
        }
    }

  // Flux aux limites (à déterminer avec les CL)
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

  // Boucle sur les arêtes intérieures
  for (int i(1) ; i < nCells ; ++i)
    {
      // Récupère les données importantes
      double hg(Sol(i-1,0));
      double qg(Sol(i-1,1));
      double hd(Sol(i,0));
      double qd(Sol(i,1));
      // Valeurs propres (pas tout à fait, mais bon...)
      double lambda1(abs(qg/hg) + sqrt(g * hg));
      double lambda2(abs(qd/hd) + sqrt(g * hd));
      double b;
      // Construit le flux
      if (hg != 0. && hd != 0.)
        {
          b = std::max(lambda1,lambda2);
          _fluxVector(i,0) = 0.5 * (qd + qg - b * (hd - hg));
          _fluxVector(i,1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) + pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (qd - qg));
        }
      else if (hg < 1e-6 && hd != 0.)
        {
          b = lambda2;
          _fluxVector(i,0) = 0.5 * (qd - b * hd);
          _fluxVector(i,1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) - b * (qd));
        }
      else if (hd < 1e-6 && hg != 0)
        {
          b = lambda1;
          _fluxVector(i,0) = 0.5 * (qg + b * hg);
          _fluxVector(i,1) = 0.5 * (pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (- qg));
        }
      else
        {
          _fluxVector.row(i) << 0., 0.;
        }
    }

  // Flux en entrée (à déterminer avec les CL)
  double hg(_function->dirichletFunction(_DF->getXmin(), 0.)(0));
  double qg(_function->dirichletFunction(_DF->getXmin(), 0.)(1) * hg);
  double hd(Sol(0,0));
  double qd(Sol(0,1));
  // Valeurs propres (pas tout à fait, mais bon...)
  double lambda1(abs(qg/hg) + sqrt(g * hg));
  double lambda2(abs(qd/hd) + sqrt(g * hd));
  double b;
  
  // Construit le flux
  if (hg != 0. && hd != 0.)
    {
      b = std::max(lambda1,lambda2);
      _fluxVector(0,0) = 0.5 * (qd + qg - b * (hd - hg));
      _fluxVector(0,1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) + pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (qd - qg));
    }
  else if (hg < 1e-6 && hd != 0.)
    {
      b = lambda2;
      _fluxVector(0,0) = 0.5 * (qd - b * hd);
      _fluxVector(0,1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) - b * (qd));
    }
  else if (hd < 1e-6 && hg != 0)
    {
      b = lambda1;
      _fluxVector(0,0) = 0.5 * (qg + b * hg);
      _fluxVector(0,1) = 0.5 * (pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (- qg));
    }
  else
    {
      _fluxVector.row(0) << 0., 0.;
    }

  // Flux en sortie (à déterminer avec les CL)
  hg = Sol(nCells-1,0);
  qg = Sol(nCells-1,1);
  hd = Sol(nCells-1,0);
  qd = Sol(nCells-1,1);
  // Valeurs propres (pas tout à fait, mais bon...)
  lambda1 = abs(qg/hg) + sqrt(g * hg);
  lambda2 = abs(qd/hd) + sqrt(g * hd);
  
  // Construit le flux
  if (hg != 0. && hd != 0.)
    {
      b = std::max(lambda1,lambda2);
      _fluxVector(nCells,0) = 0.5 * (qd + qg - b * (hd - hg));
      _fluxVector(nCells,1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) + pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (qd - qg));
    }
  else if (hg < 1e-6 && hd != 0.)
    {
      b = lambda2;
      _fluxVector(nCells,0) = 0.5 * (qd - b * hd);
      _fluxVector(nCells,1) = 0.5 * (pow(qd,2)/hd + 0.5*g*pow(hd,2) - b * (qd));
    }
  else if (hd < 1e-6 && hg != 0)
    {
      b = lambda1;
      _fluxVector(nCells,0) = 0.5 * (qg + b * hg);
      _fluxVector(nCells,1) = 0.5 * (pow(qg,2)/hg + 0.5*g*pow(hg,2) - b * (- qg));
    }
  else
    {
      _fluxVector.row(nCells) << 0., 0.;
    }
}

//--------------------------------------//
//---------------Flux HLL---------------//
//--------------------------------------//
HLL::HLL():
  FiniteVolume()
{
}

HLL::HLL(DataFile* DF, Mesh* mesh, Function* function):
  FiniteVolume(DF, mesh, function)
{
  _fluxName = "HLL";
}

void HLL::Initialize(DataFile* DF, Mesh* mesh, Function* function)
{
  _DF = DF;
  _mesh = mesh;
  _function = function;
  _fluxName = "HLL";
  _fluxVector.resize(_mesh->getNumberOfCells() + 1, 2);
}

void HLL::buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  // Récupère les données importantes
  double g(_DF->getGravityAcceleration());
  int nCells(Sol.rows());

  // Boucle sur les arêtes intérieures
  for (int i(1) ; i < nCells ; ++i)
    {
      // Récupère les données importantes
      double hg(Sol(i-1,0));
      double qg(Sol(i-1,1));
      double hd(Sol(i,0));
      double qd(Sol(i,1));
      // Valeurs propres (cette fois c'est bien les valeurs propres)
      double lambda1(qg/hg + sqrt(g * hg));
      double lambda2(qg/hg - sqrt(g * hg));
      double lambda3(qd/hd + sqrt(g * hd));
      double lambda4(qd/hd - sqrt(g * hd));
      double c1, c2;
      // Construit le flux
      if (hg != 0. && hd != 0.)
        {
          c1 = std::min(lambda2,lambda4);
          c2 = std::max(lambda1,lambda3);
          if (c1 >= 0.)
            {
              _fluxVector(i,0) = qg;
              _fluxVector(i,1) = qg*qg/hg + 0.5*g*hg*hg;
            }
          else if (c1 < 0. && 0. < c2)
            {
              _fluxVector(i,0) = (c2 * qg - c1 * qd + c1*c2 * (hd - hg))/(c2 - c1);
              _fluxVector(i,1) = (c2 * (qg*qg/hg + 0.5*g*hg*hg) - c1 * (qd*qd/hd + 0.5*g*hd*hd) + c1*c2*(qd - qg))/(c2 - c1);
            }
          else if (c2 <= 0.)
            {
              _fluxVector(i,0) = qd;
              _fluxVector(i,1) = qd*qd/hd + 0.5*g*hd*hd;              
            }
        }
      else if (hg < 1e-6 && hd != 0.)
        {
          c1 = lambda4;
          c2 = lambda3;
          if (c1 >= 0.)
            {
              _fluxVector(i,0) = 0.;
              _fluxVector(i,1) = 0.;
            }
          else if (c1 < 0. && 0. < c2)
            {
              _fluxVector(i,0) = (- c1 * qd + c1*c2 * hd)/(c2 - c1);
              _fluxVector(i,1) = (- c1 * (qd*qd/hd + 0.5*g*hd*hd) + c1*c2*qd)/(c2 - c1);
            }
          else if (c2 <= 0.)
            {
              _fluxVector(i,0) = qd;
              _fluxVector(i,1) = qd*qd/hd + 0.5*g*hd*hd;              
            }
        }
      else if (hd < 1e-6 && hg != 0)
        {
          c1 = lambda2;
          c2 = lambda1;
          if (c1 >= 0.)
            {
              _fluxVector(i,0) = qg;
              _fluxVector(i,1) = qg*qg/hg + 0.5*g*hg*hg;
            }
          else if (c1 < 0. && 0. < c2)
            {
              _fluxVector(i,0) = (c2 * qg - c1*c2 * hg)/(c2 - c1);
              _fluxVector(i,1) = (c2 * (qg*qg/hg + 0.5*g*hg*hg) - c1*c2*qg)/(c2 - c1);
            }
          else if (c2 <= 0.)
            {
              _fluxVector(i,0) = 0.;
              _fluxVector(i,1) = 0.;
            }
        }
      else
        {
          _fluxVector.row(i) << 0., 0.;
        }
    }

  // Flux aux limites (à déterminer avec les CL)
  _fluxVector.row(0) << 0., 0.;
  _fluxVector.row(nCells) << 0., 0.;
}
