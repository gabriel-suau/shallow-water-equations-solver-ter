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
  _DF(DF), _mesh(mesh), _function(function), _fluxVector(_mesh->getNumberOfEdges(), 3)
{
}

void FiniteVolume::Initialize(DataFile* DF, Mesh* mesh, Function* function)
{
  _DF = DF;
  _mesh = mesh;
  _function = function;
  _fluxVector.resize(_mesh->getNumberOfEdges(), 3);
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
  _fluxVector.resize(_mesh->getNumberOfEdges(), 3);
}

void LaxFriedrichs::buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol)
{
  // TODO
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
  _fluxVector.resize(_mesh->getNumberOfEdges(), 3);
}

void Rusanov::buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol)
{
  // TODO
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
  _fluxVector.resize(_mesh->getNumberOfEdges(), 3);
}

void HLL::buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol)
{
  // TODO
}
