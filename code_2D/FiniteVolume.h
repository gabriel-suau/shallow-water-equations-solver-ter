#ifndef FINITE_VOLUME_H
#define FINITE_VOLUME_H

#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

//--------------------------------------------------//
//--------------------Base Class--------------------//
//--------------------------------------------------//
class FiniteVolume
{
protected:
  // Pointeurs vers les trucs importants
  DataFile* _DF;
  Mesh* _mesh;
  Physics* _physics;

  // Nom du flux numérique
  std::string _fluxName;

  // Vecteur des flux
  Eigen::Matrix<double, Eigen::Dynamic, 3> _fluxVector;
  
public:
  // Constructeurs
  FiniteVolume();
  FiniteVolume(DataFile* DF, Mesh* mesh, Physics* physics);

  // Destructeur
  virtual ~FiniteVolume() = default;
  
  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Getters
  const std::string& getFluxName() const {return _fluxName;};
  const Eigen::Matrix<double, Eigen::Dynamic, 3>& getFluxVector() const {return _fluxVector;};
  
  // Fluxes
  virtual Eigen::Vector3d numFlux1D(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal) const = 0;
  virtual void buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol) = 0;
};


//--------------------------------------------------//
//-------------------Rusanov flux-------------------//
//--------------------------------------------------//
class Rusanov: public FiniteVolume
{
public:
  // Constructeur
  Rusanov();
  Rusanov(DataFile* DF, Mesh* mesh, Physics* physics);

  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Build flux vector
  void buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol);
  Eigen::Vector3d numFlux1D(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal) const;
};


//--------------------------------------------------//
//---------------------HLL flux---------------------//
//--------------------------------------------------//
class HLL: public FiniteVolume
{
public:
  // Constructeur
  HLL();
  HLL(DataFile* DF, Mesh* mesh, Physics* physics);

  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Build flux vector
  void buildFluxVector(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol);
  Eigen::Vector3d numFlux1D(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal) const;
};

#endif //FINITE_VOLUME_H
