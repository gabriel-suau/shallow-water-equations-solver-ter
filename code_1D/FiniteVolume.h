#ifndef FINITE_VOLUME_H
#define FINITE_VOLUME_H

#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

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
  Eigen::Matrix<double, Eigen::Dynamic, 2> _fluxVector;
  
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
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getFluxVector() const {return _fluxVector;};
  
  // Build the flux vector
  virtual Eigen::Vector2d numFlux(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD) const = 0;
  void buildFluxVector(const double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);
};

class LaxFriedrichs: public FiniteVolume
{
public:
  // Constructeur
  LaxFriedrichs();
  LaxFriedrichs(DataFile* DF, Mesh* mesh, Physics* physics);

  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Build flux vector
  Eigen::Vector2d numFlux(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD) const;
};

class Rusanov: public FiniteVolume
{
public:
  // Constructeur
  Rusanov();
  Rusanov(DataFile* DF, Mesh* mesh, Physics* physics);

  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Build flux vector
  Eigen::Vector2d numFlux(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD) const;
};

class HLL: public FiniteVolume
{
public:
  // Constructeur
  HLL();
  HLL(DataFile* DF, Mesh* mesh, Physics* physics);

  // Initialisation
  void Initialize(DataFile* DF, Mesh* mesh, Physics* physics);

  // Build flux vector
  Eigen::Vector2d numFlux(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD) const;
};

#endif //FINITE_VOLUME_H
