#ifndef PHYSICS_H
#define PHYSICS_H

#include "DataFile.h"
#include "Mesh.h"
#include "termcolor.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

class Physics
{
private:
  // Pointer to the DataFile (to get the initial condition,
  // the boundary conditions and the topography).
  DataFile* _DF;
  // Pointer to a Mesh object
  Mesh* _mesh;

  // Useful variables
  double _g;
  int _nCells;
  Eigen::Matrix<double, Eigen::Dynamic, 2> _cellCenters;

  // Initial condition
  Eigen::Matrix<double, Eigen::Dynamic, 3> _Sol0;
  
  // Topography and source term
  Eigen::VectorXd _topography;
  Eigen::Matrix<double, Eigen::Dynamic, 3> _source;
  
public:
  // Constructeur
  Physics();
  Physics(DataFile* DF, Mesh* mesh);

  // Initialisation
  void Initialize();
  void Initialize(DataFile* DF, Mesh* mesh);

  // Getters
  const Eigen::Matrix<double, Eigen::Dynamic, 3>& getInitialCondition() const {return _Sol0;};
  const Eigen::VectorXd& getTopography() const {return _topography;};
  const Eigen::Matrix<double, Eigen::Dynamic, 3>& getSourceTerm() const {return _source;};
  
  // Construit le terme source
  void buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 3>& Sol);

  // Conditions aux limites
  Eigen::Vector3d dirichletFunction(double x, double y, double t);
  Eigen::Vector3d neumannFunction(double x, double y, double t);

  // Compute the physical flux
  Eigen::Matrix<double, 3, 2> physicalFlux(const Eigen::Vector3d& Sol) const;

  // Compute the eigenvalues of the flux jacobian
  void computeWaveSpeed(const Eigen::Vector3d& SolG, const Eigen::Vector3d& SolD, const Eigen::Vector2d& normal, double& lambda1, double& lambda2) const;
};

#endif // PHYSICS_H
