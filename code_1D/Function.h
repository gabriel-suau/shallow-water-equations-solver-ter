#ifndef FUNCTION_H
#define FUNCTION_H

#include "DataFile.h"
#include "Mesh.h"
#include "termcolor.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

class Function
{
private:
  // Pointeur vers le fichier de paramètres pour récupérer les
  // conditions initiales, aux limites et le fichier de topographie (terme source).
  DataFile* _DF;
  Mesh* _mesh;

  // Variables pratiques
  double _xmin, _xmax;
  double _g;
  int _nCells;
  Eigen::VectorXd _cellCenters;

  // Variables utiles pour les donnees experimentales
  Eigen::Matrix<double, Eigen::Dynamic, 2> _expBoundaryData;
  int _i;

  // Condition initiale
  Eigen::Matrix<double, Eigen::Dynamic, 2> _Sol0;

  // Topographie (x,z) pour le terme source.
  Eigen::Matrix<double, Eigen::Dynamic, 2> _topography;
  Eigen::Matrix<double, Eigen::Dynamic, 2> _source;

  // Exact solution
  Eigen::Matrix<double, Eigen::Dynamic, 2> _exactSol;
  
public:
  // Constructeur
  Function();
  Function(DataFile* DF, Mesh* mesh);

  // Initialisation
  void Initialize();
  void Initialize(DataFile* DF, Mesh* mesh);

  // Getters
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getExperimentalBoundaryData() const {return _expBoundaryData;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getInitialCondition() const {return _Sol0;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getTopography() const {return _topography;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getSourceTerm() const {return _source;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getExactSolution() const {return _exactSol;};
  
  // Construit le terme source
  void buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);

  // Construit la solution exacte
  void buildExactSolution(double t);
  
  // Conditions aux limites
  Eigen::Vector2d leftBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);
  Eigen::Vector2d rightBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);
  
  // Compute the physical flux of the 1D SWE
  Eigen::Vector2d physicalFlux(const Eigen::Vector2d& Sol) const;
  // Compute the eigenvalues of the flux jacobian
  void computeWaveSpeed(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD, double* lambda1, double* lambda2) const;

  // Compute the exact solution of the case (if one exists)
  Eigen::Vector2d exactSolution(double x, double t) const;
  
protected:
  void buildTopography();
  void buildInitialCondition();
  void buildExpBoundaryData();
  // Resolution equation second ordre
  double FindRacine(double a, double b, double c);
  // On cherche le terme source en x (pour x dans le domaine)
  double FindSourceX(double x);
};

#endif // FUNCTION_H
