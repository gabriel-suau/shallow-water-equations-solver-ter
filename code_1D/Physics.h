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
  // Pointeur vers le fichier de paramètres pour récupérer les
  // conditions initiales, aux limites et le fichier de topographie (terme source).
  DataFile* _DF;
  Mesh* _mesh;

  // Variables pratiques
  double _xmin, _xmax;
  double _g;
  int _nCells;

  // Variables utiles pour les donnees experimentales
  Eigen::Matrix<double, Eigen::Dynamic, 2> _expBoundaryData;
  int _i;

  // Condition initiale
  Eigen::Matrix<double, Eigen::Dynamic, 2> _Sol0;

  // Topographie pour le terme source.
  Eigen::Matrix<double, Eigen::Dynamic, 2> _fileTopography;
  Eigen::VectorXd _topography;

  // Terme source
  Eigen::Matrix<double, Eigen::Dynamic, 2> _source;

  // Exact solution
  Eigen::Matrix<double, Eigen::Dynamic, 2> _exactSol;
  
public:
  // Constructeur
  Physics();
  Physics(DataFile* DF, Mesh* mesh);

  // Initialisation
  void Initialize();
  void Initialize(DataFile* DF, Mesh* mesh);

  // Getters
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getExperimentalBoundaryData() const {return _expBoundaryData;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getInitialCondition() const {return _Sol0;};
  const Eigen::VectorXd& getTopography() const {return _topography;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getSourceTerm() const {return _source;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getExactSolution() const {return _exactSol;};
  
  // Construit le terme source
  void buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);

  // Construit/Sauvegarde la solution exacte
  void buildExactSolution(double t);
  void saveExactSolution(std::string& fileName) const;
  
  // Conditions aux limites
  Eigen::Vector2d leftBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);
  Eigen::Vector2d rightBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);
  
  // Compute the physical flux of the 1D SWE
  Eigen::Vector2d physicalFlux(const Eigen::Vector2d& Sol) const;
  // Compute the eigenvalues of the flux jacobian
  void computeWaveSpeed(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD, double* lambda1, double* lambda2) const;
  
protected:
  void buildTopography();
  void buildInitialCondition();
  void buildExpBoundaryData();

  // Boundary conditions

  // Resolution equation second ordre
  double FindRacine(double a, double b, double c);
  // On cherche le terme source en x (pour x dans le domaine)
  double FindSourceX(double x);

  // Exact solution
  
  // Stationnary flows over a bump test cases.
  void computeCoeffabcd(double qIn, double hOut, double z, double zEnd, double* a, double* b, double* c, double *d);
  double cardanP(double a, double b, double c) const;
  double cardanQ(double a, double b, double c, double d) const;
  double cardanDet(double p, double q) const;
  double exactHeight(double p, double q, double a, double b, double hnear, double hMax) const;
  double RHJump(double hplus, double hminus, double q) const;

  // Dam breaks test cases
  double damBreakFunc(double x, double vG, double vD) const;
};

#endif // PHYSICS_H
