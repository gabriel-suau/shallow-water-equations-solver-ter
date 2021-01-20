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
  double _g;
  int _nCells;
  Eigen::Matrix<double, Eigen::Dynamic, 2> _cellCenters;

  // Condition initiale
  Eigen::Matrix<double, Eigen::Dynamic, 3> _Sol0;
  
  // Topographie pour le terme source.
  Eigen::VectorXd _topography;
  Eigen::Matrix<double, Eigen::Dynamic, 3> _source;

public:
  // Constructeur
  Function();
  Function(DataFile* DF, Mesh* mesh);

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
};

#endif // FUNCTION_H
