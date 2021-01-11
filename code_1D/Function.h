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
  int _nCells;
  Eigen::VectorXd _cellCenters;
  
  // Condition initiale
  Eigen::Matrix<double, Eigen::Dynamic, 2> _Sol0;
  
  // Topography (x,z) pour le terme source.
  Eigen::Matrix<double, Eigen::Dynamic, 2> _topography;
  Eigen::Matrix<double, Eigen::Dynamic, 2> _source;
  
public:
  // Constructeur
  Function();
  Function(DataFile* DF, Mesh* mesh);

  // Initialisation
  void Initialize();
  void Initialize(DataFile* DF, Mesh* mesh);

  // Getters
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getInitialCondition() const {return _Sol0;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getTopography() const {return _topography;};
  const Eigen::Matrix<double, Eigen::Dynamic, 2>& getSourceTerm() const {return _source;};

  // Build le terme source
  void buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol);

  // Conditions aux limites
  Eigen::Vector2d dirichletFunction(double x, double t);
  Eigen::Vector2d neumannFunction(double x, double t);
};

#endif // FUNCTION_H
