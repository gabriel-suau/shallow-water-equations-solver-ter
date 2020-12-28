#ifndef MESH_H
#define MESH_H

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"
#include "DataFile.h"
#include <fstream>

class Mesh
{
private:
  // Pointeur vers le fichier de paramètres
  DataFile* _DF;

  // Quantités importantes
  double _xmin;
  double _xmax;
  double _dx;
  int _numberOfCells;
  Eigen::VectorXd _cellCenters; // numberOfCells |-|-|-|-|-|
  Eigen::VectorXd _cellBoundaries; // numberOfCells+1

public:
  // Constructeurs
  Mesh();
  Mesh(DataFile* DF);

  // Destructeur
  ~Mesh() = default;

  // Initialisation
  void Initialize(DataFile* DF);
  void Initialize();

  // Getters
  const Eigen::VectorXd& getCellCenters() const {return _cellCenters;};
  const Eigen::VectorXd& getCellBoundaries() const {return _cellBoundaries;};
  int getNumberOfCells() const {return _numberOfCells;};
  double getSpaceStep() const {return _dx;};
  double getxMin() const {return _xmin;};
  double getxMax() const {return _xmax;};
};


#endif // MESH_H
