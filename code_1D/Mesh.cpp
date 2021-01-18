#include "Mesh.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"
#include "DataFile.h"
#include "termcolor.h"
#include <fstream>

Mesh::Mesh()
{
}

Mesh::Mesh(DataFile* DF):
  _DF(DF), _xmin(DF->getXmin()), _xmax(DF->getXmax()), _dx(DF->getDx()), _numberOfCells(DF->getNx()), _cellCenters(_numberOfCells), _cellBoundaries(_numberOfCells + 1)
{
}

void Mesh::Initialize(DataFile* DF)
{
  _DF = DF;
  _xmin = DF->getXmin();
  _xmax = DF->getXmax();
  _dx = DF->getDx();
  _numberOfCells = DF->getNx();
  _cellCenters.resize(_numberOfCells);
  _cellBoundaries.resize(_numberOfCells + 1);
  this->Initialize();
}

void Mesh::Initialize()
{
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Generating a 1D cartesian mesh..." << std::endl;
  for (int i(0) ; i < _numberOfCells ; ++i)
    {
      _cellBoundaries(i) = _xmin + i * _dx;
      _cellCenters(i) = _cellBoundaries(i) + 0.5 * _dx;
    }
  _cellBoundaries(_numberOfCells) = _xmax;
  std::cout << termcolor::green << "SUCCESS::MESH : Mesh generated succesfully !" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}
