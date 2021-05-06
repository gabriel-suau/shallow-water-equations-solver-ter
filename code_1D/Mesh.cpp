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
  _DF(DF), _xmin(DF->getXmin()), _xmax(DF->getXmax()), _dx(DF->getDx()), _numberOfCells(DF->getNx()), _cellCenters(_numberOfCells)
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
  this->Initialize();
}

void Mesh::Initialize()
{
#if VERBOSITY>0
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Generating a 1D cartesian mesh..." << std::endl;
#endif
  for (int i(0) ; i < _numberOfCells ; ++i)
    {
      _cellCenters(i) = _xmin + (i + 0.5) * _dx;
    }
#if VERBOSITY>0
  std::cout << termcolor::green << "SUCCESS::MESH : Mesh generated succesfully !" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
#endif
}
