#include "Function.h"
#include "DataFile.h"
#include "termcolor.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <algorithm>


//------------------------------------------//
//---------------Constructors---------------//
//------------------------------------------//
Function::Function()
{
}



Function::Function(DataFile* DF, Mesh* mesh):
  _DF(DF), _mesh(mesh), _xmin(mesh->getxMin()), _xmax(mesh->getxMax()), _g(_DF->getGravityAcceleration()), _nCells(mesh->getNumberOfCells()), _cellCenters(mesh->getCellCenters()), _i(0)
{
}



//--------------------------------------------//
//---------------Initialization---------------//
//--------------------------------------------//
void Function::Initialize(DataFile* DF, Mesh* mesh)
{
  _DF = DF;
  _mesh = mesh;
  _xmin = mesh->getxMin();
  _xmax = mesh->getxMax();
  _g = DF->getGravityAcceleration();
  _i = 0;
  _nCells = mesh->getNumberOfCells();
  _cellCenters = mesh->getCellCenters();
  this->Initialize();
}



void Function::Initialize()
{
  // Logs
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Building topography, initial condition, and experimental data..." << std::endl;

  // Build
  _source.resize(_nCells, 2);
  buildTopography();
  buildInitialCondition();
  if (_DF->getLeftBC() == "DataFile" || _DF->getRightBC() == "DataFile")
    buildExpBoundaryData();

  // Logs
  std::cout << termcolor::green << "SUCCESS::FUNCTION : Everything was successfully built." << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}



//----------------------------------------------//
//---------------Build Topography---------------//
//----------------------------------------------//
void Function::buildTopography()
{
  _topography.resize(_nCells, 2);
  _topography.col(0) = _cellCenters;
  // Flat botttom
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _topography.col(1).setZero();
    }
  // Thacker test case topography
  else if (_DF->getTopographyType() == "Thacker")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i,1) = 0.05 * (_cellCenters(i) - _xmin);
        }
    }
  // Bump topography
  else if (_DF->getTopographyType() == "Bump")
    {
      _topography.col(1).setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(_cellCenters(i));
          if (8 < x && x < 12)
            _topography(i,1) = 0.2 - 0.05 * pow(x - 10, 2);
        }
    }
 // Read the topography in a data file
  else if (_DF->getTopographyType() == "File")
    {
      const std::string topoFile(_DF->getTopographyFile());
      std::ifstream topoStream(topoFile);
      std::string line, properLine;
      double dummy;
      int i(0);
      if (!topoStream.is_open())
        {
          std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Unable to open the topography file : " << topoFile << std::endl;
          std::cout << termcolor::reset << "====================================================================================================" << std::endl;
          exit(-1);
        }
      else
        {
          std::cout << "Building the topography from file : " << topoFile << std::endl;
        }
      while(getline(topoStream, line))
        {
          properLine = regex_replace(line, std::regex(",") , std::string(" "));
          std::stringstream ss(properLine);
          ss >> dummy >> _topography(i,1);
          ++i;
        }
    }
  else
    {
      std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Case not implemented" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
  std::cout << termcolor::green << "SUCCESS::TOPOGRAPHY : Topography was successfully built." << std::endl;
  std::cout << termcolor::reset;
}



//-----------------------------------------------------//
//---------------Build Initial Condition---------------//
//-----------------------------------------------------//
void Function::buildInitialCondition()
{
  _Sol0.resize(_nCells, 2);
  if (_DF->getInitialCondition() == "UniformHeightAndDischarge")
    {
      double H0(_DF->getInitialHeight()), q0(_DF->getInitialDischarge());
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = std::max(H0 - _topography(i,1), 0.);
          _Sol0(i,1) = q0;
        }
    }
  else if (_DF->getInitialCondition() == "DamBreakWet")
    {
      _Sol0.col(1).setZero();
      double Hl(2.), Hr(1.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (_cellCenters(i) < 0.5*(_xmax + _xmin))
            {
              _Sol0(i,0) = std::max(Hl - _topography(i,1), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(Hr - _topography(i,1), 0.);
            }
        }
    }
  else if (_DF->getInitialCondition() == "DamBreakDry")
    {
      _Sol0.col(1).setZero();
      double Hl(2.), Hr(0.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (_cellCenters(i) < 0.5*(_xmax + _xmin))
            {
              _Sol0(i,0) = std::max(Hl - _topography(i,1), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(Hr - _topography(i,1), 0.);
            }
        }
    }
  else if (_DF->getInitialCondition() == "SinePerturbation")
    {
      _Sol0.col(1).setZero();
      double H(2.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(_cellCenters(i));
          if (-1 < x && x < 1)
            {
              _Sol0(i,0) = std::max(H + 0.2 * cos(M_PI * x) - _topography(i,1), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(1.8 - _topography(i,1),0.);
            }
        }
    }
  else
    {
      std::cout << termcolor::red << "ERROR::INITIALCONDITION : Case not implemented" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
  std::cout << termcolor::green << "SUCCESS::INITIALCONDITION : Initial Condition was successfully built." << std::endl;
  std::cout << termcolor::reset;
}



//--------------------------------------------------------------//
//---------------Build Experimental Boundary Data---------------//
//--------------------------------------------------------------//
void Function::buildExpBoundaryData()
{  
  const std::string expDataFile(_DF->getLeftBCDataFile());
  std::ifstream expDataStream(expDataFile);
  std::string line, properLine;
  int i(0);
  if (!expDataStream.is_open())
    {
      std::cout << termcolor::red << "ERROR::EXPDATA : Unable to open the experimental data file : " << expDataFile << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
  else
    {
      std::cout << "Building the experimental data from file : " << expDataFile << std::endl;
    }
  
  // First line indicates the number of experimental values
  getline(expDataStream, line);
  int size(0);
  std::stringstream ss(line);
  ss >> size;
  _expBoundaryData.resize(size,2);
  
  // Then we read the data
  while(getline(expDataStream, line))
    {
      properLine = regex_replace(line, std::regex(",") , std::string(" "));
      std::stringstream ss(properLine);
      ss >> _expBoundaryData(i,0) >> _expBoundaryData(i,1);
      ++i;
    }
  std::cout << termcolor::green << "SUCCESS::EXPDATA : Experimental data was successsfully built." << std::endl;
  std::cout << termcolor::reset;
}


//-----------------------------------------------//
//---------------Build Source Term---------------//
//-----------------------------------------------//
void Function::buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  // Construit le terme source en fonction de la topographie.
  _source.setZero();
  // Flat bottom
  if (_DF->getTopographyType() == "FlatBottom")
    {
      // Do nothing
    }
  // Bump topography
  else if (_DF->getTopographyType() == "Bump")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(_cellCenters(i));
          if (8 < x  && x < 12)
            _source(i,1) = _g * Sol(i,0) * 0.05 * 2. * (x - 10.);
        }
    }
  // Thacker test case topography
  else if (_DF->getTopographyType() == "Thacker")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * 0.05;
        }
    }
  // Topography file
  else if (_DF->getTopographyType() == "File")
    {
      double dx(_mesh->getSpaceStep());
      _source(0,1) = - _g * Sol(0,0) * (-_topography(2) + 4.*_topography(1) - 3.*_topography(0))/(2.*dx);
      for (int i(1) ; i < _nCells - 1 ; ++i)
        {
          _source(i,1) = - _g * Sol(i,0) * (_topography(i+1) - _topography(i-1))/(2. * dx);
        }
      _source(_nCells - 1, 1) = - _g * Sol(_nCells - 1,0) * (3.*_topography(_nCells - 1) - 4.*_topography(_nCells - 2) + _topography(_nCells - 3))/(2.*dx);
    }
  // Not implemented
  else
    {
      std::cout << termcolor::red << "ERROR::SOURCETERM : Case not implemented." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
}



//-------------------------------------------//
//---------------Physical flux---------------//
//-------------------------------------------//
Eigen::Vector2d Function::physicalFlux(const Eigen::Vector2d& Sol) const
{
  Eigen::Vector2d flux;
  double h(Sol(0)), qx(Sol(1));
  flux(0) = qx;
  flux(1) = qx*qx/h + 0.5*_g*h*h;
  return flux;
}


//----------------------------------------//
//---------------Wave speed---------------//
//----------------------------------------//
void Function::computeWaveSpeed(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD, double* lambda1, double* lambda2) const
{
  double hG(SolG(0)), hD(SolD(0));
  double uG(SolG(1)/hG), uD(SolD(1)/hD);
  *lambda1 = std::min(uG - sqrt(_g * hG), uD - sqrt(_g * hD));
  *lambda2 = std::max(uG + sqrt(_g * hG), uD + sqrt(_g * hD));
}


//--------------------------------------------//
//---------------Exact Solution---------------//
//--------------------------------------------//
// Compute the exact solutions corresponding to the test cases
Eigen::Vector2d Function::exactSolution(double x, double t) const
{
  Eigen::Vector2d exactSol;
  if (_DF->getInitialCondition() == "UniformHeightAndDischarge")
    {
      exactSol(0) = _DF->getInitialHeight();
      exactSol(1) = _DF->getInitialDischarge();
    }
  else if (_DF->getInitialCondition() == "DamBreakWet")
    {
      // TODO
    }
  else if (_DF->getInitialCondition() == "DamBreakDry")
    {
      // TODO
    }
  return exactSol;
}



//------------------------------------------------------//
//---------------Left Boundary Conditions---------------//
//------------------------------------------------------//
Eigen::Vector2d Function::leftBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  Eigen::Vector2d SolG(0.,0.);

  // Calcul du nombre de Froude au bord
  double h(Sol(0,0)), q(Sol(0,1));
  double Fr(abs(q)/(h * sqrt(_g * h)));
  
  // Choix entre les differentes CL
  if (_DF->getLeftBC() == "Neumann")
    {
      SolG(0) = Sol(0,0);
      SolG(1) = Sol(0,1);
    }
  else if (_DF->getLeftBC() == "Wall")
    {
      SolG(0) = Sol(0,0);
      SolG(1) = 0;
    }
  else if (_DF->getLeftBC() == "ImposedConstantDischarge")
    {
      // Entrée/sortie fluviale
      if (Fr < 1)
        {
          SolG(0) = Sol(0,0);
          SolG(1) = _DF->getLeftBCImposedDischarge();
        }
      // Sortie torrentielle (sortie libre, on n'impose rien)
      else if (Fr > 1 && q < 0)
        {
          SolG(0) = Sol(0,0);
          SolG(1) = Sol(0,1);
        }
      // Entrée torrentielle (on impose une hauteur et un debit)
      else if (Fr > 1 && q > 0)
        {
          SolG(0) = _DF->getLeftBCImposedHeight();
          SolG(1) = _DF->getLeftBCImposedDischarge();
        }
    }
  else if (_DF->getLeftBC() == "PeriodicWaves" || _DF->getLeftBC() == "DataFile" || _DF->getLeftBC() == "ImposedConstantHeight")
    {
      // Recupere la solution dans les mailles de centre x1 et x2 ainsi que dx et dt
      double h1(Sol(0,0)), h2(Sol(1,0));
      double u1(Sol(0,1)/h1), u2(Sol(1,1)/h2);
      double dx(_DF->getDx()), dt(_DF->getTimeStep());
      double x1(_DF->getXmin() + 0.5*dx);
      double a(pow(1 + dt/dx * (u2 - u1), 2));
      double b(2*dt*(u1 - x1/dx * (u2 - u1)) * (1 + dt/dx * (u2 - u1)) - dt*dt*_g*(h2 - h1)/dx);
      double c(pow(dt*u1 - dt/dx * x1 * (u2 - u1), 2) - dt*dt * _g * (h1 - x1/dx * (h2 - h1)));
      double xe(FindRacine(a, b, c));
      double uXe(u1 + (xe - x1)*(u2 - u1)/dx);
      double hXe(h1 + (xe - x1)*(h2 - h1)/dx);
      double source_terme_xe(FindSourceX(xe));
      double beta_moins_xe_tn(uXe - 2*sqrt(_g*hXe));
      double beta_moins_0_tnplus1(beta_moins_xe_tn - _g*dt*source_terme_xe);
      if (_DF->getLeftBC() == "ImposedConstantHeight")
        {
          // Entrée/sortie fluviale
          if (Fr < 1)
            {
              SolG(0) = _DF->getLeftBCImposedHeight();
              SolG(1) = SolG(0)*(beta_moins_0_tnplus1 + 2*sqrt(_g*SolG(0)));
            }
          // Sortie torrentielle (sortie libre, on n'impose rien)
          else if (Fr > 1 && q < 0)
            {
              SolG(0) = Sol(0,0);
              SolG(1) = Sol(0,1);
            }
          // Entrée torrentielle (on impose une hauteur et un debit)
          else if (Fr > 1 && q > 0)
            {
              SolG(0) = _DF->getLeftBCImposedHeight();
              SolG(1) = _DF->getLeftBCImposedDischarge();
            }
        }
      if (_DF->getLeftBC() == "PeriodicWaves")
        {
          SolG(0) = 3. + 0.1*sin(5 * M_PI * t);
          SolG(1) = SolG(0)*(beta_moins_0_tnplus1 + 2*sqrt(_g*SolG(0)));
        }
      else if (_DF->getLeftBC() == "DataFile")
        {
          int i_max;
          i_max = _expBoundaryData.rows() - _i;
          double temps1(_expBoundaryData(_i,0)), temps2(_expBoundaryData(_i+1,0));
          double hauteur1(_expBoundaryData(_i,1)), hauteur2(_expBoundaryData(_i+1,1));
          while ((!(temps1 < t && t <= temps2)) && (_i < i_max)) // On cherche a trouver les temps connus des capteurs tq temps1 < t <= temps2
            {
              _i++;
              temps1 = _expBoundaryData(_i,0);
              temps2 = _expBoundaryData(_i+1,0);
            }
          if (_i == i_max)
            {
              std::cout << "Aucun pas de temps ne correspond --> pas trop eleve ?" << std::endl;
              return SolG;
            }
          hauteur1 = _expBoundaryData(_i,1);
          hauteur2 = _expBoundaryData(_i+1,1);
          SolG(0) = hauteur1 + (t - temps1)*(hauteur2 - hauteur1)/(temps2 - temps1);
          SolG(1) = SolG(0)*(beta_moins_0_tnplus1 + 2*sqrt(_g*SolG(0)));
        }
    }
  return SolG;
}



//-------------------------------------------------------//
//---------------Right Boundary Conditions---------------//
//-------------------------------------------------------//
Eigen::Vector2d Function::rightBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  Eigen::Vector2d SolD(0.,0.);

  // Calcul du nombre de Froude au bord
  double h(Sol(_nCells - 1,0)), q(Sol(_nCells - 1,1));
  double Fr(abs(q)/(h * sqrt(_g * h)));
  
  // Choix entre les differentes CL
  if (_DF->getRightBC() == "Neumann")
    {
      SolD(0) = Sol(_nCells - 1,0);
      SolD(1) = Sol(_nCells - 1,1);
    }
  else if (_DF->getRightBC() == "Wall")
    {
      SolD(0) = Sol(_nCells - 1,0);
      SolD(1) = 0;
    }
  else if (_DF->getRightBC() == "ImposedConstantDischarge")
    {
      // Entrée/sortie fluviale
      if (Fr < 1)
        {
          SolD(0) = Sol(_nCells - 1,0);
          SolD(1) = _DF->getRightBCImposedDischarge();
        }
      // Sortie torrentielle (sortie libre, on n'impose rien)
      else if (Fr > 1 && q > 0)
        {
          SolD(0) = Sol(_nCells - 1,0);
          SolD(1) = Sol(_nCells - 1,1);
        }
      // Entrée torrentielle (on impose une hauteur et un debit)
      else if (Fr > 1 && q < 0)
        {
          SolD(0) = _DF->getRightBCImposedHeight();
          SolD(1) = _DF->getRightBCImposedDischarge();
        }
    }
  else if (_DF->getRightBC() == "PeriodicWaves" || _DF->getRightBC() == "DataFile" || _DF->getRightBC() == "ImposedConstantHeight")
    {
      // Recupere la solution dans la maille de bord
      double h1(Sol(_nCells - 1,0)), u1(Sol(_nCells - 1,1)/h1);
      if (_DF->getRightBC() == "ImposedConstantHeight")
        {
          // Entrée/sortie fluviale
          if (Fr < 1)
            {
              SolD(0) = _DF->getRightBCImposedHeight();
              SolD(1) = SolD(0) * (u1 + 2. * sqrt(_g * h1) - 2. * sqrt(_g * SolD(0)));
            }
          // Sortie torrentielle (sortie libre, on n'impose rien)
          else if (Fr > 1 && q > 0)
            {
              SolD(0) = Sol(_nCells - 1,0);
              SolD(1) = Sol(_nCells - 1,1);
            }
          // Entrée torrentielle (on impose une hauteur et un debit)
          else if (Fr > 1 && q < 0)
            {
              SolD(0) = _DF->getRightBCImposedHeight();
              SolD(1) = _DF->getRightBCImposedDischarge();
            }
        }
      else if (_DF->getRightBC() == "PeriodicWaves")
        {
          SolD(0) = 3. + 0.1*sin(5 * M_PI * t);
          SolD(1) = SolD(0) * (u1 + 2. * sqrt(_g * h1) - 2. * sqrt(_g * SolD(0)));
        }
      else if (_DF->getRightBC() == "DataFile")
        {
          int i_max;
          i_max = _expBoundaryData.size() - _i;
          double temps1(_expBoundaryData(_i,0)), temps2(_expBoundaryData(_i+1,0));
          double hauteur1(_expBoundaryData(_i,1)), hauteur2(_expBoundaryData(_i+1,1));
          while ((!(temps1 < t <= temps2)) && (_i < i_max)) // On cherche a trouver les temps connus des capteurs tq temps1 < t <= temps2
            {
              _i++;
              temps1 = _expBoundaryData(_i,0);
              temps2 = _expBoundaryData(_i+1,0);
            }
          if (_i == i_max)
            {
              std::cout << "Aucun pas de temps ne correspond --> pas trop eleve ?" << std::endl;
              return SolD;
            }
          hauteur1 = _expBoundaryData(_i,1);
          hauteur2 = _expBoundaryData(_i+1,1);
          SolD(0) = (t - temps1)*(hauteur2 - hauteur1)/(temps2 - temps1) - hauteur1;
          SolD(1) = SolD(0) * (u1 + 2. * sqrt(_g * h1) - 2. * sqrt(_g * SolD(0)));
        }
    }
  return SolD;
}



// Other
double Function::FindRacine(double a, double b, double c)
{
  double delta;
  delta = b*b - 4*a*c;
  if (delta < 0)
    {
      std::cout << "Pas de racine reelle" << std::endl;
      exit(1);
    }
  else if (delta == 0)
    {
      double r(-b/(2*a));
      return r;
    }
  else
    {
      double r1, r2;
      r1 = (-b - sqrt(delta))/(2*a);
      r2 = (-b + sqrt(delta))/(2*a);
      return r2;
    }
}



// Donne le terme source en x par interpolation
double Function::FindSourceX(double x)
{
  int i(0);
  double dx(_DF->getDx());
  double x1(_xmin + (i + 0.5) * dx), x2(x1 + dx);
  double source1, source2, source;
  while (x < x1)
    {
      i++;
      x1 += dx;
    }
  x2 = x1 + dx;
  source1 = _source(i,1);
  source2 = _source(i+1,1);
  source = source1 + (x - x1)*(source2 - source1)/(x2 - x1);
  return source;
}
