#include "Physics.h"
#include "DataFile.h"
#include "termcolor.h"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <algorithm>
#include <complex>

//------------------------------------------//
//---------------Constructors---------------//
//------------------------------------------//
Physics::Physics()
{
}



Physics::Physics(DataFile* DF, Mesh* mesh):
  _DF(DF), _mesh(mesh), _xmin(mesh->getxMin()), _xmax(mesh->getxMax()), _g(_DF->getGravityAcceleration()), _nCells(mesh->getNumberOfCells()), _i(0)
{
}



//--------------------------------------------//
//---------------Initialization---------------//
//--------------------------------------------//
void Physics::Initialize(DataFile* DF, Mesh* mesh)
{
  _DF = DF;
  _mesh = mesh;
  _xmin = mesh->getxMin();
  _xmax = mesh->getxMax();
  _g = DF->getGravityAcceleration();
  _i = 0;
  _nCells = mesh->getNumberOfCells();
  this->Initialize();
}



void Physics::Initialize()
{
  // Logs
#if VERBOSITY>0
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Building topography, initial condition, and experimental data..." << std::endl;
#endif

  // Build
  _source.resize(_nCells, 2);
  buildTopography();
  buildInitialCondition();
  if (_DF->getLeftBC() == "DataFile" || _DF->getRightBC() == "DataFile")
    buildExpBoundaryData();

  // Logs
#if VERBOSITY>0
  std::cout << termcolor::green << "SUCCESS::FUNCTION : Everything was successfully built." << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
#endif
}



//----------------------------------------------//
//---------------Build Topography---------------//
//----------------------------------------------//
void Physics::buildTopography()
{
  _topography.resize(_nCells);
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
  
  // Flat botttom
  if (_DF->getTopographyType() == "FlatBottom")
    {
      _topography.setZero();
    }
  // Thacker test case topography
  else if (_DF->getTopographyType() == "Thacker")
    {
      _topography.setZero();
      double xmin(_DF->getXmin()), xmax(_DF->getXmax()), L(xmax - xmin);
      double a(1.), h0(0.5);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          _topography(i) = h0 * (1. / pow(a,2) * pow(x - 0.5 * L, 2) - 1.);
        }
    }
  // Bump topography
  else if (_DF->getTopographyType() == "Bump")
    {
      _topography.setZero();
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          if (8 < x && x < 12)
            _topography(i) = 0.2 - 0.05 * pow(x - 10, 2);
        }
    }
 // Read the topography in a data file
  else if (_DF->getTopographyType() == "File")
    {
      const std::string topoFile(_DF->getTopographyFile());
      std::ifstream topoStream(topoFile);
      std::string line, properLine;
      double size;
      int i(0);
      if (!topoStream.is_open())
        {
          std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Unable to open the topography file : " << topoFile << std::endl;
          std::cout << termcolor::reset << "====================================================================================================" << std::endl;
          exit(-1);
        }
#if VERBOSITY>0
      else
        {
          std::cout << "Building the topography from file : " << topoFile << std::endl;
        }
#endif
      topoStream >> size;
      _fileTopography.resize(size + 1, 2);
      while(getline(topoStream, line))
        {
          properLine = regex_replace(line, std::regex(",") , std::string(" "));
          std::stringstream ss(properLine);
          ss >> _fileTopography(i,0) >> _fileTopography(i,1);
          ++i;
        }

      // Ajuste la topographie au domaine (interpolation lineaire)
      double x1(_fileTopography(0,0)), z1(_fileTopography(0,1));
      double x2(_fileTopography(1,0)), z2(_fileTopography(1,1));
      int j(0);
      for (int k(0) ; k < _nCells ; ++k)
        {
          double x(cellCenters(k));
          while (!(x1 < x && x <= x2))
            {
              ++j;
              x1 = _fileTopography(j,0);
              x2 = _fileTopography(j+1,0);
            }
          z1 = _fileTopography(j,1);
          z2 = _fileTopography(j+1,1);
          _topography(k) = z1 + (x - x1) * (z2 - z1) / (x2 - x1);
        }

      // Fixe la plus petite valeur de la topographie a 0.
      double topoMin(_topography.minCoeff());
      for (int i(0) ; i < _nCells ; ++i)
        {
          _topography(i) += abs(topoMin);
        }
    }
  else
    {
      std::cout << termcolor::red << "ERROR::TOPOGRAPHY : Case not implemented" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
#if VERBOSITY>0
  std::cout << termcolor::green << "SUCCESS::TOPOGRAPHY : Topography was successfully built." << std::endl;
  std::cout << termcolor::reset;
#endif
}



//-----------------------------------------------------//
//---------------Build Initial Condition---------------//
//-----------------------------------------------------//
void Physics::buildInitialCondition()
{
  _Sol0.resize(_nCells, 2);
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
  if (_DF->getInitialCondition() == "UniformHeightAndDischarge")
    {
      double H0(_DF->getInitialHeight()), q0(_DF->getInitialDischarge());
      for (int i(0) ; i < _nCells ; ++i)
        {
          _Sol0(i,0) = std::max(H0 - _topography(i), 0.);
          _Sol0(i,1) = q0;
        }
    }
  else if (_DF->getInitialCondition() == "DamBreakWet")
    {
      _Sol0.col(1).setZero();
      double Hl(2.), Hr(1.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (cellCenters(i) < 0.5*(_xmax + _xmin))
            {
              _Sol0(i,0) = std::max(Hl - _topography(i), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(Hr - _topography(i), 0.);
            }
        }
    }
  else if (_DF->getInitialCondition() == "DamBreakDry")
    {
      _Sol0.col(1).setZero();
      double Hl(2.), Hr(0.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          if (cellCenters(i) < 0.5*(_xmax + _xmin))
            {
              _Sol0(i,0) = std::max(Hl - _topography(i), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(Hr - _topography(i), 0.);
            }
        }
    }
  else if (_DF->getInitialCondition() == "Thacker")
    {
      _Sol0.col(1).setZero();
      double xmin(_DF->getXmin()), xmax(_DF->getXmax()), L(xmax - xmin);
      double a(1.), h0(0.5);
      double x1(- 0.5 - a + 0.5 * L), x2(- 0.5 + a + 0.5 * L);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          if (x1 <= x && x <= x2)
            {
              _Sol0(i,0) = - h0 * (pow(1. / a * (x - 0.5 * L) + 1. / (2. * a), 2) - 1);
            }
          else
            _Sol0(i,0) = 0.;
        }
    }
  else if (_DF->getInitialCondition() == "SinePerturbation")
    {
      _Sol0.col(1).setZero();
      double H(2.);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          if (-1 < x && x < 1)
            {
              _Sol0(i,0) = std::max(H + 0.2 * cos(M_PI * x) - _topography(i), 0.);
            }
          else
            {
              _Sol0(i,0) = std::max(1.8 - _topography(i),0.);
            }
        }
    }
  else if (_DF->getInitialCondition() == "File")
    {
      const std::string initFile(_DF->getInitFile());
      std::ifstream initStream(initFile);
      std::string line, trash;
      int i(0);
      // First line is comments.
      getline(initStream, line);
      // Now we build the initial condition
      while(getline(initStream, line))
        {
          std::stringstream ss(line);
          ss >> trash >> trash >> _Sol0(i,0) >> trash >> _Sol0(i,1) >> trash ;
          ++i;
        }
    }
  else
    {
      std::cout << termcolor::red << "ERROR::INITIALCONDITION : Case not implemented" << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
#if VERBOSITY>0
  std::cout << termcolor::green << "SUCCESS::INITIALCONDITION : Initial Condition was successfully built." << std::endl;
  std::cout << termcolor::reset;
#endif
}



//--------------------------------------------------------------//
//---------------Build Experimental Boundary Data---------------//
//--------------------------------------------------------------//
void Physics::buildExpBoundaryData()
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
#if VERBOSITY>0
  else
    {
      std::cout << "Building the experimental data from file : " << expDataFile << std::endl;
    }
#endif
  
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
#if VERBOSITY>0
  std::cout << termcolor::green << "SUCCESS::EXPDATA : Experimental data was successsfully built." << std::endl;
  std::cout << termcolor::reset;
#endif
}


//-----------------------------------------------//
//---------------Build Source Term---------------//
//-----------------------------------------------//
void Physics::buildSourceTerm(const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
{
  // Construit le terme source en fonction de la topographie.
  _source.setZero();
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
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
          double x(cellCenters(i));
          if (8 < x  && x < 12)
            _source(i,1) = _g * Sol(i,0) * 0.05 * 2. * (x - 10.);
        }
    }
  // Thacker test case topography
  else if (_DF->getTopographyType() == "Thacker")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          double xmin(_DF->getXmin()), xmax(_DF->getXmax()), L(xmax - xmin);
          double a(1.), h0(0.5);
          double x(cellCenters(i));
          _source(i,1) = - _g * Sol(i,0) * h0 * (2. / pow(a,2) * (x - 0.5 * L));
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



//--------------------------------------------------//
//---------------Build Exact Solution---------------//
//--------------------------------------------------//
void Physics::buildExactSolution(double t)
{
  _exactSol.resize(_nCells, 2);
  const std::string& testCase(_DF->getTestCase());
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
  // Resting lake solutions
  if (testCase == "RestingLake")
    {
      for (int i(0) ; i < _nCells ; ++i)
        {
          _exactSol(i,0) = std::max(_DF->getLeftBCImposedHeight() - _topography(i), 0.);
          _exactSol(i,1) = 0.;
        }
    }
  else if (testCase == "DamBreakWet" || testCase == "DamBreakDry")
    {
      // Mesh parameters
      double xmin(_DF->getXmin()), xmax(_DF->getXmax());
      double xdam(0.5 * (xmax - xmin));
      // Parameters for the dichotomy
      double eps(1e-6);
      int nmax(1000);
      int iter(0);
      // Heights
      double hG(0.), hD(0.), hMid(0.);
      // Middle discharge and velocity
      double qMid(0.), uMid(0.);
      // Values for the dichotomy
      double xA(0.), xB(0.), func(0.), mid(0.);
      // Velocity of the shock wave
      double v(0.);
      // Dam break on a wet domain
      if (testCase == "DamBreakWet")
        {
          hG = 2.; hD = 1.;
        }
      // Dam break on a dry domain
      else if (testCase == "DamBreakDry")
        {
          hG = 2.; hD = 0.;
        }
      // Wave velocities.
      double cG(sqrt(_g * hG)), cD(sqrt(_g * hD)), cMid(0.);
      func = damBreakFunc(cG, cG, cD);
      if (func < 0.)
        xA = cG;
      else
        xB = cG;
      func = damBreakFunc(cD, cG, cD);
      if (func < 0.)
        xA = cD;
      else
        xB = cD;
      // Dichotomy
      while (abs(xA - xB) > eps && iter < nmax)
        {
          ++iter;
          mid = 0.5 * (xA + xB);
          func = damBreakFunc(mid, cG, cD);
          if (func < 0.)
            xA = mid;
          else
            xB = mid;
        }
      // cMid
      cMid = 0.5 * (xA + xB);
      if (hD == 0.)
        cMid = 0.;
      hMid = cMid * cMid / _g;
      uMid = 2. * (cG - cMid);
      qMid = hMid * uMid;
      v = qMid / (hMid - hD);
      // Compute the solution depending on where we are.
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          if (x <= xdam - cG * t)
            {
              _exactSol(i,0) = hG;
              _exactSol(i,1) = 0.;
            }
          else if (x <= xdam + (2. * cG - 3. * cMid) * t)
            {
              _exactSol(i,0) = 4. / (9. * _g) * pow(sqrt(_g * hG) - 0.5 * (x - xdam) / t, 2);
              _exactSol(i,1) = _exactSol(i,0) * 2. / 3. * ((x - xdam) / t + sqrt(_g * hG));
            }
          else if (x <= xdam + v * t)
            {
              _exactSol(i,0) = hMid;
              _exactSol(i,1) = qMid;
            }
          else
            {
              _exactSol(i,0) = hD;
              _exactSol(i,1) = 0.;
            }
        }
    }
  // Thacker test case
  else if (testCase == "Thacker")
    {
      double xmax(_DF->getXmax()), xmin(_DF->getXmin());
      double L(xmax - xmin);
      double a(1.);
      double h0(0.5);
      double x1(- 0.5 * cos(sqrt(2. * _g * h0) * t / a) - a + 0.5 * L);
      double x2(- 0.5 * cos(sqrt(2. * _g * h0) * t / a) + a + 0.5 * L);
      for (int i(0) ; i < _nCells ; ++i)
        {
          double x(cellCenters(i));
          if (x1 <= x && x <= x2)
            {
              _exactSol(i,0) = - h0 * (pow(1. / a * (x - 0.5 * L) + 1. / (2. * a) * cos(sqrt(2. * _g * h0) * t / a), 2) - 1);
              _exactSol(i,1) = sqrt(2. * _g * h0) / (2. * a) * sin(sqrt(2. * _g * h0) * t / a);
            }
          else
            {
              _exactSol(i,0) = 0.;
              _exactSol(i,1) = 0.;
            }
        }
    }
  // Non hydrostatic stationnary solutions
  else if (testCase == "SubcriticalFlow" || testCase == "TranscriticalFlowWithoutShock" || testCase == "TranscriticalFlowWithShock")
    {
      double epsilon(1.0 / _nCells);
      double qIn(_DF->getLeftBCImposedDischarge());
      double hOut(_DF->getRightBCImposedHeight());
      double hMiddle(pow(pow(qIn,2)/_g, 1./3.)); // = critical height
      double hMax(3.), zMax(0.2);
      double a(0.), b(0.), c(0.), d(0.);
      double p, q;
      // Discharge must be constant in the whole domain
      for (int i(0) ; i < _nCells ; ++i)
        {
          _exactSol(i,1) = qIn;
        }
      // Subcritical flow
      if (testCase == "SubcriticalFlow")
        {
          for (int i(_nCells - 1) ; i >= 0  ; --i)
            {
              double z(_topography(i)), zEnd(_topography(_nCells - 1));
              computeCoeffabcd(qIn, hOut, z, zEnd, &a, &b, &c, &d);
              p = cardanP(a, b, c);
              q = cardanQ(a, b, c, d);
              double hnear;
              if (i == _nCells - 1)
                hnear = hOut;
              else
                hnear = _exactSol(i+1, 0);
              _exactSol(i,0) = exactHeight(p, q, a, b, hnear, hMax);
            }
        }
      else if (testCase == "TranscriticalFlowWithoutShock")
        {
          // Subcritical part (before the bump)
          for (int i(2. * _nCells / 5. - 1) ; i >= 0 ; --i)
            {
              double z(_topography(i));
              computeCoeffabcd(qIn, hMiddle, z, zMax, &a, &b, &c, &d);
              p = cardanP(a, b, c);
              q = cardanQ(a, b, c, d);
              double hnear;
              if (i == 2. * _nCells / 5. - 1)
                hnear = hMiddle;
              else
                hnear = _exactSol(i+1, 0);
              _exactSol(i,0) = exactHeight(p, q, a, b, hnear*(1+epsilon), hMax);
            }
          // Critical part (middle of the bump)
          double z(_topography(2. * _nCells / 5.));
          computeCoeffabcd(qIn, hMiddle, z, zMax, &a, &b, &c, &d);
          p = cardanP(a, b, c); q = cardanQ(a, b, c, d);
          _exactSol(2. * _nCells / 5., 0) = exactHeight(p, q, a, b, hMiddle, hMax);
          // Supercritical part (after the bump)
          for (int i(2. * _nCells / 5. + 1) ; i < _nCells ; ++i)
            {
              double z(_topography(i));
              computeCoeffabcd(qIn, hMiddle, z, zMax, &a, &b, &c, &d);
              p = cardanP(a, b, c); q = cardanQ(a, b, c, d);
              _exactSol(i,0) = exactHeight(p, q, a, b, _exactSol(i-1, 0)*(1-epsilon), hMax);
            }
        }
      else if (testCase == "TranscriticalFlowWithShock")
        {
          // Search for the limit between sub-super-sub
          double test(100.);
          double hplus(0.), hminus(0.);
          int abslim(2 * _nCells / 5);
          double epsi(10.0/_nCells);
          double zEnd(_topography(_nCells - 1));
          while(test > epsi && abslim < _nCells)
            {
              computeCoeffabcd(qIn, hOut, _topography(abslim), zEnd, &a, &b, &c, &d);
              p = cardanP(a, b, c);
              q = cardanQ(a, b, c, d);
              hplus = exactHeight(p, q, a, b, hOut, hMax);
              computeCoeffabcd(qIn, hMiddle, _topography(abslim), zMax, &a, &b, &c, &d);
              p = cardanP(a, b, c);
              q = cardanQ(a, b, c, d);
              hminus = exactHeight(p, q, a, b, hMiddle, hMax);
              test = RHJump(hplus, hminus, qIn);
              ++abslim;
            }
          _exactSol(abslim, 0) = hminus;

          for (int i(abslim - 1) ; i >= 0 ; --i)
            {
              computeCoeffabcd(qIn, hMiddle, _topography(i), zMax, &a, &b, &c, &d);
              p = cardanP(a, b, c);
              q = cardanQ(a, b, c, d);
              _exactSol(i,0) = exactHeight(p, q, a, b, _exactSol(i+1,0) + epsilon, hMax);
            }
          for (int i(_nCells - 1) ; i > abslim ; --i)
            {
              computeCoeffabcd(qIn, hOut, _topography(i), zEnd, &a, &b, &c, &d);
              p = cardanP(a, b, c);
              q = cardanQ(a, b, c, d);
              double hnear;
              if (i == _nCells - 1)
                hnear = hOut;
              else
                hnear = _exactSol(i+1,0);
              _exactSol(i,0) = exactHeight(p, q, a, b, hnear, hMax);
            }
        }
    }
}


// Function to solve by dichotomy for the dam break solutions
double Physics::damBreakFunc(double x, double vG, double vD) const
{
  return pow(x, 6.) - 9. * pow(vD, 2.) * pow(x, 4.) + 16. * vG * pow(vD, 2.) * pow(x, 3.) - pow(vD, 2.) * (pow(vD, 2.) + 8. * pow(vG, 2.)) * pow(x, 2.) + pow(vD, 6.);
}

// Methode de cardan
void Physics::computeCoeffabcd(double qIn, double hOut, double z, double zEnd, double* a, double* b, double* c, double *d)
{
  *a = 1.;
  *b = - (qIn * qIn / (2. * _g * hOut * hOut) + hOut - (z - zEnd));
  *c = 0.;
  *d = qIn * qIn / (2. * _g);
}

double Physics::cardanP(double a, double b, double c) const
{
  return (-b * b / (3.0 * a * a) + c / a);
}

double Physics::cardanQ(double a, double b, double c, double d) const
{
  return (b / (27.0 * a) * (2.0 * b * b / (a * a) - 9.0 * c / a) + d / a);
}

double Physics::cardanDet(double p, double q) const
{
  return (pow(q, 2.0) + 4.0 / 27.0 * pow(p, 3.0));
}

double Physics::RHJump(double hplus, double hminus, double q) const
{
  return abs((q * q * (1.0 / hplus - 1.0 / hminus) + 0.5 * _g * ((hplus - hminus) * (hplus + hminus))));
}

double Physics::exactHeight(double p, double q, double a, double b, double hnear, double hMax) const
{
  double det(cardanDet(p, q));
  double h(0.);
  double h1, h2;
  double l1, l2, l3;
  
  if (det > 0)
    {
      h1 = 0.5 * (- q + sqrt(det));
      h2 = 0.5 * (- q - sqrt(det));
      h = h1 / abs(h1) * pow(abs(h1), 1.0/3.0) + h2 / abs(h2) * pow(abs(h2), 1.0/3.0) - b / (3.0 * a);
    }
  else if (det == 0)
    {
      l1 = 3. * q / p;
      l2 = - 3. * q / (2. * p);
      if (l1 > 0)
        h = l1 - b / (3. * a);
      else
        h = l2 - b / (3. * a);
    }
  else
    {
      std::complex<double> u3(-q / 2.0, sqrt(abs(det)) / 2.);
      std::complex<double> u = pow(u3, 1.0 / 3.0);
      std::complex<double> j(- 1.0 / 2.0, sqrt(3.0) / 2.0);
      l1 = 2.0 * real(u) - b / (3.0 * a);
      l2 = 2.0 * real(j * u) - b / (3.0 * a);
      l3 = 2.0 * real(j * j * u) - b / (3.0 * a);
      if(l1 <= 0 && l2 <= 0 && l3 <= 0)
        std::cerr << "Error: no positive height" << std::endl;
      else if(l1 >= hMax && l2 >= hMax && l3 >= hMax)
        std::cerr << "Error: Probably irregular solution" << std::endl;
      else
        {
          double m(std::min(std::min(abs(l1 - hnear), abs(l2 - hnear)), abs(l3 - hnear)));
          if (m == abs(l1 - hnear))
            h = l1;
          else if (m == abs(l2 - hnear))
            h = l2;
          else
            h = l3;
        }
    }
  return h;
}


// Save the exact solution in a file
void Physics::saveExactSolution(std::string& fileName) const
{
#if VERBOSITY>0
  std::cout << "Saving exact solution" << std::endl;
#endif
  std::ofstream outputFile(fileName, std::ios::out);
  const Eigen::VectorXd& cellCenters(_mesh->getCellCenters());
  outputFile << "# x  H=h+z   h       u       q       Fr=|u|/sqrt(gh)" << std::endl;
  for (int i(0) ; i < _exactSol.rows() ; ++i)
    {
      outputFile << cellCenters(i) << " " <<
        _exactSol(i,0) + _topography(i) << " " <<
        _exactSol(i,0) << " " <<
        _exactSol(i,1)/_exactSol(i,0) << " " <<
        _exactSol(i,1) << " " <<
        abs(_exactSol(i,1)/_exactSol(i,0))/sqrt(_g * _exactSol(i,0)) << std::endl;
    }
}


//-------------------------------------------//
//---------------Physical flux---------------//
//-------------------------------------------//
Eigen::Vector2d Physics::physicalFlux(const Eigen::Vector2d& Sol) const
{
  Eigen::Vector2d flux;
  double h(Sol(0)), qx(Sol(1));
  if (h <= 0.)
    qx = 0.;
  flux(0) = qx;
  flux(1) = qx*qx/h + 0.5*_g*h*h;
  return flux;
}


//----------------------------------------//
//---------------Wave speed---------------//
//----------------------------------------//
void Physics::computeWaveSpeed(const Eigen::Vector2d& SolG, const Eigen::Vector2d& SolD, double* lambda1, double* lambda2) const
{
  double hG(SolG(0)), hD(SolD(0));
  double uG(SolG(1)/hG), uD(SolD(1)/hD);
  if (hG < 1e-6)
    uG = 0.;
  if (hD < 1e-6)
    uD = 0.;
  *lambda1 = std::min(uG - sqrt(_g * hG), uD - sqrt(_g * hD));
  *lambda2 = std::max(uG + sqrt(_g * hG), uD + sqrt(_g * hD));
}


//------------------------------------------------------//
//---------------Left Boundary Conditions---------------//
//------------------------------------------------------//
Eigen::Vector2d Physics::leftBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
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
      // std::cout << a << " " << b << " " << c << std::endl;
      double xe(FindRacine(a, b, c));
      // std::cout << xe << std::endl;
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
          SolG(0) = 3. + 0.1*sin(5 * M_PI * t) - _topography(0);
          SolG(1) = SolG(0)*(beta_moins_0_tnplus1 + 2*sqrt(_g*SolG(0)));
        }
      else if (_DF->getLeftBC() == "DataFile")
        {
          int i_max;
          i_max = _expBoundaryData.rows();
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
Eigen::Vector2d Physics::rightBoundaryFunction(double t, const Eigen::Matrix<double, Eigen::Dynamic, 2>& Sol)
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
          SolD(0) = _DF->getRightBCImposedHeight() - _topography(_nCells - 1);
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
              SolD(0) = _DF->getRightBCImposedHeight() - _topography(_nCells - 1);
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
              SolD(0) = _DF->getRightBCImposedHeight() - _topography(_nCells - 1);
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
double Physics::FindRacine(double a, double b, double c)
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
double Physics::FindSourceX(double x)
{
  int i(0);
  double dx(_DF->getDx());
  double x1(_xmin + (i + 0.5) * dx), x2(x1 + dx);
  double source1, source2, source;
  while (x1 < x)
    {
      ++i;
      x1 += dx;
    }
  x2 = x1 + dx;
  source1 = _source(i,1);
  source2 = _source(i+1,1);
  source = source1 + (x - x1)*(source2 - source1)/(x2 - x1);
  return source;
}
