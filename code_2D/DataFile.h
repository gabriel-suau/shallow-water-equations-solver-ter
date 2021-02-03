#ifndef DATA_FILE_H
#define DATA_FILE_H

#include "Eigen/Eigen/Dense"

#include <iostream>
#include <string>
#include <vector>

class DataFile
{
private:
  // Name of the data file
  std::string _fileName;

  // General parameters
  std::string _scenario;
  std::string _resultsDir;

  // Mesh parameters
  std::string _meshFile;

  // Numerical Flux
  std::string _numericalFlux;

  // Time parameters
  std::string _timeScheme;
  double _initialTime;
  double _finalTime;
  double _timeStep;
  double _CFL;

  // Gravity Acceleration
  double _g;

  // Saving frequency (number of iterations) of the solution
  int _saveFrequency;

  // Topography
  bool _isTopography;
  std::string _topographyType;
  std::string _topographyFile;

  // Boundary conditions
  int _nBoundaries;
  Eigen::VectorXi _boundaryConditionReference;
  std::vector<std::string> _boundaryConditionType;
  
public:
  // Constructors
  DataFile();
  DataFile(const std::string& fileName);

  // Destructor
  ~DataFile() = default;

  // Initializer
  void Initialize(const std::string& fileName);

  // Read the data file
  void readDataFile();

  // Nettoyer une ligne du fichier
  std::string cleanLine(std::string &line);

  // Getters
  const std::string& getFileName() const {return _fileName;};
  const std::string& getScenario() const {return _scenario;};
  const std::string& getResultsDirectory() const {return _resultsDir;};
  const std::string& getMeshFile() const {return _meshFile;};
  const std::string& getNumericalFlux() const {return _numericalFlux;};
  const std::string& getTimeScheme() const {return _timeScheme;};
  double getInitialTime() const {return _initialTime;};
  double getFinalTime() const {return _finalTime;};
  double getTimeStep() const {return _timeStep;};
  double getCFL() const {return _CFL;};
  double getGravityAcceleration() const {return _g;};
  int getSaveFrequency() const {return _saveFrequency;};
  bool isTopography() const {return _isTopography;};
  const std::string& getTopographyType() const {return _topographyType;};
  const std::string& getTopographyFile() const {return _topographyFile;};
  const Eigen::VectorXi& getBoundaryConditionReference() const {return _boundaryConditionReference;};
  const std::vector<std::string>& getBoundaryConditionType() const {return _boundaryConditionType;};

  // Printer (for information purposes)
  void printData() const;
};

#endif // DATA_FILE_H
