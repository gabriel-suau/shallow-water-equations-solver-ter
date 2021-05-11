#ifndef DATA_FILE_H
#define DATA_FILE_H

#include <iostream>
#include <string>
#include <vector>



class DataFile
{
private:
  // Nom du fichier de paramètres
  std::string _fileName;

  // Sauvegarde des resultats
  std::string _resultsDir;
  bool _isSaveFinalTimeOnly;
  int _saveFrequency;
  // Water height probes (to follow the water height at a precise position with time).
  int _nProbes;
  std::vector<int> _probesReferences;
  std::vector<double> _probesPositions;
  
  // Test cases
  bool _isTestCase;
  std::string _testCase;
  
  // Conditions initiales
  std::string _initialCondition;
  std::string _initFile;
  double _initialHeight, _initialDischarge;
  
  // Mesh parameters
  double _xmin, _xmax;
  double _dx;
  int _Nx;
  
  // Numerical Flux
  std::string _numericalFlux;
  int _schemeOrder;
  
  // Time parameters
  std::string _timeScheme;
  double _initialTime;
  double _finalTime;
  double _timeStep;
  double _CFL;

  // Gravity Acceleration
  double _g;

  // Boundary conditions
  std::string _leftBC, _rightBC;
  std::string _leftBCDataFile, _rightBCDataFile;
  double _leftBCImposedHeight, _leftBCImposedDischarge, _rightBCImposedHeight, _rightBCImposedDischarge;
  
  // Topography
  bool _isTopography;
  std::string _topographyType;
  std::string _topographyFile;
  
public:
  // Constructeurs
  DataFile();
  DataFile(const std::string& fileName);

  // Destructeur
  ~DataFile() = default;

  // Initialise l'objet
  void Initialize(const std::string& fileName);
  
  // Lit le fichier
  void readDataFile();

  // Nettoyer une ligne du fichier
  std::string cleanLine(std::string &line);

  // Getters
  // Data file name
  const std::string& getFileName() const {return _fileName;};
  // Solution saving 
  const std::string& getResultsDirectory() const {return _resultsDir;};
  bool isSaveFinalTimeOnly() const {return _isSaveFinalTimeOnly;};
  int getSaveFrequency() const {return _saveFrequency;};
  int getNumberOfProbes() const {return _nProbes;};
  const std::vector<int>& getProbesReferences() const {return _probesReferences;};
  const std::vector<double>& getProbesPositions() const {return _probesPositions;};
  // Test cases
  bool isTestCase() const {return _isTestCase;};
  const std::string& getTestCase() const {return _testCase;};
  // Initial Conditions
  const std::string& getInitialCondition() const {return _initialCondition;};
  const std::string& getInitFile() const {return _initFile;};
  double getInitialHeight() const {return _initialHeight;};
  double getInitialDischarge() const {return _initialDischarge;};
  // Mesh related
  double getXmin() const {return _xmin;};
  double getXmax() const {return _xmax;};
  double getDx() const {return _dx;};
  int getNx() const {return _Nx;};
  // Numerical flux related
  const std::string& getNumericalFlux() const {return _numericalFlux;};
  int getSchemeOrder() const {return _schemeOrder;};
  // Time scheme related
  const std::string& getTimeScheme() const {return _timeScheme;};
  double getInitialTime() const {return _initialTime;};
  double getFinalTime() const {return _finalTime;};
  double getTimeStep() const {return _timeStep;};
  double getCFL() const {return _CFL;};
  // Gravity related
  double getGravityAcceleration() const {return _g;};
  // Boundary conditions related
  const std::string& getLeftBC() const {return _leftBC;};
  const std::string& getRightBC() const {return _rightBC;};
  const std::string& getLeftBCDataFile() const {return _leftBCDataFile;};
  const std::string& getRightBCDataFile() const {return _rightBCDataFile;};
  double getLeftBCImposedHeight() const {return _leftBCImposedHeight;};
  double getLeftBCImposedDischarge() const {return _leftBCImposedDischarge;};
  double getRightBCImposedHeight() const {return _rightBCImposedHeight;};
  double getRightBCImposedDischarge() const {return _rightBCImposedDischarge;};
  // Topography related
  bool isTopography() const {return _isTopography;};
  const std::string& getTopographyType() const {return _topographyType;};
  const std::string& getTopographyFile() const {return _topographyFile;};
  
  // Affichage des paramètres
  void printData() const;
};

#endif // DATA_FILE_H
