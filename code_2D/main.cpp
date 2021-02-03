#include "termcolor.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Physics.h"
#include "FiniteVolume.h"
#include "TimeScheme.h"

#include <iostream>

int main(int argc, char** argv)
{
  //-------------------------------------------------------//
  //---------------------Vérifications---------------------//
  //-------------------------------------------------------//
  if (argc < 2)
    {
      std::cout << termcolor::red << "ERROR::MAIN : Please, enter the name of your data file." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
  else if (argc > 2)
    {
      std::cout << termcolor::yellow << "WARNING::MAIN : Too many arguments : only 1 was expected but " << argc - 1 << " were given..."<< std::endl;
      std::cout << termcolor::reset;
    }

  //-------------------------------------------------------//
  //---------------------Logs de début---------------------//
  //-------------------------------------------------------//
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Solving 2D St-Venant equations for you !" << std::endl;
  std::cout << "====================================================================================================" << std::endl << std::endl;

  //-------------------------------------------------------//
  //---------------------Fichier de paramètres-------------//
  //-------------------------------------------------------//
  DataFile* DF = new DataFile(argv[1]);
  DF->readDataFile();
  DF->printData();

  //--------------------------------------------------//
  //---------------------Maillage---------------------//
  //--------------------------------------------------//
  Mesh* mesh = new Mesh(DF);
  mesh->Initialize();
  mesh->printParameters();





  //vérification aire
  double AireTotale(0);
  for (int i = 0; i < mesh->getNumberOfCells(); i++)
  {
    AireTotale += mesh->getCellsArea()[i];
    std::cout << "Aire n. " << i << "\n" << mesh->getCellsArea()[i] << "\n";
  }
  std::cout << "Aire Totale" << "\n" << AireTotale << "\n";
  
               
  //----------------------------------------------------------------//
  //---------------------CI, CL, Termes sources---------------------//
  //----------------------------------------------------------------//
  Physics* physics = new Physics(DF, mesh);
  physics->Initialize();
  
  //--------------------------------------------------------//
  //---------------------Flux numérique---------------------//
  //--------------------------------------------------------//
  FiniteVolume* finVol;
  if (DF->getNumericalFlux() == "Rusanov")
    {
      finVol = new Rusanov(DF, mesh, physics);
    }
  else if (DF->getNumericalFlux() == "HLL")
    {
      finVol = new HLL(DF, mesh, physics);
    }
  else
    {
      std::cout << termcolor::red << "ERROR::FINITEVOLUME : Case not implemented." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }

  //---------------------------------------------------------//
  //---------------------Schéma en temps---------------------//
  //---------------------------------------------------------//
  TimeScheme* TS;
  if (DF->getTimeScheme() == "ExplicitEuler")
    {
      TS = new ExplicitEuler(DF, mesh, physics, finVol);
    }
  else
    {
      std::cout << termcolor::red << "ERROR::TIMESCHEME : Case not implemented." << std::endl;
      std::cout << termcolor::reset;
      exit(-1);
    }
  

  //----------------------------------------------------//
  //---------------------Résolution---------------------//
  //----------------------------------------------------//
  TS->solve();
  
  //-----------------------------------------------------------//
  //---------------------Libère la mémoire---------------------//
  //-----------------------------------------------------------//
  delete DF;
  delete mesh;
  delete physics;
  delete finVol;
  delete TS;

  //-----------------------------------------------------//
  //---------------------Logs de fin---------------------//
  //-----------------------------------------------------//
  std::cout << "====================================================================================================" << std::endl;
  std::cout << termcolor::green << "SUCCESS : Successfully solved the 2D St-Venant equations for you !" << std::endl;
  std::cout << termcolor::reset << "Let me terminate myself now..." << std::endl;
  std::cout << "====================================================================================================" << std::endl << std::endl;

  // Fin
  return 0;
}
