#include "Mesh.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"
#include "DataFile.h"
#include "termcolor.h"
#include <fstream>

//--------------------------------------------------//
//---------------------Vertices---------------------//
//--------------------------------------------------//
Vertex::Vertex():
  _index(-1), _coordinates(0.,0.)
{
}

Vertex::Vertex(double x, double y, int index):
  _index(index), _coordinates(x,y)
{
}

void Vertex::print() const
{
  std::cout << "(x,y,index) = (" << _coordinates(0)  << "," << _coordinates(1) << "," << _index << ")" << std::endl;
}

//-----------------------------------------------//
//---------------------Edges---------------------//
//-----------------------------------------------//
Edge::Edge():
  _index(-1), _verticesIndex(-1,-1)
{
}

Edge::Edge(int vertex1, int vertex2, int index, const std::string& boundaryCondition):
  _index(index), _t1(-1), _t2(-1), _boundaryCondition(boundaryCondition)
{
  if (vertex1 > vertex2)
    {
      _verticesIndex(0) = vertex2;
      _verticesIndex(1) = vertex1;
    }
  else
    {
      _verticesIndex(0) = vertex1;
      _verticesIndex(1) = vertex2;
    }
};

//---------------------------------------------------//
//---------------------Triangles---------------------//
//---------------------------------------------------//
Triangle::Triangle()
{
}

Triangle::Triangle(int vertex1,int vertex2,int vertex3,int index):
  _index(index), _verticesIndex(vertex1, vertex2, vertex3)
{
}

//----------------------------------------------//
//---------------------Mesh---------------------//
//----------------------------------------------//
Mesh::Mesh()
{
}

Mesh::Mesh(DataFile* DF):
  _DF(DF), _meshFile(_DF->getMeshFile())
{
}

void Mesh::Initialize(DataFile* DF)
{
  _DF = DF;
  _meshFile = _DF->getMeshFile();
  this->Initialize();
}

void Mesh::Initialize()
{
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Generating a 2D unstructured mesh from file" << _meshFile << "..." << std::endl;

  std::ifstream meshStream(_meshFile);
  std::string line;
  int dimension(3), numberOfVertices(0);

  // Vérifie que le fichier est bien ouvert.
  if (!meshStream.is_open())
    {
      std::cout << termcolor::red << "ERROR::MESH : Unable to open the mesh file : " << _meshFile << std::endl;
      std::cout << termcolor::reset << "====================================================================================================" << std::endl;
      exit(-1);
    }
  else
    {
      std::cout << "Building the mesh from file : " << _meshFile << std::endl;
    }

  // Parcours les lignes du fichier de maillage
  while (getline(meshStream, line))
    {
      // Dimension du maillage
      if (line.find("Dimension") != std::string::npos)
        {
          meshStream >> dimension;
        }
      // Création des sommets du maillage
      else if (line.find("Vertices") != std::string::npos)
        {
          meshStream >> numberOfVertices;
          _vertices.resize(numberOfVertices);
          for (int i(0) ; i < numberOfVertices ; ++i)
            {
              double x, y, z;
              int index;
              meshStream >> x >> y >> z >> index;
              _vertices[i] = Vertex(x, y, index);
            }
        }
      else if (line.find("Edges") != std::string::npos)
        {
          // TODO
        }
      else if (line.find("Triangles") != std::string::npos)
        {
          // TODO
        }
    }
  std::cout << termcolor::green << "SUCCESS::MESH : Mesh generated succesfully !" << std::endl;
  std::cout << termcolor::reset << "====================================================================================================" << std::endl << std::endl;
}

void Mesh::printParameters() const
{
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Printing parameters of the mesh..." << std::endl;
  std::cout << "Mesh file           = " << _meshFile << std::endl;
  std::cout << "Number of vertices  = " << _vertices.size() << std::endl;
  std::cout << "Number of edges     = " << _edges.size() << std::endl;
  std::cout << "Number of triangles = " << _triangles.size() << std::endl;
  std::cout << "====================================================================================================" << std::endl;
}
