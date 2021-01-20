#include "Mesh.h"
#include "DataFile.h"
#include "termcolor.h"

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

#include <fstream>
#include <vector>

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

Edge::Edge(int vertex1, int vertex2, int index, std::string boundaryCondition):
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
}

void Edge::print() const
{
  std::cout << "(vertex1,vertex2,index) = (" << _verticesIndex(0)  << "," << _verticesIndex(1) << "," << _index << ")" << std::endl;
}


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

void Triangle::print() const
{
  std::cout << "(vertex1,vertex2,vertex3,index) = (" << _verticesIndex(0)  << "," << _verticesIndex(1) << "," << _verticesIndex(2) << "," << _index << ")" << std::endl;
}

//----------------------------------------------//
//---------------------Mesh---------------------//
//----------------------------------------------//
Mesh::Mesh()
{
}

Mesh::Mesh(DataFile* DF):
  _DF(DF), _meshFile(_DF->getMeshFile()), _boundaryConditionReference(_DF->getBoundaryConditionReference()), _boundaryConditionType(_DF->getBoundaryConditionType())
{
}

void Mesh::Initialize(DataFile* DF)
{
  _DF = DF;
  _meshFile = _DF->getMeshFile();
  _boundaryConditionReference = _DF->getBoundaryConditionReference();
  _boundaryConditionType = _DF->getBoundaryConditionType();
  this->Initialize();
}

void Mesh::Initialize()
{
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Generating a 2D unstructured mesh from file" << _meshFile << "..." << std::endl;

  std::ifstream meshStream(_meshFile);
  std::string line;
  int dimension(3), nBoundaryEdges(0);
  std::vector<Edge> boundaryEdges;

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
          meshStream >> _numberOfVertices;
          _vertices.resize(_numberOfVertices);
          for (int i(0) ; i < _numberOfVertices ; ++i)
            {
              double x(0.), y(0.), z(0.);
              int index(0);
              meshStream >> x >> y >> z >> index;
              _vertices[i] = Vertex(x, y, index);
            }
        }
      // Création des arêtes extérieures du maillage
      else if (line.find("Edges") != std::string::npos)
        {
          meshStream >> nBoundaryEdges;
          boundaryEdges.resize(nBoundaryEdges);
          for (int i(0) ; i < nBoundaryEdges ; ++i)
            {
              int vertex1(0), vertex2(0);
              int index(0);
              meshStream >> vertex1 >> vertex2 >> index;
              std::string BCType(_boundaryConditionType[index - 1]);
              boundaryEdges[i] = Edge(vertex1, vertex2, index, BCType);
            }
        }
      // Création des triangles du maillage
      else if (line.find("Triangles") != std::string::npos)
        {
          meshStream >> _numberOfTriangles;
          _triangles.resize(_numberOfTriangles);
          for (int i(0) ; i < _numberOfTriangles ; ++i)
            {
              int vertex1, vertex2, vertex3;
              int index;
              meshStream >> vertex1 >> vertex2 >> vertex3 >> index;
              _triangles[i] = Triangle(vertex1, vertex2, vertex3, index);
            }
        }
    }

  // Création du vecteur _edges.
  _numberOfEdges = _numberOfTriangles + _numberOfVertices - 1;
  _edges.resize(_numberOfEdges);
  
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
