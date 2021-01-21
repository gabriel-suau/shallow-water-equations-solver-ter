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

// Ajoute une arête
void Mesh::addEdge(const Edge& edge, int nt, std::vector<int>& headMinv, std::vector<int>& nextEdge, int& nbEdge)
{
  int vertex1 = edge.getVerticesIndex()(0);
  int vertex2 = edge.getVerticesIndex()(1);
  int ref = edge.getIndex();
  std::string BCType = edge.getBoundaryCondition();

  bool exist = false;
  // we look at the list of edges leaving from n1
  // if we find the same edge than n1->n2 we add the edge
  for (int e(headMinv[vertex1]) ; e != -1 ; e = nextEdge[e])
    {
      if (_edges[e].getVerticesIndex()(1) == vertex2)
        {
          if (nt >= 0)
            {
              _edges[e].addTriangle(nt);
            }
          exist = true;
        }
    }

  // if the edge has not been found, we create it
  if (!exist)
    {
      // we initialize the edge
      _edges[nbEdge] = Edge(vertex1, vertex2, ref, BCType);
      if (nt >= 0)
        {
          _edges[nbEdge].addTriangle(nt);
        }
      // we update the arrays next_edge and head_minv
      nextEdge[nbEdge] = headMinv[vertex1];
      headMinv[vertex1] = nbEdge;
      nbEdge++;
    }
}

// Calcule les centres et les aires des triangles
void Mesh::buildTrianglesCenterAndArea()
{
  _trianglesCenter.resize(_numberOfTriangles,2);
  _trianglesArea.resize(_numberOfTriangles);
  // Boucle sur les triangles
  for (int i(0) ; i < _numberOfTriangles ; ++i)
    {
      // Récupération des coordonnées des sommets du triangle
      int vertex1(_triangles[i].getVerticesReference()[0]);
      int vertex2(_triangles[i].getVerticesReference()[1]);
      int vertex3(_triangles[i].getVerticesReference()[2]);
      
      double x1(_vertices[vertex1].getCoordinates()[0]);
      double y1(_vertices[vertex1].getCoordinates()[1]);
      double x2(_vertices[vertex2].getCoordinates()[0]);
      double y2(_vertices[vertex2].getCoordinates()[1]);
      double x3(_vertices[vertex3].getCoordinates()[0]);
      double y3(_vertices[vertex3].getCoordinates()[1]);

      // Calcul du centre
      _trianglesCenter(i,0) = (x1 + x2 + x3)/3.;
      _trianglesCenter(i,1) = (y1 + y2 + y3)/3.;

      // Calcul de l'aire
      double l12(sqrt(pow(x1-x2,2) + pow(y1-y2,2)));
      double l13(sqrt(pow(x1-x3,2) + pow(y1-y3,2)));
      double l23(sqrt(pow(x2-x3,2) + pow(y2-y3,2)));
      double p(0.5 * (l12 + l13 + l23));
      _trianglesArea(i) = sqrt(p*(p-l12)*(p-l13)*(p-l23));
    }
}

// Calcule le centre, la longueur et la normale de chaque arête
void Mesh::buildEdgesNormalAndLengthAndCenter()
{
  _edgesCenter.resize(_numberOfEdges,2);
  _edgesLength.resize(_numberOfEdges);
  _edgesNormal.resize(_numberOfEdges,2);

  for (int i(0) ; i < _numberOfEdges ; ++i)
    {
      // Calcul de la longueur
      int vertex1(_edges[i].getVerticesIndex()(0));
      int vertex2(_edges[i].getVerticesIndex()(1));
      double x1(_vertices[vertex1].getCoordinates()(0));
      double y1(_vertices[vertex1].getCoordinates()(1));
      double x2(_vertices[vertex2].getCoordinates()(0));
      double y2(_vertices[vertex2].getCoordinates()(1));
      _edgesLength(i) = sqrt(pow(x1-x2,2) + pow(y1-y2,2));
      
      // Calcul du centre
      _edgesCenter(i,0) = 0.5 * (x1 + x2);
      _edgesCenter(i,1) = 0.5 * (y1 + y2);
      
      // Calcul du vecteur (centre du triangle t1 to centre de l'arête)
      int t1(_edges[i].getT1());
      Eigen::Vector2d diff(_edgesCenter.row(i) - _trianglesCenter.row(t1));
      // Calcul de la normale dans un sens arbitraire
      _edgesNormal(i,0) = y1 - y2;
      _edgesNormal(i,1) = x2 - x1;
      // Produit scalaire entre la normale et le vecteur diff
      double scalar(_edgesNormal.row(i).dot(diff));
      // Forcer la normale à sortir de T1
      if (scalar < 0.)
        {
          _edgesNormal.row(i) = -_edgesNormal.row(i);
        }
      // Normalisation
      _edgesNormal.row(i) /= _edgesLength(i);
    }
}

void Mesh::Initialize()
{
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "Generating a 2D unstructured mesh from file" << _meshFile << "..." << std::endl;

  std::ifstream meshStream(_meshFile, std::ios::in);

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

  std::string line;
  int dimension(3), nBoundaryEdges(0);
  std::vector<Edge> boundaryEdges;

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
              --vertex1; --vertex2;
              std::string BCType("none");
              for (int i(0) ; i < _boundaryConditionReference.size() ; ++i)
                {
                  if (index == _boundaryConditionReference[i])
                    {
                      BCType = _boundaryConditionType[i];
                    }
                }
              if (BCType == "none")
                {
                  std::cout << termcolor::red << "ERROR::MESH : Problem with boundary conditions in your mesh (reference or types are wrong)" << std::endl;
                  std::cout << termcolor::reset;
                  exit(-1);
                }
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
              --vertex1; --vertex2; --vertex3;
              _triangles[i] = Triangle(vertex1, vertex2, vertex3, index);
            }
        }
    }

  // Création du vecteur _edges.
  _numberOfEdges = (3*_numberOfTriangles + nBoundaryEdges)/2;
  // _numberOfEdges = _numberOfTriangles + _numberOfVertices - 1;
  _edges.resize(_numberOfEdges);

  std::vector<int> headMinv(_numberOfVertices, -1);
  std::vector<int> nextEdge(_numberOfEdges, -1);

  // Ajout des arêtes intérieures
  int nbEdges(0);
  for (int i(0) ; i < nBoundaryEdges ; ++i)
    {
      addEdge(boundaryEdges[i], -1, headMinv, nextEdge, nbEdges);
    }

  // Ajout des arêtes extérieures
  for (int i(0); i < _numberOfTriangles; i++)
    {
      const Eigen::Vector3i& nv = _triangles[i].getVerticesReference();
      for (int j = 0; j < 3; j++)
        {
          Edge edge(nv(j), nv((j+1)%3), 0, "none");
          addEdge(edge, i, headMinv, nextEdge, nbEdges);
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
  std::cout << "Number of vertices  = " << _numberOfVertices << std::endl;
  std::cout << "Number of edges     = " << _numberOfEdges << std::endl;
  std::cout << "Number of triangles = " << _numberOfTriangles << std::endl;
  std::cout << "====================================================================================================" << std::endl << std::endl;
}
