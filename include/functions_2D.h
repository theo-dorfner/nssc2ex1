#include<iostream>
#include<cmath>
#include<vector>

using namespace std;

vector<int> BiggestDivisors(int procs)
{
      vector<int>DIV{1,procs};

      double i = sqrt(procs);

      int l = std::floor(i);
      for(int j=procs; j>0; j--)
      {
          double m = j*l;
          if(abs(m - procs) < 1e-5)
          {
              DIV[0]=j;
              DIV[1]=l;
          }
      }

    return DIV;
}

vector<int> NX_fun(int DIVx, int prec) //separation of nodes in x_direction
{
    vector<int> Nx;
    int min_prec = floor(prec / DIVx);
    int rest = prec % DIVx;

    for(int i=0; i<DIVx; i++)
    {
        Nx.push_back(min_prec);
    }

    for(int j=0; j<rest; j++)
    {
        Nx[j]++;
    }

    return Nx;
}

vector<int> NY_fun(int DIVy, int prec) //separation of nodes in x_direction
{
    vector<int> Ny;
    int min_prec = floor(prec / DIVy);
    int rest = prec % DIVy;

    for(int i=0; i<DIVy; i++)
    {
        Ny.push_back(min_prec);
    }

    for(int j=0; j<rest; j++)
    {
        Ny[j]++;
    }

    return Ny;
}

vector<int>XBegin(vector<int>Nx, int prec)
{
    int N = Nx.size();
    vector<int>x_begin;
    x_begin.push_back(1);

    for(int j = 0; j<N-1; j++)
    {
        x_begin.push_back(x_begin[j] + Nx[j]);
    }
    x_begin.push_back(prec+1);

    return x_begin;
}

vector<int>YBegin(vector<int>Ny, int prec)
{
    int N = Ny.size();
    vector<int>y_begin;
    y_begin.push_back(1);

    for(int j = 0; j<N-1; j++)
    {
        y_begin.push_back(y_begin[j] + Ny[j]);
    }
    y_begin.push_back(prec+1);

    return y_begin;
}

vector<int>UPP_fun(vector<int>Nx, vector<int>Ny) //gives out number of unknowns for each rank!
{
    vector<int>UPP;
    
    int ymax = Ny.size();
    int xmax = Nx.size();
    for(int j=0; j<ymax; j++)
    {
        for(int i=0; i<xmax ; i++)
        {
            UPP.push_back(Ny[j]*Nx[i]);
        }
    }

    return UPP;
}

void vectorprinter(vector<int>B)
{
    int max = B.size();
    for(int i=0; i<max; i++)
    {
        std::cout << B[i] << std::endl;
    }
    std::cout << std::endl;
}

void vectorprinter_d(vector<double>B)
{
    int max = B.size();
    for(int i=0; i<max; i++)
    {
        std::cout << B[i] << std::endl;
    }
    std::cout << std::endl;
}

double H(int res)
{
      double h = 1.0 / (res - 1);
      return h;
}


vector<double>Initialize_A0(vector<double>A, int N, int Nx, double h)
{
    double alpha = 4 + 4 * M_PI * M_PI * h * h;
    vector<double> stencil {-1, alpha,-1};
    
    for(int j=0; j<N; j++) //loop creates banded matrix 
    {
        if(j == 0)
        {
            for(int i=0; i<2; i++)
            {
                A.at(j*N + i) = stencil.at(i+1);
            }
        }   
        else if(j == N-1)
        {
            for(int i=N-2; i<N; i++)
            {
                A.at(j*N + i) = stencil.at(j-i+1);
            }
        }
        else
        {
            for(int i=j-1; i<j+2; i++)
            {
                A.at(j*N + i) = stencil.at(j-i+1);
            }
        }
    }

    for(int i=1; i<N; i++) 
    {
        if(i%Nx == 0)
        {
            A[i*N + i - 1] = 0;
            A[(i-1)*N + i] = 0;
        }
    }
        
    for(int j = 0; j < N; j++)
    {
        if(j < Nx)
        {
            A[j*N + j + Nx] = -1;
        }
        else if(j >= Nx && j < N-Nx)
        {
            A[j*N + j + Nx] = -1;
            A[j*N + j - Nx] = -1;
        }  
        else if(j >= N-Nx)
        {
            A[j*N + j - Nx] = -1;
        }
    }
    
    return A;
}

void matrix_printer(vector<double>A, int N)
{
  for(int j=0; j<N; j++)
  {
      for(int i=0; i<N; i++)
      {
          cout << A[j*N + i] << " ";
      }
  cout<<endl;
  }
  cout<<endl;
}

double f(double X, double Y)
{
    return M_PI * M_PI * 4 * sin(2 * M_PI * X) * sinh(2 * M_PI * Y);
}

double BC(double X)
{
    return sin(2 * M_PI * X) * sinh(2*M_PI);
}


vector<double>Initialize_b0(vector<double>b, int x_coord, int y_coord, vector<int>x_begin, vector<int> y_begin, double h)
{
    
    vector<double>X;
    vector<double>Y;

    for(int j=y_begin[y_coord]; j<y_begin[y_coord+1]; j++) //die Y-Koordinate lauft werte bis zum naechsten rank durch
    {
        for(int i=x_begin[x_coord]; i<x_begin[x_coord+1]; i++)
        {
            X.push_back(i*h);
            Y.push_back(j*h);
        }
    }

    int max = X.size();
    for(int k=0; k<max; k++)
    {
        if(abs(Y[k] - (1-h)) < 1e-5)
        {
            b.push_back(f(X[k], Y[k]) + BC(X[k]));   //u_p including the BC
        }
        else
        {
            b.push_back(f(X[k], Y[k])); //only u_p without BC
        }
    }
        
    return b;
}


