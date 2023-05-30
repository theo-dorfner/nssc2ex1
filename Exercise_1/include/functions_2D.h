#include<iostream>
#include<cmath>
#include<vector>

using namespace std;


vector<int> BiggestDivisors(int procs)
{
    vector<int>DIV{procs, 1};

    double i = sqrt(procs);
    int l = std::floor(i);

    for(int j=0; j<procs-l-1; j++)
    {
        if(procs%(l+j) == 0)
        {
            DIV[0]=l+j;
            DIV[1]=procs/DIV[0];
            break;
        }
    }

    if(DIV[0] < DIV[1]) 
    //i want the first entry in DIV to be bigger than the second!
    {
        int a = DIV[0];
        DIV[0] = DIV[1];
        DIV[1] = a;
    }

    return DIV;
}


vector<int> UnknownsPerProc(vector<int>PPP, int precs, int proc)
//this function declares how big the A, u and f are in the "Blocks"
{
    vector<int>A;
    int min_prec = floor(precs / proc);
    int rest = precs % proc;

    for(int i=0; i<proc; i++)
    {
      A.push_back(min_prec);
    }

    for(int j=0; j<rest; j++)
    {
      A[j]++;
    }

    for(int l=0; l<proc; l++){
      PPP.push_back(A[l]*precs);
    }

    return PPP;
}


vector<int> UPPtoYBegin(vector<int>UPP, int prec, int proc) //transforms f.e. (12,12,6,6) to (1, 3, 5, 6, ...); helpfunction for b0 initialization
{
    vector<int>Y_begin;
    int N = UPP.size();

    for(int i = 0; i < N; i++)
    {
        UPP[i] = UPP[i] / prec;
    }

    Y_begin.push_back(1);

    for(int j = 0; j<N-1; j++)
    {
        Y_begin.push_back(Y_begin[j] + UPP[j]);
    }
    Y_begin.push_back(prec+1);

    return Y_begin;
}


vector<int> N_DIST(int DIVxy, int prec) //separation of nodes in x_direction
{
    vector<int> N;
    int min_prec = floor(prec / DIVxy);
    int rest = prec % DIVxy;

    for(int i=0; i<DIVxy; i++)
    {
        N.push_back(min_prec);
    }

    for(int j=0; j<rest; j++)
    {
        N[j]++;
    }

    return N;
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


vector<double>Initialize_b0_1D(vector<double>b, vector<int>Y_begin, int prec, double h, int my_rank, int proc)
//this function initializes b0
{
    double X;
    double Y;

    for(int j=Y_begin[my_rank]; j<Y_begin[my_rank+1]; j++) //die Y-Koordinate lauft werte bis zum naechsten rank durch
    {
        for(int i=1; i<prec+1; i++)
        {
            X = i*h;
            Y = j*h;
            if(abs(Y - (1-h)) < 1e-5)
            {
                b.push_back(f(X, Y)*h*h + BC(X));   //u_p including the BC
            }
            else
            {
                b.push_back(f(X, Y)*h*h); //only u_p without BC
            }
        }
    }

    //if(my_rank == 1)vector_printer(b);
    return b;
}

vector<double>Initialize_b0_2D(vector<double>b, int x_coord, int y_coord, vector<int>x_begin, vector<int> y_begin, double h)
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
            b.push_back(f(X[k], Y[k])*h*h + BC(X[k]));   //u_p including the BC
        }
        else
        {
            b.push_back(f(X[k], Y[k])*h*h); //only u_p without BC
        }
    }
        
    return b;
}

double u_p(double X, double Y)
{
    return sin(2 * M_PI * X) * sinh(2 * M_PI * Y);
} 

vector<double>Initialize_up(vector<double>b, vector<int>Y_begin, int prec, double h, int my_rank, int proc)
{
    double X;
    double Y;

    for(int j=Y_begin[my_rank]; j<Y_begin[my_rank+1]; j++) //die Y-Koordinate lauft werte bis zum naechsten rank durch
    {
        for(int i=1; i<prec+1; i++)
        {
            X = i*h;
            Y = j*h;
            b.push_back(u_p(X, Y)); 
        }
    }

    return b;
}


