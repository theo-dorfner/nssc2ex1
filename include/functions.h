#include<cmath>
#include<vector>
#include<random>

using namespace std;

vector<int> UnknownsPerProc(vector<int>PPP, int precs, int proc) //this function declares how big the A, u and f are in the "Blocks"
{
    int A[proc];
    int min_prec = floor(precs / proc);
    int rest = precs % proc;

    for(int i=0; i<proc; i++)
    {
      A[i] = min_prec;
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

vector<double>Initialize_u0(vector<double>u, int N)
{
    for(int i=0; i<N; i++)
    {
        u.push_back(0);
    }

    return u;
}

vector<double>Initialize_b0(vector<double>b, int N)
{
    for(int i=0; i<N; i++)
    {
        b.push_back(0);
    }

    return b;
}

double H(int res)
{
      double h = 1.0 / (res - 1);
      return h;
}

vector<double>Initialize_A0(vector<double>A, int N, int width, double h)
{
    double alpha = 4 + 4 * M_PI * M_PI * h * h;
    vector<double> stencil {-1, alpha,-1};

    for(int j=0; j<N; j++)
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

    if(N/width > 1) //if its not an 1D Problem we need to include the "north" and "south" neigbours
    {
        for(int j = 0; j < N; j++)
        {
            if(j < 3)
            {
                A[j*N + j + 3] = -1;
            }
            else if(j >= 3 && j < N-3)
            {
                A[j*N + j + 3] = -1;
                A[j*N + j - 3] = -1;
            }
            else if(j >= N-3)
            {
                A[j*N + j - 3] = -1;
            }
        }
    }
    return A;
}


vector<double>Initialize_v_random(vector<double>b, int N)
{
    random_device dev; //random number generator from 1 to 100
    mt19937 rng(dev());
    uniform_int_distribution<mt19937::result_type> dist100(1,100);

    for(int i=0; i<N; i++)
    {
        b.push_back(dist100(rng));
    }

    return b;
}

vector<double>Initialize_A_random(vector<double>A, int N)
{
    random_device dev; //random number generator from 1 to 100
    mt19937 rng(dev());
    uniform_int_distribution<mt19937::result_type> dist100(1,100);

    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            A.push_back(dist100(rng));
        }
    }

    return A;
}

void vector_printer(vector<double>b)
{
    for(int i=0; i<b.size(); i++)
    {
        cout << b[i] << endl;
    }
    cout << endl;
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


/*
mpiCC meister1d.cpp -o meister
mpirun -n 6 meister1d
*/
