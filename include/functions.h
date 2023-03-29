#include<cmath>
#include<vector>
#include<random>
#include<algorithm>

using namespace std;

void vector_printer(vector<double>b);

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

double f(double X, double Y)
{
    return M_PI * M_PI * 4 * sin(2 * M_PI * X) * sinh(2 * M_PI * Y);
}

double BC(double X)
{
    return sin(2 * M_PI * X) * sinh(2*M_PI);
}

vector<double>Initialize_b0(vector<double>b, vector<int>Y_begin, int prec, double h, int my_rank, int proc)
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


vector<double>Initialize_u0(vector<double>u, int N)
{
    for(int i=0; i<N; i++)
    {
        u.push_back(0);
    }

    return u;
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


double H(int res)
//this function calculates h
{
      double h = 1.0 / (res - 1);
      return h;
}

double A_at(vector<double>A, vector<int>X, vector<int>Y, int x_value, int y_value)
{
    int max = X.size();
    //int meister = 0;
    
    for(int i=0; i<max; i++)
    {
        if(X[i]==x_value)
        {
            if(Y[i]==y_value)
            {
                return A[i];
            }
        }
    }
    return 0;
}

void vector_printer(vector<double>b)
//function prints vector<double>
{
    for(int i=0; i<static_cast<int>(b.size()); i++)
    {
        cout << b[i] << endl;
    }
    cout << endl;
}

void vector_printer_int(vector<int>b)
//function prints vector<int>
{
    for(int i=0; i<static_cast<int>(b.size()); i++)
    {
        cout << b[i] << endl;
    }
    cout << endl;
}

void matrix_printer(vector<double>A, int N)
//function prints vector<double> as a matrix
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
