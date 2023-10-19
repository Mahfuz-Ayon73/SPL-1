#include <iostream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <random>
#include<cmath>

#define EPS 10e-6

#define Lamda 0.5

using namespace std;



vector<vector<double>> matrixMultiplication(int row, int rank, int col, vector<vector<double>> &matrix1, vector<vector<double>> &matrix2)
{

  vector<vector<double>> multiple_matrix;

  for (int i = 0; i < row; i++)
  {
    vector<double> tempvec;
    for (int j = 0; j < col; j++)
    {
      tempvec.push_back(0);
    }
    multiple_matrix.push_back(tempvec);
  }
  for (int k = 0; k < row; k++)
  {
    for (int i = 0; i < col; i++)
    {
      double sum = 0;

      for (int j = 0; j < rank; j++)
      {
        sum += matrix1[k][j] * matrix2[j][i];
      }
      multiple_matrix[k][i] = sum;
      sum = 0.0;
    }
  }

  return multiple_matrix;
}

vector<vector<double>> transpose(int row, int col, vector<vector<double>> &matrix)
{
  vector<vector<double>> transpose_matrix;

  for (int i = 0; i < col; i++)
  {
    vector<double> tempvec;
    for (int j = 0; j < row; j++)
    {
      tempvec.push_back(0);
    }
    transpose_matrix.push_back(tempvec);
  }

  for (int i = 0; i < col; i++)
    for (int j = 0; j < row; j++)
      transpose_matrix[i][j] = matrix[j][i];

  return transpose_matrix;
}

double error_calculation(vector<vector<double>> &mat1,vector<vector<double>> &mat2)
{
  int i,j;

  double error_square;

  for(i=0;i<mat1.size();i++)
   {
     for(j=0;j<mat1[0].size();j++)
     {
       error_square += ((mat1[i][j] - mat2[i][j]) * (mat1[i][j] - mat2[i][j]));
     }

   }
   return sqrt(error_square);
}

void print(int row, int col, vector<vector<double>> &matrix)
{
  for (int i = 0; i < row; i++)
  {
    for (int j = 0; j < col; j++)
    {
      cout << matrix[i][j] << " ";
    }
    cout << endl;
  }

  cout << endl;
}

vector<vector<double>> add_elementwise(vector<vector<double>> &mat1,vector<vector<double>> &mat2,int row,int col)
{
    int i , j;
    vector<vector<double>> added_matrix;

  for (int i = 0; i < row; i++)
  {
    vector<double> tempvec;
    for (int j = 0; j < col; j++)
    {
      tempvec.push_back(0);
    }
    added_matrix.push_back(tempvec);
  }

  for(i=0;i<row;i++)
   {
     for(j=0;j<col;j++)
     {
       added_matrix[i][j] = mat1[i][j] + mat2[i][j] ;
     }
   }

   return added_matrix;

}

vector<vector<double>> multiply_elementwise(vector<vector<double>> &mat1,vector<vector<double>> &mat2,int row,int col)
{
    int i , j;
    vector<vector<double>> mul_matrix;

  for (int i = 0; i < row; i++)
  {
    vector<double> tempvec;
    for (int j = 0; j < col; j++)
    {
      tempvec.push_back(0);
    }
    mul_matrix.push_back(tempvec);
  }

  for(i=0;i<row;i++)
   {
     for(j=0;j<col;j++)
     {
       mul_matrix[i][j] = mat1[i][j] * mat2[i][j] ;
     }
   }

   return mul_matrix;

}

vector<vector<double>> divide_elementwise(vector<vector<double>> &mat1,vector<vector<double>> &mat2,int row,int col)
{
    int i , j;
    vector<vector<double>> div_matrix;

  for (int i = 0; i < row; i++)
  {
    vector<double> tempvec;
    for (int j = 0; j < col; j++)
    {
      tempvec.push_back(0);
    }
    div_matrix.push_back(tempvec);
  }

  for(i=0;i<row;i++)
   {
     for(j=0;j<col;j++)
     {
       div_matrix[i][j] = mat1[i][j] / mat2[i][j] ;
     }
   }

   return div_matrix;

}

vector<vector<double>> constant_multiply(vector<vector<double>> &mat,int row,int col,int constant)
{
    int i , j;

  for(i=0;i<row;i++)
   {
     for(j=0;j<col;j++)
     {
       mat[i][j] = constant * mat[i][j];
     }
   }

   return mat;

}

vector<vector<double>> constant_add(vector<vector<double>> &mat,int row,int col,double constant)
{
    int i , j;

  for(i=0;i<row;i++)
   {
     for(j=0;j<col;j++)
     {
       mat[i][j] = constant + mat[i][j];
     }
   }

   return mat;

}




vector<vector<double>> Adjacency_matrix(vector<vector<double>> &original_matrix,int row,int col)
{
   double dis = 0.0;

   int i , j ,k;

   vector<vector<double>> adjacency_matrix;

   for(i=0;i<row;i++)
  {
    vector<double> temp;
    for(j=0;j<row;j++)
    {
      temp.push_back(0);
    }
    adjacency_matrix.push_back(temp);
  }

   for(k=0;k<row-1;k++)
   {
     for(i=1+k;i<row;i++)
    {
        for(j=0;j<col;j++)
         {
           dis += pow((original_matrix[k][j] - original_matrix[i][j]),2);
           adjacency_matrix[k][i] = adjacency_matrix[i][k] = sqrt(dis);
         }
         //cout << sqrt(dis) << " ";
         dis = 0.0;
    }
    cout << endl;
   }

   return adjacency_matrix;
}

vector<vector<double>> Degree_matrix(vector<vector<double>> &adjacency_matrix,int row)
{

   int i , j;

   vector<vector<double>> degree_matrix;

    for(i=0;i<row;i++)
   {
    vector<double> temp;
    for(j=0;j<row;j++)
    {
      temp.push_back(0);
    }
    degree_matrix.push_back(temp);
   }

   for(i=0;i<row;i++)
    {
     for(j=0;j<row;j++)
      {
        degree_matrix[i][i] += adjacency_matrix[i][j];
      }
    }

    return degree_matrix;
}

vector<vector<double>> Z_calculation(vector<vector<double>> &X,vector<vector<double>> &UVt ,int N , int M)
{
     vector<vector<double>> Z;

     vector<double> D;

     for(int i=0;i<N;i++)
     {
         vector<double> tempvec;
         for(int j=0;j<M;j++)
         {
             tempvec.push_back(0);
         }
         Z.push_back(tempvec);
     }

     for(int i=0;i<N;i++)
     {
         double dis = 0.0;

         for(int j=0;j<M;j++)
         {
             dis += pow((X[i][j] - UVt[i][j]),2);
         }

         D[i] = sqrt(dis);
     }

     for(int i=0;i<N;i++)
     {
         Z[i][i] = 1/D[i];
     }

     return Z;
}

double H_calculation(vector<vector<double>> &mat,int row , int col)
{
   vector<vector<double>> V = mat;

   double H = 0.0;

   for(int i=0;i<row;i++)
     {
         vector<double> tempvec;
         for(int j=0;j<col;j++)
         {
             tempvec.push_back(0);
         }
         V.push_back(tempvec);
     }

   for(int i=0;i<row;i++)
   {
     for(int j=0;j<col;j++)
     {
        H += (V[i][j]*V[i][j]);
     }
   }
   H = sqrt(H);

   return H;
}

vector<vector<double>> E_matrix(int row)
{
  int N = row;

  vector<vector<double>> E;

  for(int i=0;i<N;i++)
     {
         vector<double> tempvec;
         for(int j=0;j<N;j++)
         {
             tempvec.push_back(0);
         }
         E.push_back(tempvec);
     }

   for(int i=0;i<N;i++)
     {
         for(int j=0;j<N;j++)
         {
             E[i][j] = 1/(double)N;
         }
         
     }
  
  return E;
}


