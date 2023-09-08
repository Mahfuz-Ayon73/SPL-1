#include <iostream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <random>
#include<cmath>

#define MAX_ITERATION 2000

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

vector<vector<double>> update_matrix_U(vector<vector<double>> &mat_A,vector<vector<double>> &mat_B,vector<vector<double>> &original_mat,int M ,int K , int N)
{
   vector<vector<double>> X = original_mat , U = mat_A, V = mat_B, W = Adjacency_matrix(original_mat,M,N), L , D = Degree_matrix(W,M) , Vt = transpose(N,K,V);

   vector<vector<double>> numerator1 , numerator2 , numerator , denominator1 , denominator2 , denominator ,ratio ,updated_U;

   //numrator operations
   numerator1 = matrixMultiplication(M,N,K,X,V);

   numerator2 = matrixMultiplication(M,M,K,W,U);

   numerator2 = constant_multiply(numerator2,M,K,Lamda);

   numerator = add_elementwise(numerator1,numerator2,M,K);

   //print(M,K,numerator);

   //denominator operations
   denominator1 = matrixMultiplication(M,K,N,U,Vt);

   denominator1 = matrixMultiplication(M,N,K,denominator1,V);

   denominator2 = matrixMultiplication(M,M,K,D,U);

   denominator2 = constant_multiply(denominator2,M,K,Lamda);

   denominator = add_elementwise(denominator1,denominator2,M,K);

   denominator = constant_add(denominator,M,K,EPS);

   //print(M,K,denominator);

   ratio = divide_elementwise(numerator,denominator,M,K);

   //print(M,K,ratio);

   updated_U = multiply_elementwise(U,ratio,M,K);

   //print(M,K,update_U);

   return updated_U;
   
   }

vector<vector<double>> update_matrix_V(vector<vector<double>> &mat_A, vector<vector<double>> &mat_B, vector<vector<double>> &original_matrix,int M,int K,int N)
{
   vector<vector<double>> Xt = transpose(M,N,original_matrix), U = mat_A , Ut = transpose(M,K,U) , V = mat_B ;

   vector<vector<double>> numerator  , denominator , ratio , updated_V;  

   numerator = matrixMultiplication(N,M,K,Xt,U);

   //print(N,K,numerator);

   denominator = matrixMultiplication(N,K,M,V,Ut);

   denominator = matrixMultiplication(N,M,K,denominator,U);

   denominator = constant_add(denominator,N,K,EPS);

   //print(N,K,denominator);

   ratio = divide_elementwise(numerator,denominator,N,K);

  // print(N,K,ratio);

   updated_V= multiply_elementwise(V,ratio,N,K);

 //  print(N,K,updated_V);

   return updated_V;
   
}


int main()
{
   freopen("GNMF_input.txt","r",stdin);

   srand(time(0));

   int M , N , K;

   int i,j,k;

   cin >> M >> N >> K;

   vector<vector<double>> X, U, V ,update_U , update_V ,update_X , adj , deg ;

   double error = 0.0 , min = INT32_MAX;

   for(i=0;i<M;i++)
   {
     vector<double> vec;
     for(j=0;j<N;j++)
     {
       double temp;
       cin >> temp;
       vec.push_back(temp);
     }
     X.push_back(vec);
   }


  // initializing A and B randomly
  for (i = 0; i < M; i++)
  {

    vector<double> tempvec1;
    for (j = 0; j < K; j++)
    {
      double temp = rand() % 20; // unif(re);
      tempvec1.push_back(temp);
    }
    U.push_back(tempvec1);
    update_U.push_back(tempvec1);
  }

  for (i = 0; i < N; i++)
  {
    vector<double> tempvec2;
    for (j = 0; j < K; j++)
    {
      double temp = rand() % 20; // unif(re);
      tempvec2.push_back(temp);
    }
    V.push_back(tempvec2);
    update_V.push_back(tempvec2);
  }
 
   cout << "matrix X" << endl;

   print(M,N,X);

   cout << "matrix U" << endl;

   print(M,K,U);

   cout << "matrix V" << endl;

   print(N,K,V);

 
    for(i=1;i<=MAX_ITERATION;i++)
  {

     if(i % 2 != 0)
   {
       U = update_matrix_U(U,V,X,M,K,N);

       V = transpose(N,K,V);
       
       cout << "After Update U and V" << endl;
       
       print(M, K, U);
       print(K, N, V);
  
       update_X = matrixMultiplication(M,K,N,U,V);

       print(M,N,update_X);
  }
 
    else{
      
      V = transpose(K,N,V);

      V = update_matrix_V(U,V,X,M,K,N);
   
      V = transpose(K,N,V);

      cout << "After Update A and B" << endl;
      print(M, K, U);
      print(K,N,V);
  
      update_X = matrixMultiplication(M,K,N,U,V);

      print(M,N,update_X);
   
  }
  
   error = error_calculation(X,update_X);
   cout << "Error:" << error << endl;
    if(error <= min)
      min = error;
  }

  cout << "Factorisation complete!"<<endl;

  cout << "total iteration:" << i-1 << endl << "Minimum Error: " << min << endl;

  

   
}