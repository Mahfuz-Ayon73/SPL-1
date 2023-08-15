#include <iostream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <random>

#define MAX_ITERATION 100

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
   return error_square;
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

vector<vector<double>> update_matrix_A(vector<vector<double>> &mat_A,vector<vector<double>> &mat_B,vector<vector<double>> &original_mat,int rowA,int colA,int rowB,int colB,int rowX,int colX)
{
   vector<vector<double>> numerator,denominator ,temp ,ratio;
   
   vector<vector<double>> updated_matA = mat_A,trans_B = transpose(rowB,colB,mat_B);

   int i,j,k;

   numerator = matrixMultiplication(rowX,colX,colB,original_mat,trans_B);

   temp = matrixMultiplication(rowB,colB,rowB,mat_B,trans_B);

   denominator = matrixMultiplication(rowA,colA,rowB,mat_A,temp);
   
   for (i = 0; i < rowA; i++)
  {
    vector<double> temp;
    for (j = 0; j < colA; j++)
    {
      temp.push_back(numerator[i][j] / denominator[i][j]);
    }
    ratio.push_back(temp);
  }

  for (i = 0; i <rowA; i++)
  {
    for (j = 0; j < colA; j++)
    {
      updated_matA[i][j] = mat_A[i][j] * ratio[i][j];
    }
  }

  return updated_matA;
}

vector<vector<double>> update_matrix_B(vector<vector<double>> &mat_A,vector<vector<double>> &mat_B,vector<vector<double>> &original_mat,int rowA,int colA,int rowB,int colB,int rowX,int colX)
{
   vector<vector<double>> numerator,denominator ,temp ,ratio;
   
   vector<vector<double>> updated_matB = mat_B,trans_A = transpose(rowA,colA,mat_A);

   int i,j,k;

   numerator = matrixMultiplication(colA,rowA,colX,trans_A,original_mat);

   temp = matrixMultiplication(colA,rowA,colA,trans_A,mat_A);

   denominator = matrixMultiplication(colA,colA,colB,temp,mat_B);
   
   for (i = 0; i < rowB; i++)
  {
    vector<double> temp;
    for (j = 0; j < colB; j++)
    {
      temp.push_back(numerator[i][j] / denominator[i][j]);
    }
    ratio.push_back(temp);
  }

  for (i = 0; i <rowB; i++)
  {
    for (j = 0; j < colB; j++)
    {
      updated_matB[i][j] = mat_B[i][j] * ratio[i][j];
    }
  }

  return updated_matB;
}

int main()
{

   freopen("input.txt","r",stdin); 
   //freopen("output.txt","w",stdout);

  int N = 3, K , M = 3;

  cin >> K;

  int i, j, k;

  srand(time(0));

  vector<vector<double>> X, A, B ,update_A , update_B ,update_X;

  for(i=0;i<N;i++)
  {
    vector<double> temp;
    for(j=0;j<M;j++)
    {
      double a;
      cin >> a;
      temp.push_back(a);
    }
    X.push_back(temp);
  }
  // initializing A and B randomly
  for (i = 0; i < N; i++)
  {

    vector<double> tempvec1;
    for (j = 0; j < K; j++)
    {
      double temp = rand() % 20; // unif(re);
      tempvec1.push_back(temp);
    }
    A.push_back(tempvec1);
    update_A.push_back(tempvec1);
  }

  for (i = 0; i < K; i++)
  {
    vector<double> tempvec2;
    for (j = 0; j < M; j++)
    {
      double temp = rand() % 20; // unif(re);
      tempvec2.push_back(temp);
    }
    B.push_back(tempvec2);
    update_B.push_back(tempvec2);
  }

  cout << "Matrix A:" << endl;
  print(N, K, A);
  cout << "Matrix B:" << endl;
  print(K, M, B);

  vector<vector<double>> mul;
  mul = matrixMultiplication(N, K, M, A, B);
  cout << "A matrix * B matrix:" << endl;
  print(N, M, mul);

 // Iteration
  for(i=0;i<MAX_ITERATION;i++)
  {

  // updating process of A matrix
  //  A = A * (X * Bt) / (A * B * Bt)

  update_A = update_matrix_A(A,B,X,N,K,K,M,N,M);

  cout << "Updated A" << endl;
  print(N, K, update_A);

  // updating B matrix:
  //  B = B * (At * X) / (At * A * B) 

  update_B = update_matrix_B(A,B,X,N,K,K,M,N,M);

  cout << "updated B" << endl;
  print(K, M, update_B);

  //update A * update B
 
   update_X = matrixMultiplication(N,K,M,update_A,update_B);

  cout << "Updated A * B matrix" << endl;
  print(N,M,update_X);

  cout << "Error:" << error_calculation(X,update_X);

  A = update_A , B = update_B;

  cout << endl << endl;
  }
  
 
  return 0;
}