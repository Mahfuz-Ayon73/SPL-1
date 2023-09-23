#include<bits/stdc++.h>
#include "Matrix_Basic_Operations.h"

#define lemda 0.5

#define alpha 0.1

#define beta 0.2

#define MAX_ITERATION 500

using namespace std;

vector<vector<double>> update_matrix_U(vector<vector<double>> &original_mat,vector<vector<double>> &mat_U,vector<vector<double>> &mat_V,vector<vector<double>> &adj , vector<vector<double>> &deg , int N,int K,int M)
{
   vector<vector<double>> U = mat_U , V = mat_V , X = original_mat , V_transpose ,W, D;
 
   vector<vector<double>> numerator , denominator ,X_Vt , W_U , D_U ,U_Vt ,U_Vt_V , lemda_WU ,lemda_DU;
  
   W = adj , D = deg;

   X_Vt = matrixMultiplication(M,N,K,X,V);

   W_U = matrixMultiplication(M,M,K,W,U);

   lemda_WU = constant_multiply(W_U,M,K,lemda);

   numerator = add_elementwise(X_Vt,lemda_WU,M,K);

   U_Vt = matrixMultiplication(M,K,N,U,V_transpose);

   U_Vt_V = matrixMultiplication(M,N,K,U_Vt,V);

   
}

vector<vector<double>> update_matrix_V(vector<vector<double>> &original_mat,vector<vector<double>> &mat_U,vector<vector<double>> &mat_V , int N,int K,int M)
{
   vector<vector<double>> U_transpose , U = mat_U , V = mat_V , X = original_mat, Z ,V_transpose , U_transposeX;

   vector<vector<double>> numerator , denominator , U_transpose_U , Ut_U_Vt , Ratio ,updated_V;

   V_transpose = transpose(N,K,V);//K*N

   U_transpose = transpose(M,K,U);//K*M

   U_transposeX = matrixMultiplication(K,M,N,U_transpose,X);

   //Z = Z_calculation(X,);

   numerator = U_transposeX;

   U_transpose_U  = matrixMultiplication(K,M,K,U_transpose,U);

   Ut_U_Vt = matrixMultiplication(K,K,N,U_transpose_U,V_transpose);

   denominator = Ut_U_Vt;

   Ratio = divide_elementwise(numerator,denominator,K,N);

   updated_V = multiply_elementwise(V_transpose,Ratio,K,N);
}
int main()
{
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

   adj = Adjacency_matrix(X,M,N);

   deg = Degree_matrix(X,M);

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

      cout << "After Update U and V" << endl;
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
