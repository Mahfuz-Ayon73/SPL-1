#include<bits/stdc++.h>
#include "Matrix_Basic_Operations.h"

#define lemda 0.01

#define alpha 0.1

#define beta 0.2

#define MAX_ITERATION 100

using namespace std;

vector<vector<double>> update_matrix_U(vector<vector<double>> &original_mat,vector<vector<double>> &mat_U,vector<vector<double>> &mat_V,vector<vector<double>> &adj , vector<vector<double>> &deg , int N,int K,int M)
{
   vector<vector<double>> U = mat_U , V = mat_V , X = original_mat , V_transpose ,W, D;
 
   vector<vector<double>> numerator , denominator ,X_V , W_U , D_U ,U_Vt ,U_Vt_V , lemda_WU ,lemda_DU , Ratio , updated_U;
  
   W = adj , D = deg;

   V_transpose  = transpose(N,K,V);

   X_V = matrixMultiplication(M,N,K,X,V);//M*K

  // print(M,K,X_V);


   W_U = matrixMultiplication(M,M,K,W,U);//M*K

   lemda_WU = constant_multiply(W_U,M,K,lemda);

   numerator = add_elementwise(X_V,lemda_WU,M,K);//M*K

   //print(M,K,numerator);


   U_Vt = matrixMultiplication(M,K,N,U,V_transpose);

   cout << "X " << endl;

   //print(M,N,U_Vt);

   U_Vt_V = matrixMultiplication(M,N,K,U_Vt,V);

   //print(M,K,U_Vt_V);

   D_U = matrixMultiplication(M,M,K,D,U);

   lemda_DU = constant_multiply(D_U,M,K,lemda);

   //print(M,K,lemda_DU);

   denominator = add_elementwise(U_Vt_V,lemda_DU,M,K);

   //print(M,K,denominator);

   Ratio = divide_elementwise(numerator,denominator,M,K);

   //print(M,K,Ratio);

   updated_U = multiply_elementwise(U,Ratio,M,K);

  // print(M,K,updated_U);

   return updated_U;
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

   //print(K,M,numerator);

   U_transpose_U  = matrixMultiplication(K,M,K,U_transpose,U);

   Ut_U_Vt = matrixMultiplication(K,K,N,U_transpose_U,V_transpose);

   denominator = Ut_U_Vt;

   //print(K,M,denominator);

   Ratio = divide_elementwise(numerator,denominator,K,N);

   updated_V = multiply_elementwise(V_transpose,Ratio,K,N);

   return updated_V;
}
void RNMF()
{
   freopen("RNMF_input.txt","r",stdin);
   int M , N , K;

   int i,j,k;

   cout << "Enter dimensions:M , N , K\n";
   cin >> M >> N >> K;

   vector<vector<double>> X, U, V ,Vt,update_U , update_V ,update_X , adj , deg ;

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

  // initializing U and V randomly
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
        U = update_matrix_U(X,U,V,adj,deg,N,K,M);

        V = transpose(N,K,V);

        update_X = matrixMultiplication(M,K,N,U,V);

        print(M,N,update_X);

        V = transpose(K,N,V);

     }
     else
     {

        V = update_matrix_V(X,U,V,N,K,M);

        update_X = matrixMultiplication(M,K,N,U,V);

         print(M,N,update_X);
        
        V = transpose(K,N,V);
     }

   }
}
