#include<iostream>
#include"Matrix_Basic_Operations.h"

#define MAX_ITERATION 1000

#define alpha 0.1

#define beta 0.2

using namespace std;

vector<vector<double>> update_matrix_U(vector<vector<double>> &original_mat,vector<vector<double>> &mat_U,vector<vector<double>> &mat_V, vector<vector<double>> &adj , vector<vector<double>> &deg , int N,int K,int M)
{
   vector<vector<double>> U = mat_U , V = mat_V , X = original_mat , V_transpose = transpose(K,N,V) , U_transpose = transpose(M,K,U) , W = adj , D  = deg;

   vector<vector<double>> XZ , XZVt , WU , alphaWU , betaU , UV , UVZ , UVZVt , DU , alphaDU , betaUt , E , betaUE ,Ratio , update_U;

   vector<vector<double>> numerator , denominator ;

   double Z = 0.0;

   Z = Z_calculation(X,U,V,M,K,N);

  //  cout << "Z " << endl;
  //  cout << Z << endl;

   XZ = constant_multiply(X,M,N,Z);

  //  cout << "XZ" << endl;
  //  print(M,N,XZ);

   XZVt = matrixMultiplication(M,N,K,XZ,V_transpose);

  //  cout << "XZVt" << endl;
  //  print(M,K,XZVt);

   WU = matrixMultiplication(M,M,K,W,U);

   alphaWU = constant_multiply(WU,M,K,alpha);

   betaU = constant_multiply(U,M,K,beta);

   numerator = add_elementwise(XZVt,alphaWU,M,K);//XZV + alphaWU

  //  cout << "numerator" << endl;
  //  print(M,K,numerator);

   numerator = add_elementwise(numerator,betaU,M,K);//XZV + alphaWU + betaU

  //  cout << "numerator" << endl;
  //  print(M,K,numerator);

   UV = matrixMultiplication(M,K,N,U,V);

   UVZ = constant_multiply(UV,M,N,Z);

   UVZVt  = matrixMultiplication(M,N,K,UVZ,V_transpose);

   DU = matrixMultiplication(M,M,K,D,U);

   alphaDU = constant_multiply(DU,M,K,alpha);

   betaU = constant_multiply(U,M,K,beta);

   E = E_matrix(K);

  //  cout << "E "<<endl;
  //  print(K,K,E);

   betaUE = matrixMultiplication(M,K,K,betaU,E);

  //  cout << "betaUE" << endl;
  //  print(M,K,betaUE);

   denominator = add_elementwise(UVZVt,alphaDU,M,K);

   denominator = add_elementwise(denominator,betaUE,M,K);

  //  cout << "denominator" << endl;
  //  print(M,K,denominator);

   denominator = constant_add(denominator,M,K,EPS);

  //  cout << "denominator" << endl;
  //  print(M,K,denominator);

   Ratio = divide_elementwise(numerator,denominator,M,K);
 
   update_U = multiply_elementwise(U,Ratio,M,K);

   return update_U;
   
}

vector<vector<double>> update_matrix_V(vector<vector<double>> &original_mat,vector<vector<double>> &mat_U,vector<vector<double>> &mat_V, int N,int K,int M)
{
  vector<vector<double>> U = mat_U , V = mat_V , X = original_mat , U_transpose = transpose(M,K,U);

  vector<vector<double>> numerator , denominator ,UtX , UtXZ , UtU , UtUV , UtUVZ ,VH ,  lemdaVH , Ratio , updated_V;

  double Z = 0.0 , H = 0.0;

  UtX = matrixMultiplication(K,M,N,U_transpose,X);

  // cout << "UtX" << endl;
  // print(K,N,UtX);

  Z = Z_calculation(X,U,V,M,K,N);

  // cout << "Z" << endl;
  // cout << Z << endl;

  UtXZ = constant_multiply(UtX,K,N,Z);

  // cout << "UtXZ" << endl;
  // print(K,N,UtXZ);

  numerator = UtXZ;

  // cout << "numeratorV" << endl;
  // print(K,N,numerator);

  UtU = matrixMultiplication(K,M,K,U_transpose,U);

  // cout << "UtU" << endl;
  // print(K,K,UtU);

  UtUV = matrixMultiplication(K,K,N,UtU,V);

  UtUVZ = constant_multiply(UtUV,K,N,Z);

  H = H_calculation(V,K,N);

  VH = constant_multiply(V,K,N,H);

  lemdaVH = constant_multiply(VH,K,N,Lamda);

  denominator = add_elementwise(UtUVZ,lemdaVH,K,N);

  // cout << "denominatorV" << endl;
  // print(K,N,denominator);

  denominator = constant_add(denominator,K,N,EPS);

  // cout << "denominatorV" << endl;
  // print(K,N,denominator);

  Ratio = divide_elementwise(numerator,denominator,K,N);

  updated_V = multiply_elementwise(V,Ratio,K,N);

  
  //print(K,N,updated_V);

  return updated_V;

}
void RSNMF()
{
   freopen("input.txt","r",stdin);

   int M , N , K;

   int i , j , k;

   cout << "Enter dimensions:M , N , K \n" ;
   cin >> M >> N >> K;

   vector<vector<double>> X, U, V ,Vt,update_U , update_V ,update_X , adj , deg ;

   double error = 0.0 , min = INT32_MAX ,current_error = 0.0;

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
      double temp = rand() % 20;
      tempvec1.push_back(temp);
    }
    U.push_back(tempvec1);
    update_U.push_back(tempvec1);
  }

  for (i = 0; i < K; i++)
  {
    vector<double> tempvec2;
    for (j = 0; j < N; j++)
    {
      double temp = rand() % 20;
      tempvec2.push_back(temp);
    }
    V.push_back(tempvec2);
    update_V.push_back(tempvec2);
  }

   //cout << "matrix X" << endl;

   //print(M,N,X);

  // cout << "matrix U" << endl;

   //print(M,K,U);

   //cout << "matrix V" << endl;

   //print(K,N,V);

   update_V = update_matrix_V(X,U,V,N,K,M);

  // print(K,N,update_V);

   update_U = update_matrix_U(X,U,V,adj,deg,N,K,M);

   //print(M,K,update_U);

   update_X = matrixMultiplication(M,K,N,U,V);

   current_error = error_calculation(X,update_X);

   for(i=0;i<MAX_ITERATION;i++)
  {

  // updating process of A matrix
  //  A = A * (X * Bt) / (A * B * Bt)

  if(i % 2 == 0)
  {
   update_U = update_matrix_U(X,U,V,adj,deg,N,K,M);

   U = update_U;

   //cout << "After Update A and B" << endl;
   //print(M, K, U);
   //print(K, N, V);

   update_X = matrixMultiplication(M,K,N,U,V);


  }


  // updating B matrix:
  //  B = B * (At * X) / (At * A * B)
  else{
   update_V = update_matrix_V(X,U,V,N,K,M);

   V = update_V;

   //cout << "After Update A and B" << endl;
   //print(M, K, U);
   //print(K, N, V);

   update_X = matrixMultiplication(M,K,N,U,V);


  }

   error = error_calculation(X,update_X);
   cout << "Error:" << error << endl;

  if(error <= 1e-6)
    break;
  else if(error > current_error)
  {
    break;
  }

  current_error = error;
   
  }

  cout << "Factorisation complete!"<<endl;

  cout << "total iteration:" << i << endl;

  cout << "Matrix U and V are:" << endl;

  print(M, K, U);
  print(K, N, V);

  cout << "U * V :" << endl;

  print(M,N,update_X);

 
}
