#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "Los.h"

using namespace std;


void LUsq(vector<double> di,
   vector<double> a_l,
   vector<double> a_u,
   vector<int> ig,
   vector<int> jg,
   vector<double>& q,
   vector<double> B,
   int N)
{

   vector<double> di_f;
   vector<double> ggl_f;
   vector<double> ggu_f;

   di_f = di;
   ggl_f = a_l;
   ggu_f = a_u;

   for (int i = 0; i < N; i++)
   {
      double sumdi = 0.0;

      int i0 = ig[i];
      int i1 = ig[i + 1];


      for (int k = i0; k < i1; k++)
      {
         int j = jg[k];
         int j0 = ig[j];

         int j1 = ig[j + 1];


         int ik = i0;
         int kj = j0;

         double suml = 0.0;
         double sumu = 0.0;

         while (ik < k)
         {

            if (jg[ik] == jg[kj])
            {

               suml += ggl_f[ik] * ggu_f[kj];
               sumu += ggu_f[ik] * ggl_f[kj];
               ik++;
               kj++;
            }

            else
               jg[ik] > jg[kj] ? kj++ : ik++;
         }

         ggl_f[k] = (ggl_f[k] - suml) / di_f[j];
         ggu_f[k] = (ggu_f[k] - sumu) / di_f[j];
         sumdi += ggl_f[k] * ggu_f[k];
      }

      di_f[i] = sqrt(di_f[i] - sumdi);
   }
   ofstream Test;
   Test.open("Test.txt", std::ios_base::app);
   Test << "\n Factorai\nvar Di = new[] {";
   for (size_t i = 0; i < di_f.size(); i++)
   {
      Test << di_f[i] << ", ";
   }
   Test << "};\nvar Al = new[] {";
   for (size_t i = 0; i < ggl_f.size(); i++)
   {
      Test << ggl_f[i] << ", ";
   }
   Test << "};\nvar Au = new[] {";
   for (size_t i = 0; i < ggu_f.size(); i++)
   {
      Test << ggu_f[i] << ", ";
   }
   Test << "};\nvar Au = new[] {";
   for (size_t i = 0; i < ggu_f.size(); i++)
   {
      Test << ggu_f[i] << ", ";
   }
   Test << "};\n";
   Test.close();
   LoS_precond(di_f, ggu_f, ggl_f, di, a_u, a_l, ig, jg, q, B, N);
}

void LoS_precond(vector<double>& dif, vector<double>& a_u_f,  vector<double>& a_l_f, vector<double>& di, vector<double>& a_u, vector<double>& a_l, vector<int>& ig, vector<int>& jg, vector<double>& q, vector<double> B, int N)
{
   int k = 0;
   vector<double> buf = Mult(q, di, a_u,a_l, ig,jg);

   //ofstream Test;
   //Test.open("Vectora.txt", std::ios_base::app);
   //Test << "\n buf = new[] {";
   //for (size_t i = 0; i < buf.size(); i++)
   //{
   //   Test << buf[i] << ", ";
   //}
   //Test << "};\n";

   for (int i = 0; i < N; i++)
   {
      buf[i] = B[i] - buf[i];
   }

   /*Test << "\n buf_b_buf = new[] {";
   for (size_t i = 0; i < buf.size(); i++)
   {
      Test << buf[i] << ", ";
   }
   Test << "};\n";*/
   vector<double> r = LUDirect(buf, dif, a_l_f,ig, jg);
   /*Test << "\n r = new[] {";
   for (size_t i = 0; i < r.size(); i++)
   {
      Test << r[i] << ", ";
   }
   Test << "};\n";*/
   double error = scalar_prod(r, r);
   double error_ = 0;
   vector<double> z = LUReverse(r, dif, a_u_f, ig, jg, N);
   /*Test << "\n z = new[] {";
   for (size_t i = 0; i < z.size(); i++)
   {
      Test << z[i] << ", ";
   }
   Test << "};\n";*/
   buf = Mult(z, di, a_u,a_l, ig, jg);
  /* Test << "\n buf_z = new[] {";
   for (size_t i = 0; i < buf.size(); i++)
   {
      Test << buf[i] << ", ";
   }
   Test << "};\n";*/
   vector<double> p = LUDirect(buf, dif, a_l_f, ig,jg);
   /*Test << "\n p = new[] {";
   for (size_t i = 0; i < p.size(); i++)
   {
      Test << p[i] << ", ";
   }
   Test << "};\n";
   Test.close();*/
   while (error > 1e-15 && k < 1000 && abs(error_ - error)>=1e-15)
   {
      double pp = scalar_prod(p, p);
      double pr = scalar_prod(p, r);
      double alpha = pr / pp;
      error_ = error;
      error -= alpha * alpha * pp;
      for (int i = 0; i < N; i++)
      {
         q[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }
      vector<double> Ur = LUReverse(r, dif, a_u_f,ig,jg,N);
      buf = Mult(Ur, di, a_u,a_l, ig, jg);
      buf = LUDirect(buf, dif, a_l_f,ig,jg);
      double betta = -(scalar_prod(p, buf) / pp);
      for (int i = 0; i < N; i++)
      {
         z[i] = Ur[i] + betta * z[i];
         p[i] = buf[i] + betta * p[i];
      }
      k++;
   }
   cout << "k:" << k << "\nerror: " << error << "\n";
  // print_(dif, a_u_f, a_l_f,ig,jg,B,N);
};

vector<double> Mult(vector<double>& v, vector<double> di, vector<double> a_u, vector<double> a_l, vector<int> ig, vector<int> jg)
{
   vector<double> res(v.size());
   for (int i = 0; i < v.size(); i++)
   {
      res[i] = di[i] * v[i];
      for (int j = ig[i]; j < ig[i + 1]; j++)
      {
         res[i] += a_l[j] * v[jg[j]];
         res[jg[j]] += a_u[j] * v[i];
      }
   }
   return res;
}

vector<double> LUDirect(const  vector<double>& b, vector<double>& _dif, vector<double>& _alf, vector<int>& ig, vector<int>& jg)
{
   vector<double> res = b;

   for (size_t i = 0; i < res.size(); i++)
   {
      double sum = 0.0;
      for (size_t j = ig[i]; j < ig[i + 1]; j++)
         sum += _alf[j] * res[jg[j]];
      res[i] -= sum;
      res[i] /= _dif[i];
   }
   return res;
}

vector<double> LUReverse(const  vector<double>& b, vector<double>& _dif, vector<double>& _alf, vector<int>& ig, vector<int>& jg, int N)
{
   vector<double> res = b;

   for (int i = N - 1; i >= 0; i--)
   {
      res[i] /= _dif[i];
      for (size_t j = ig[i]; j < ig[i + 1]; j++)
         res[jg[j]] -= _alf[j] * res[i];
   }
   return res;
}

double scalar_prod(vector<double>& x, vector<double>& y)
{
   double res = 0.0;
   if (x.size() == y.size())
   {
      for (int i = 0; i < y.size(); i++)
      {
         res += x[i] * y[i];
      }
      return res;
   }
   else
   {
      cout << "Error!";
      return res;
   }
}

void print_(vector<double>& di, vector<double>& ggu, vector<double>& ggl, vector<int>& ig, vector<int>& jg, vector<double>& B, int N)
{
   vector < vector<double>> mat;
   mat.resize(N);
   for (int i = 0; i < mat.size(); i++)
   {
      mat[i].resize(N);
   }

   for (int i = 0; i < mat.size(); i++)
   {
      mat[i][i] = di[i];
      for (int j = ig[i]; j < ig[i + 1]; j++)
      {
         mat[i][jg[j]] = ggl[j];
         mat[jg[j]][i] = ggu[j];
      }
   }
   ofstream matrix;
   matrix.open("matrix.txt", std::ios_base::app);
   matrix << endl;
   for (int i = 0; i < mat.size(); i++)
   {
      for (int j = 0; j < mat[i].size(); j++)
         matrix << setw(12) << mat[i][j] << " ";
      matrix <<"    |    " << B[i] << endl;
   }
   matrix.close();
}