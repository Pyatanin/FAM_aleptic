#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "parmetrs.h"
#include "Los.h"


using namespace std;
typedef double type;

class MKA
{
private:
   int n_r; // количество точек по r
   int n_z;// количество точек по z
   vector<int> ig; 
   vector<int> jg;
   vector<type> al; // верхняя диагональная
   vector<type> au;// нижняя диагональная
   vector<type> di; 
   vector<vector<type>> loc_m; // локальная матрица масс
   vector<vector< vector<type>>> loc_g;// локальная матрица жесткости
   vector<vector<type>> grid;// сетка
   vector<vector< vector<int>>> el_pints;// связи в элементе фиктивной нумерации и обычной
   vector<vector<int>> connection; // связи точек
   vector<type> B;// вектор правой части
   vector<type> q; // решение
public:
   MKA()
   {
      this->n_r = 0;
      this->n_z = 0;
      this->ig = vector<int>();
      this->jg = vector<int>();
      this->al = vector<type>();
      this->au = vector<type>();
      this->di = vector<type>();
      this->loc_m = vector< vector<type>>();
      this->loc_g = vector< vector< vector<type>>>();
      this->grid = vector< vector<type>>();
      this->el_pints = vector< vector< vector<int>>>();
      this->connection = vector< vector<int>>();
      this->B = vector<type>();
      this->q = vector<type>();
   }
   MKA(int nx,
      int ny,
      vector<int> ig_,
      vector<int> jg_,
      vector<type> al_,
      vector<type> au_,

      vector<type> di_,
      vector< vector<type>> LocMass,
      vector< vector< vector<type>>> LocG,
      vector< vector<type>> Grid,
      vector< vector< vector<int>>> ElBonds,
      vector< vector<int>> NodeBonds,
      vector<type> GlobB,
      vector<type> Q
   )
   {
      n_r = nx;
      n_z = ny;
      ig = ig_;
      jg = jg_;
      al = al_;
      al = au_;
      di = di_;
      loc_m = LocMass;
      loc_g = LocG;
      grid = Grid;
      el_pints = ElBonds;
      connection = NodeBonds;
      B = GlobB;
      q = Q;
   }
   // создает сетку по шагам и связи элементов
   void makingGrid()
   {
      fstream input;
      input.open("input.txt");
      type xmin, xmax, ymin, ymax;
      input >> xmin;
      input >> xmax;
      input >> ymin;
      input >> ymax;
      input >> n_r;
      input >> n_z;
      type hr = (xmax - xmin) / (n_r - 1);
      type hz = (ymax - ymin) / (n_z - 1);
      q.resize((2 * n_r - 1) * (2 * n_z - 1));
      for (int i=0;i<q.size();i++)
      {
         q[i] = 1;
      }

      grid.resize(n_r * n_z);
      int k = 0;
      for (int i = 0; i < n_z; i++)
         for (int j = 0; j < n_r; j++)
         {
            grid[k].resize(2);
            grid[k][0] = xmin + j * hr;
            grid[k][1] = ymin + i * hz;
            k++;
         }
      input.close();
      el_pints.resize(2);
      el_pints[0].resize((n_r - 1) * (n_z - 1));
      el_pints[1].resize((n_r - 1) * (n_z - 1));
      for (int i = 0; i < (n_z - 1); i++)
      {
         for (int j = 0; j < (n_r - 1); j++)
         {
            el_pints[0][i * (n_r - 1) + j].resize(9);
            el_pints[0][i * (n_r - 1) + j][0] = 2 * i * (2 * n_r - 1) + j * 2;
            el_pints[0][i * (n_r - 1) + j][1] = 2 * i * (2 * n_r - 1) + j * 2 + 1;
            el_pints[0][i * (n_r - 1) + j][2] = 2 * i * (2 * n_r - 1) + j * 2 + 2;
            el_pints[0][i * (n_r - 1) + j][3] = 2 * i * (2 * n_r - 1) + j * 2 + 2 * n_r - 1;
            el_pints[0][i * (n_r - 1) + j][4] = 2 * i * (2 * n_r - 1) + j * 2 + 2 * n_r;
            el_pints[0][i * (n_r - 1) + j][5] = 2 * i * (2 * n_r - 1) + j * 2 + 2 * n_r + 1;
            el_pints[0][i * (n_r - 1) + j][6] = 2 * (i + 1) * (2 * n_r - 1) + j * 2;
            el_pints[0][i * (n_r - 1) + j][7] = 2 * (i + 1) * (2 * n_r - 1) + j * 2 + 1;
            el_pints[0][i * (n_r - 1) + j][8] = 2 * (i + 1) * (2 * n_r - 1) + j * 2 + 2;
            el_pints[1][i * (n_r - 1) + j].resize(4);
            el_pints[1][i * (n_r - 1) + j][0] = i * n_r + j;
            el_pints[1][i * (n_r - 1) + j][1] = i * n_r + j + 1;
            el_pints[1][i * (n_r - 1) + j][2] = (i + 1) * n_r + j;
            el_pints[1][i * (n_r - 1) + j][3] = (i + 1) * n_r + j + 1;
         }
      }

      connection.resize((2 * n_r - 1) * (2 * n_z - 1));
      for (int i = 0; i < (n_r - 1) * (n_z - 1); i++)
      {
         for (int j = 0; j < 9; j++)
         {
            for (int k = 0; k < 9; k++)
            {
               if (el_pints[0][i][j] > el_pints[0][i][k])
               {
                  connection[el_pints[0][i][j]].push_back(el_pints[0][i][k]);
               }
            }
         }
      }

      for (int i = 0; i < ((2 * n_r - 1) * (2 * n_z - 1)); i++)
      {
         sort(connection[i].begin(), connection[i].end());
         auto last =  unique(connection[i].begin(), connection[i].end());
         connection[i].erase(last, connection[i].end());
      }
   }
   // создает профиль используя связи элементов
   void profile()
   {
      ig.resize((2 * n_r - 1) * (2 * n_z - 1) + 1);
      di.resize((2 * n_r - 1) * (2 * n_z - 1));
      ig[0] = 0;
      ig[1] = 0;
      for (int i = 0; i < (2 * n_r - 1) * (2 * n_z - 1); i++)
      {
         int k = 0;
         for (int j = 0; j < connection[i].size(); j++)
         {
            k++;
            jg.push_back(connection[i][j]);
         }
         ig[i + 1] = ig[i] + k;
      }
   }
   // Получает локальные матрицы и вставляет на нужные места
   void global_build()
   {
      type gamma = 0;
       vector<type> LocB =  vector<type>(9);
       vector< vector<type>> LocG;
      B.resize((2 * n_r - 1) * (2 * n_z - 1));
      al.resize(ig[(2 * n_r - 1) * (2 * n_z - 1)]);
      for (int i = 0; i < (n_r - 1) * (n_z - 1); i++)
      {
         type rk = grid[el_pints[1][i][0]][0];
         type hr =  grid[el_pints[1][i][1]][0] - rk;
         type hz = grid[el_pints[1][i][2]][1]-grid[el_pints[1][i][0]][1];

         loc_build(rk, hr, hz);
         loc_gamma(i, gamma);
         loc_f(LocB, i);
         m_mult_v(LocB);
         for (int k = 0; k < 9; k++)
         {
            B[el_pints[0][i][k]] += LocB[k];
         }
         loc_G(i, LocG);
         for (int k = 0; k < 9; k++)
         {
            di[el_pints[0][i][k]] += gamma * loc_m[k][k] + LocG[k][k];
         }

         int ii = 0;
         for (int k = 1; k < 9; k++)
         {
            for (int j = 0; j < k; j++)
            {
               for (ii = ig[el_pints[0][i][k]]; jg[ii] != el_pints[0][i][j];)
               {
                  ii++;
               }
               al[ii] += gamma * loc_m[k][j] + LocG[k][j];
            }
         }
      }
      au = al;
   }
   //Делает локальну. матрицу основываясь на шагах и rk
   void loc_build(type rk, type hr, type hz)
   {
       vector< vector< vector<type>>> Loc_M_Old =  vector< vector< vector<type>>>(3);

      Loc_M_Old[0].resize(3);
      // старая
      Loc_M_Old[0][0].resize(3);
      Loc_M_Old[0][1].resize(3);
      Loc_M_Old[0][2].resize(3);
      Loc_M_Old[0][0][0] = 2.0 / 15.0;
      Loc_M_Old[0][0][1] = 1.0 / 15.0;
      Loc_M_Old[0][0][2] = -1.0 / 30.0;
      Loc_M_Old[0][1][0] = 1.0 / 15.0;
      Loc_M_Old[0][1][1] = 8.0 / 15.0;
      Loc_M_Old[0][1][2] = 1.0 / 15.0;
      Loc_M_Old[0][2][0] = -1.0 / 30.0;
      Loc_M_Old[0][2][1] = 1.0 / 15.0;
      Loc_M_Old[0][2][2] = 2.0 / 15.0;

      //x
      Loc_M_Old[1].resize(3);
      Loc_M_Old[1][0].resize(3);
      Loc_M_Old[1][1].resize(3);
      Loc_M_Old[1][2].resize(3);
      Loc_M_Old[1][0][0] = 1.0 / 60.0;
      Loc_M_Old[1][0][1] = 0.;
      Loc_M_Old[1][0][2] = -1.0 / 60.0;
      Loc_M_Old[1][1][0] = 0.;
      Loc_M_Old[1][1][1] = 4.0 / 15.0;
      Loc_M_Old[1][1][2] = 1.0 / 15.0;
      Loc_M_Old[1][2][0] = -1.0 / 60.0;
      Loc_M_Old[1][2][1] = 1.0 / 15.0;
      Loc_M_Old[1][2][2] = 7.0 / 60.0;

      //1-x
      Loc_M_Old[2].resize(3);
      Loc_M_Old[2][0].resize(3);
      Loc_M_Old[2][1].resize(3);
      Loc_M_Old[2][2].resize(3);
      Loc_M_Old[2][0][0] = 7.0 / 60.0;
      Loc_M_Old[2][0][1] = 1.0 / 15.0;
      Loc_M_Old[2][0][2] = -1.0 / 60.0;
      Loc_M_Old[2][1][0] = 1.0 / 15.0;
      Loc_M_Old[2][1][1] = 4.0 / 15.0;
      Loc_M_Old[2][1][2] = 0.;
      Loc_M_Old[2][2][0] = -1.0 / 60.0;
      Loc_M_Old[2][2][1] = 0.;
      Loc_M_Old[2][2][2] = 1.0 / 60.0;

       vector< vector< vector<type>>> Loc_M_New =  vector< vector< vector<type>>>(3);

      // Новая
      Loc_M_New[0].resize(3);
      Loc_M_New[0][0].resize(3);
      Loc_M_New[0][1].resize(3);
      Loc_M_New[0][2].resize(3);
      Loc_M_New[0][0][0] = 1.0 / 60.0;
      Loc_M_New[0][0][1] = 0.;
      Loc_M_New[0][0][2] = -1.0 / 60.0;
      Loc_M_New[0][1][0] = 0.;
      Loc_M_New[0][1][1] = 4.0 / 15.0;
      Loc_M_New[0][1][2] = 1.0 / 15.0;
      Loc_M_New[0][2][0] = -1.0 / 60.0;
      Loc_M_New[0][2][1] = 1.0 / 15.0;
      Loc_M_New[0][2][2] = 7.0 / 60.0;
      //Новая x
      Loc_M_New[1].resize(3);
      Loc_M_New[1][0].resize(3);
      Loc_M_New[1][1].resize(3);
      Loc_M_New[1][2].resize(3);
      Loc_M_New[1][0][0] = 1.0 / 210.0;
      Loc_M_New[1][0][1] = -1./105.;
      Loc_M_New[1][0][2] = -1.0 / 84.0;
      Loc_M_New[1][1][0] = -1. / 105.;
      Loc_M_New[1][1][1] = 16.0 / 105.0;
      Loc_M_New[1][1][2] = 2.0 / 35.0;
      Loc_M_New[1][2][0] = -1.0 / 84.0;
      Loc_M_New[1][2][1] = 2.0 / 35.0;
      Loc_M_New[1][2][2] = 11.0 / 105.0;
      //Новая 1-x
      Loc_M_New[2].resize(3);
      Loc_M_New[2][0].resize(3);
      Loc_M_New[2][1].resize(3);
      Loc_M_New[2][2].resize(3);
      Loc_M_New[2][0][0] = 1.0 / 84.0;
      Loc_M_New[2][0][1] = 1.0 / 105.0;
      Loc_M_New[2][0][2] = -1.0 / 210.0;
      Loc_M_New[2][1][0] = 1.0 / 105.0;
      Loc_M_New[2][1][1] = 23.0 / 140.0;
      Loc_M_New[2][1][2] = -17./420.;
      Loc_M_New[2][2][0] = -1.0 / 210.0;
      Loc_M_New[2][2][1] = -17. / 420.;
      Loc_M_New[2][2][2] = 1.0 / 84.0;

       vector< vector< vector<type>>> Loc_G_Old =  vector< vector< vector<type>>>(3);

      for (int i = 0; i < 3; i++)
      {
         Loc_G_Old[i].resize(3);

      }
      // старая
      Loc_G_Old[0][0].resize(3);
      Loc_G_Old[0][1].resize(3);
      Loc_G_Old[0][2].resize(3);
      Loc_G_Old[0][0][0] = 7.0 / 3.0;
      Loc_G_Old[0][0][1] = -8. / 3.;
      Loc_G_Old[0][0][2] = 1.0 / 3.0;
      Loc_G_Old[0][1][0] = -8. / 3.;
      Loc_G_Old[0][1][1] = 16.0 / 3.0;
      Loc_G_Old[0][1][2] = -8. / 3.;
      Loc_G_Old[0][2][0] = 1. / 3.;
      Loc_G_Old[0][2][1] = -8. / 3.;
      Loc_G_Old[0][2][2] = 7. / 3.;
      //x
      Loc_G_Old[1][0].resize(3);
      Loc_G_Old[1][1].resize(3);
      Loc_G_Old[1][2].resize(3);
      Loc_G_Old[1][0][0] = 1. / 2.;
      Loc_G_Old[1][0][1] = -2.0 / 3.;
      Loc_G_Old[1][0][2] = 1.0 / 6.0;
      Loc_G_Old[1][1][0] = -2.0 / 3.0;
      Loc_G_Old[1][1][1] = 8.0 / 3.0;
      Loc_G_Old[1][1][2] = -2.;
      Loc_G_Old[1][2][0] = 1. / 6.;
      Loc_G_Old[1][2][1] = -2.;
      Loc_G_Old[1][2][2] = 11. / 6.;
      //1-x
      Loc_G_Old[2][0].resize(3);
      Loc_G_Old[2][1].resize(3);
      Loc_G_Old[2][2].resize(3);
      Loc_G_Old[2][0][0] = 11.0 / 6.0;
      Loc_G_Old[2][0][1] = -2.;
      Loc_G_Old[2][0][2] = 1.0 / 6.0;
      Loc_G_Old[2][1][0] = -2.;
      Loc_G_Old[2][1][1] = 8.0 / 3.0;
      Loc_G_Old[2][1][2] = -2. / 3.;
      Loc_G_Old[2][2][0] = 1. / 6.;
      Loc_G_Old[2][2][1] = -2. / 3.;
      Loc_G_Old[2][2][2] = 1. / 2.;

       vector< vector< vector<type>>> Loc_G_New =  vector< vector< vector<type>>>(3);

      for (int i = 0; i < 3; i++)
      {
         Loc_G_New[i].resize(3);
      }
      // Новая
      Loc_G_New[0][0].resize(3);
      Loc_G_New[0][1].resize(3);
      Loc_G_New[0][2].resize(3);
      Loc_G_New[0][0][0] = 1. / 2.;
      Loc_G_New[0][0][1] = -2.0 / 3.;
      Loc_G_New[0][0][2] = 1.0 / 6.0;
      Loc_G_New[0][1][0] = -2.0 / 3.0;
      Loc_G_New[0][1][1] = 8.0 / 3.0;
      Loc_G_New[0][1][2] = -2.;
      Loc_G_New[0][2][0] = 1. / 6.;
      Loc_G_New[0][2][1] = -2.;
      Loc_G_New[0][2][2] = 11. / 6.;
      //x
      Loc_G_New[1][0].resize(3);
      Loc_G_New[1][1].resize(3);
      Loc_G_New[1][2].resize(3);
      Loc_G_New[1][0][0] = 1.0 / 5.0;
      Loc_G_New[1][0][1] = -2./5.;
      Loc_G_New[1][0][2] = 1.0 / 5.0;
      Loc_G_New[1][1][0] = -2. / 5.;
      Loc_G_New[1][1][1] = 32.0 / 15.0;
      Loc_G_New[1][1][2] = -26. / 15.;
      Loc_G_New[1][2][0] = 1. / 5.;
      Loc_G_New[1][2][1] = -26. / 15.;
      Loc_G_New[1][2][2] = 23. / 15.;
      // 1-x
      Loc_G_New[2][0].resize(3);
      Loc_G_New[2][1].resize(3);
      Loc_G_New[2][2].resize(3);
      Loc_G_New[2][0][0] = 3.0 / 10.0;
      Loc_G_New[2][0][1] = -4. / 15.;
      Loc_G_New[2][0][2] = -1.0 / 30.0;
      Loc_G_New[2][1][0] = -4. / 15.;
      Loc_G_New[2][1][1] = 8.0 / 15.0;
      Loc_G_New[2][1][2] = -4. / 15.;
      Loc_G_New[2][2][0] = -1. / 30.;
      Loc_G_New[2][2][1] = -4. / 15.;
      Loc_G_New[2][2][2] = 3. / 10.;

      loc_m.resize(9);
      for (int i = 0; i < 9; i++)
      {
         loc_m[i].resize(9);
         for (int j = 0; j < 9; j++)
         {
            loc_m[i][j] = hr*hz*(hr*Loc_M_New[0][i % 3][j % 3] + rk*Loc_M_Old[0][i % 3][j % 3])* Loc_M_Old[0][i / 3][j / 3];
         }
      }
      loc_g.resize(4);
      for (int k = 0; k < 4; k++)
      {
         loc_g[k].resize(9);
      }
      for (int i = 0; i < 9; i++)
      {
         loc_g[0][i].resize(9);
         loc_g[1][i].resize(9);
         loc_g[2][i].resize(9);
         loc_g[3][i].resize(9);
         for (int j = 0; j < 9; j++)
         {
            loc_g[0][i][j] = hz * (hr * Loc_G_New[1][i % 3][j % 3] + rk * Loc_G_Old[1][i % 3][j % 3])
               * Loc_M_Old[1][i / 3][j / 3] / hr
               + hr * (hr * Loc_M_New[1][i % 3][j % 3] + rk * Loc_M_Old[1][i % 3][j % 3])
               * Loc_G_Old[1][i / 3][j / 3] / hz;

            loc_g[1][i][j] = hz * (hr * Loc_G_New[2][i % 3][j % 3] + rk * Loc_G_Old[2][i % 3][j % 3])
               * Loc_M_Old[1][i / 3][j / 3] / hr
               + hr * (hr * Loc_M_New[2][i % 3][j % 3] + rk * Loc_M_Old[2][i % 3][j % 3])
               * Loc_G_Old[1][i / 3][j / 3] / hz;

            loc_g[2][i][j] = hz * (hr * Loc_G_New[1][i % 3][j % 3] + rk * Loc_G_Old[1][i % 3][j % 3])
               * Loc_M_Old[2][i / 3][j / 3] / hr
               + hr * (hr * Loc_M_New[1][i % 3][j % 3] + rk * Loc_M_Old[1][i % 3][j % 3])
               * Loc_G_Old[2][i / 3][j / 3] / hz; 

            loc_g[3][i][j] = hz * (hr * Loc_G_New[2][i % 3][j % 3] + rk * Loc_G_Old[2][i % 3][j % 3])
               * Loc_M_Old[2][i / 3][j / 3] / hr
               + hr * (hr * Loc_M_New[2][i % 3][j % 3] + rk * Loc_M_Old[2][i % 3][j % 3])
               * Loc_G_Old[2][i / 3][j / 3] / hz;
         };
      }
   }
   // локальная матрица жесткости с разложеной лямбдой
   void loc_G(int num,  vector< vector<type>>& LocG)
   {
       vector<type> LocPhiLam =  vector<type>(4);
      loc_lam(num, LocPhiLam);
      LocG.resize(9);
      for (int i = 0; i < 9; i++)
      {
         LocG[i].resize(9);
         for (int j = 0; j < 9; j++)
         {
            LocG[i][j] = LocPhiLam[0] * loc_g[0][i][j]
               + LocPhiLam[1] * loc_g[1][i][j]
               + LocPhiLam[2] * loc_g[2][i][j]
               + LocPhiLam[3] * loc_g[3][i][j];
         }
      }
   }
   // считает лямбду в узлах элемента
   void loc_lam(int num,  vector<type>& LocPhiLam)
   {
      for (int i = 0; i < 4; i++)
      {
         LocPhiLam[i] = Lambda(grid[el_pints[1][num][i]][0], grid[el_pints[1][num][i]][1]);
      }
   }
   // считает локальный вектор правой части
   void loc_f( vector<type>& LocB, int num)
   {
      type hr, hz;
      hr = (grid[el_pints[1][num][1]][0] - grid[el_pints[1][num][0]][0]) / 2.;
      hz = (grid[el_pints[1][num][2]][1] - grid[el_pints[1][num][0]][1]) / 2.;
       vector<type> res =  vector<type>(9);
      for (int i = 0; i < 9; i++)
      {
         res[i] = F(grid[el_pints[1][num][0]][0] + hr * (i % 3), grid[el_pints[1][num][0]][1] + hz * (i / 3));
      }
      LocB = res;
   }
   // счет усредненной гаммы
   void loc_gamma(int num, type& gamma)
   {
      type hr, hz;
      type g = 0;
      hr = (grid[el_pints[1][num][1]][0] - grid[el_pints[1][num][0]][0]) / 2.;
      hz = (grid[el_pints[1][num][2]][1] - grid[el_pints[1][num][0]][1]) / 2.;
      for (int i = 0; i < 9; i++)
      {
         g += Gamma(grid[el_pints[1][num][0]][0] + hr * (i % 3), grid[el_pints[1][num][0]][1] + hz * (i / 3));
      }
      gamma = g / 9.;
   }
   // умножение вектора правой части на матрицу масс
   void m_mult_v( vector<type>& f)
   {
       vector<type> res =  vector<type>(9);
      for (int i = 0; i < 9; i++)
      {
         type sum = 0.0;
         for (int j = 0; j < 9; j++)
         {
            sum += f[i] * loc_m[i][j];
         }
         res[i] = sum;
      }
      f = res;
   }

   // первое краевое
   void one()
   {
      int metod = 2;
       fstream input;
      input.open("firstboundary.txt");
       vector<int> BoundaryBorder;
      int variable;
      while (input >> variable)
      {
         BoundaryBorder.push_back(variable);
      }
      input.close();
      for (int Border : BoundaryBorder)
      {
         double C = 1e+20;
         switch (Border)
         {
         case 0:
         {
            int num;
            for (int i = 0; i < (n_z - 1); i++)
            {
               num = i * (n_r - 1);
               type h = (grid[el_pints[1][num][2]][1] - grid[el_pints[1][num][0]][1]) / 2.0;
               type Ugi;
               for (int j = 0; j < 3; j++)
               {
                  Ugi = Ug(grid[el_pints[1][num][0]][0], grid[el_pints[1][num][0]][1] + h * j, Border);
                  if (metod == 2)
                  {
                     del_str(el_pints[0][num][3 * j], Ugi);
                     di[el_pints[0][num][3 * j]] = 1.;
                     B[el_pints[0][num][3 * j]] = Ugi;
                  }
                  if (metod == 1)
                  {
                     del_str_col(el_pints[0][num][3 * j], Ugi);
                     di[el_pints[0][num][3 * j]] = 1.;
                     B[el_pints[0][num][3 * j]] = Ugi;
                  }
                  if (metod == 0)
                  {
                     B[el_pints[0][num][3 * j]] = C * Ugi;
                     di[el_pints[0][num][3 * j]] = C * 1.;
                  }

               }
            }
            break;
         }
         case 1:
         {
            int num;
            for (int i = 0; i < (n_z - 1); i++)
            {
               num = i;
               type h = (grid[el_pints[1][num][1]][0] - grid[el_pints[1][num][0]][0]) / 2.0;
               type Ugi;
               for (int j = 0; j < 3; j++)
               {
                  Ugi = Ug(grid[el_pints[1][num][0]][0] + h * j, grid[el_pints[1][num][0]][1], Border);
                  if (metod == 2)
                  {
                     del_str(el_pints[0][num][j], Ugi);
                     B[el_pints[0][num][j]] = Ugi;
                     di[el_pints[0][num][j]] = 1.;
                  }if (metod == 1)
                  {
                     del_str_col(el_pints[0][num][j], Ugi);
                     B[el_pints[0][num][j]] = Ugi;
                     di[el_pints[0][num][j]] = 1.;
                  }
                  if (metod == 0)
                  {
                     B[el_pints[0][num][j]] = C * Ugi;
                     di[el_pints[0][num][j]] = C * 1.;
                  }
               }
            }
            break;
         }
         case 2:
         {
            int num;
            for (int i = 0; i < (n_z - 1); i++)
            {
               num = i * (n_r - 1) + n_r - 2;
               type h = (grid[el_pints[1][num][2]][1] - grid[el_pints[1][num][0]][1]) / 2.0;
               type Ugi;
               for (int j = 0; j < 3; j++)
               {
                  Ugi = Ug(grid[el_pints[1][num][1]][0], grid[el_pints[1][num][1]][1] + h * j, Border);
                  if (metod == 2)
                  {
                     del_str(el_pints[0][num][3 * j + 2], Ugi);
                     B[el_pints[0][num][3 * j + 2]] = Ugi;
                     di[el_pints[0][num][3 * j + 2]] = 1.;
                  } 
                  if (metod == 1)
                  {
                     del_str_col(el_pints[0][num][3 * j + 2], Ugi);
                     B[el_pints[0][num][3 * j + 2]] = Ugi;
                     di[el_pints[0][num][3 * j + 2]] = 1.;
                  }
                  if (metod == 0)
                  {
                     B[el_pints[0][num][3 * j + 2]] = C * Ugi;
                     di[el_pints[0][num][3 * j + 2]] = C * 1.;
                  }
               }
            }
            break;
         }
         case 3:
         {
            int num;
            for (int i = 0; i < (n_z - 1); i++)
            {
               num = i + (n_r - 1) * (n_z - 2);
               type h = (grid[el_pints[1][num][1]][0] - grid[el_pints[1][num][0]][1]) / 2.0;
               type Ugi;
               for (int j = 0; j < 3; j++)
               {
                  Ugi = Ug(grid[el_pints[1][num][2]][0] + h * j, grid[el_pints[1][num][2]][1], Border);
                  if (metod == 2)
                  {
                     del_str(el_pints[0][num][6 + j], Ugi);
                     B[el_pints[0][num][6 + j]] = Ugi;
                     di[el_pints[0][num][6 + j]] = 1.;
                  }          
                  if (metod == 1)
                  {
                     del_str_col(el_pints[0][num][6 + j], Ugi);
                     B[el_pints[0][num][6 + j]] = Ugi;
                     di[el_pints[0][num][6 + j]] = 1.;
                  }
                  if (metod == 0)
                  {
                     B[el_pints[0][num][6 + j]] = C * Ugi;
                     di[el_pints[0][num][6 + j]] = C * 1.;
                  }

               }
            }
            break;
         }
         }
      }
         if (metod == 0)
            au = al;
         if (metod == 1)
            au = al;
   }
   //зануляет нужные строки и столбцы и ставит 1 на диагональ
   void del_str_col(int Node, type Ugi)
   {
      int CurNode;
      for (int i = ig[Node]; i < ig[Node + 1]; i++)
      {
         B[jg[i]] -= Ugi * al[i];
         al[i] = 0.;
      }

      for (int i = Node; i < (2 * n_r - 1) * (2 * n_z - 1); i++)
      {
         CurNode = jg[ig[i]];
         int k = 0;
         while ((CurNode <= Node) && (k < (ig[i + 1] - ig[i])))
         {
            if (CurNode == Node)
            {
               B[i] -= al[ig[i] + k] * Ugi;
               al[ig[i] + k] = 0.;
            }
            k++;
            if (k < (ig[i + 1] - ig[i]))
               CurNode = jg[ig[i] + k];
         }
      }
   }
   //зануляет нужные строки и ставит 1 на диагональ
   void del_str(int Node, type Ugi)
   {
      int CurNode;
      for (int i = ig[Node]; i < ig[Node + 1]; i++)
      {
         al[i] = 0.;
      }

      for (int i = Node; i < (2 * n_r - 1) * (2 * n_z - 1); i++)
      {
         CurNode = jg[ig[i]];
         int k = 0;
         while ((CurNode <= Node) && (k < (ig[i + 1] - ig[i])))
         {
            if (CurNode == Node)
            {
               au[ig[i] + k] = 0.;
            }
            k++;
            if (k < (ig[i + 1] - ig[i]))
               CurNode = jg[ig[i] + k];
         }
      }
   }
   // второе краевое
   void two()
   {
       fstream input;
      input.open("secondboundary.txt");
       vector<int> BoundaryBorder;
      int variable;
      while (input >> variable)
      {
         BoundaryBorder.push_back(variable);
      }
      input.close();
      for (int Border : BoundaryBorder)
      {
         switch (Border)
         {
         case 0:
         {
             vector<type> LocTheta =  vector<type>(3);
            int num;
            for (int i = 0; i < (n_z - 1); i++)
            {
               num = i * (n_r - 1);
               type h = (grid[el_pints[1][num][2]][1] - grid[el_pints[1][num][0]][1]);
               loc_teta(num, h / 2.0, Border, LocTheta);
               B[el_pints[0][num][0]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
               B[el_pints[0][num][3]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
               B[el_pints[0][num][6]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
            }
            break;
         }
         case 1:
         {
             vector<type> LocTheta =  vector<type>(3);
            int num;
            for (int i = 0; i < (n_r - 1); i++)
            {
               num = i;
               type h = (grid[el_pints[1][num][1]][0] - grid[el_pints[1][num][0]][0]);
               loc_teta(num, h / 2., Border, LocTheta);
               B[el_pints[0][num][0]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
               B[el_pints[0][num][1]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
               B[el_pints[0][num][2]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
            }
            break;
         }
         case 2:
         {
             vector<type> LocTheta =  vector<type>(3);
            int num;
            for (int i = 0; i < (n_z - 1); i++)
            {
               num = i * (n_r - 1) + n_z - 2;
               type h = (grid[el_pints[1][num][2]][1] - grid[el_pints[1][num][0]][1]);
               loc_teta(num, h / 2., Border, LocTheta);
               B[el_pints[0][num][2]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
               B[el_pints[0][num][5]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
               B[el_pints[0][num][8]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
            }
            break;
         }
         case 3:
         {
             vector<type> LocTheta =  vector<type>(3);
            int num;
            for (int i = 0; i < (n_r - 1); i++)
            {
               num = i + (n_r - 1) * (n_z - 2);
               type h = (grid[el_pints[1][num][1]][0] - grid[el_pints[1][num][0]][0]);
               loc_teta(num, h / 2., Border, LocTheta);
               B[el_pints[0][num][6]] += h * (4. * LocTheta[0] + 2. * LocTheta[1] - LocTheta[2]) / 30.;
               B[el_pints[0][num][7]] += h * (2. * LocTheta[0] + 16. * LocTheta[1] + 2. * LocTheta[2]) / 30.;
               B[el_pints[0][num][8]] += h * (-1. * LocTheta[0] + 2. * LocTheta[1] + 4. * LocTheta[2]) / 30.;
            }
            break;
         }
         }
      }
   }
    // высчитвает значение тетты на границе
   void loc_teta(int num, type h, int Border,  vector<type>& LocTheta)
   {

      for (int i = 0; i < 3; i++)
      {
         switch (Border)
         {
         case 0:
         {
            LocTheta[i] = Theta(grid[el_pints[1][num][0]][0], grid[el_pints[1][num][0]][1] + h * i, Border);
            break;
         }
         case 1:
         {
            LocTheta[i] = Theta(grid[el_pints[1][num][0]][0] + h * i, grid[el_pints[1][num][0]][1], Border);
            break;
         }
         case 2:
         {
            LocTheta[i] = Theta(grid[el_pints[1][num][1]][0], grid[el_pints[1][num][1]][1] + h * i, Border);
            break;
         }
         case 3:
         {
            LocTheta[i] = Theta(grid[el_pints[1][num][2]][0] + h * i, grid[el_pints[1][num][2]][1], Border);
            break;
         }
         default:
            break;
         }
      }
   }
   // третье краевое
   void three()
   {
       fstream input;
      input.open("thirdboundary.txt");
       vector<int> BoundaryBorder;
      int variable;
      while (input >> variable)
      {
         BoundaryBorder.push_back(variable);
      }
      input.close();
       vector<type> LocUbeta =  vector<type>(3);
       vector< vector<type>> LocBondA =  vector< vector<type>>(3);
      for (int i = 0; i < 3; i++)
      {
         LocBondA[i].resize(3);
      }
      LocBondA[0][0] = 2. / 15.;
      LocBondA[0][1] = 1. / 15.;
      LocBondA[0][2] = -1. / 30.;
      LocBondA[1][0] = 1. / 15.;
      LocBondA[1][1] = 8. / 15.;
      LocBondA[1][2] = 1. / 15.;
      LocBondA[2][0] = -1. / 30.;
      LocBondA[2][1] = 1. / 15.;
      LocBondA[2][2] = 2. / 15.;
      for (int Border : BoundaryBorder)
      {
         switch (Border)
         {
         case 0:
         {
            int num;
            for (int i = 0; i < (n_z - 1); i++)
            {
               num = i * (n_r - 1);
               type beta = Beta();
               type h = (grid[el_pints[1][num][2]][1] - grid[el_pints[1][num][0]][1]);
               insert_loc_a(num, h, Border, LocBondA);
               loc_u_beta(num, h / 2., Border, LocUbeta);
               B[el_pints[0][num][0]] += h * beta * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
               B[el_pints[0][num][3]] += h * beta * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
               B[el_pints[0][num][6]] += h * beta * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
            }
            break;
         }
         case 1:
         {
            int num;
            for (int i = 0; i < (n_r - 1); i++)
            {
               num = i;
               type beta = Beta();
               type h = (grid[el_pints[1][num][1]][0] - grid[el_pints[1][num][0]][0]);
               insert_loc_a(num, h, Border, LocBondA);
               loc_u_beta(num, h / 2., Border, LocUbeta);
               B[el_pints[0][num][0]] += h * beta * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
               B[el_pints[0][num][1]] += h * beta * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
               B[el_pints[0][num][2]] += h * beta * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
            }
            break;
         }
         case 2:
         {
            int num;
            for (int i = 0; i < (n_z - 1); i++)
            {
               num = i * (n_r - 1) + n_r - 2;
               type beta = Beta();
               type h = (grid[el_pints[1][num][2]][1] - grid[el_pints[1][num][0]][1]);
               insert_loc_a(num, h, Border, LocBondA);
               loc_u_beta(num, h / 2., Border, LocUbeta);
               B[el_pints[0][num][2]] += h * beta * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
               B[el_pints[0][num][5]] += h * beta * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
               B[el_pints[0][num][8]] += h * beta * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
            }
            break;
         }
         case 3:
         {
            int num;
            for (int i = 0; i < (n_r - 1); i++)
            {
               num = i + (n_r - 1) * (n_z - 2);
               type beta = Beta();
               type h = (grid[el_pints[1][num][1]][0] - grid[el_pints[1][num][0]][0]);
               insert_loc_a(num, h, Border, LocBondA);
               loc_u_beta(num, h / 2., Border, LocUbeta);
               B[el_pints[0][num][6]] += h * beta * (4. * LocUbeta[0] + 2. * LocUbeta[1] - LocUbeta[2]) / 30.;
               B[el_pints[0][num][7]] += h * beta * (2. * LocUbeta[0] + 16. * LocUbeta[1] + 2. * LocUbeta[2]) / 30.;
               B[el_pints[0][num][8]] += h * beta * (-1. * LocUbeta[0] + 2. * LocUbeta[1] + 4. * LocUbeta[2]) / 30.;
            }
            break;
         }
         }
      }
   }
   // считает бетту на границу
   void loc_u_beta(int num, type h, int Border,  vector<type>& LocUbeta)
   {
      for (int i = 0; i < 3; i++)
      {
         switch (Border)
         {
         case 0:
         {
            LocUbeta[i] = Ubeta(grid[el_pints[1][num][0]][0], grid[el_pints[1][num][0]][1] + h * i, Border);
            break;
         }
         case 1:
         {
            LocUbeta[i] = Ubeta(grid[el_pints[1][num][0]][0] + h * i, grid[el_pints[1][num][0]][1], Border);
            break;
         }
         case 2:
         {
            LocUbeta[i] = Ubeta(grid[el_pints[1][num][1]][0], grid[el_pints[1][num][1]][1] + h * i, Border);
            break;
         }
         case 3:
         {
            LocUbeta[i] = Ubeta(grid[el_pints[1][num][2]][0] + h * i, grid[el_pints[1][num][2]][1], Border);
            break;
         }
         default:
            break;
         }
      }
   }
   void GetIndex(int i, int j, int& Index)
   {
      Index = ig[i];
      while (jg[Index] != j)
      {
         Index++;
      }
   }
   // вставляет локальную матрицу от краевого условия
   void insert_loc_a(int num, type h, int Border,  vector< vector<type>> LocBondA)
   {
      switch (Border)
      {
      case 0:
      {
         type beta = Beta();
         di[el_pints[0][num][0]] += h * beta * LocBondA[0][0];
         di[el_pints[0][num][3]] += h * beta * LocBondA[1][1];
         di[el_pints[0][num][6]] += h * beta * LocBondA[2][2];
         int Index;
         for (int i = 1; i < 3; i++)
         {
            for (int j = 0; j < i; j++)
            {
               GetIndex(el_pints[0][num][i * 3], el_pints[0][num][i * 3], Index);
               al[Index] += h * beta * LocBondA[i][j];
               au[Index] += h * beta * LocBondA[i][j];
            }
         }
         break;
      }
      case 1:
      {
         type beta = Beta();
         di[el_pints[0][num][0]] += h * beta * LocBondA[0][0];
         di[el_pints[0][num][1]] += h * beta * LocBondA[1][1];
         di[el_pints[0][num][2]] += h * beta * LocBondA[2][2];
         int Index;
         for (int i = 1; i < 3; i++)
         {
            for (int j = 0; j < i; j++)
            {
               GetIndex(el_pints[0][num][i], el_pints[0][num][j], Index);
               al[Index] += h * beta * LocBondA[i][j];
               au[Index] += h * beta * LocBondA[i][j];
            }
         }
         break;
      }
      case 2:
      {
         type beta = Beta();
         di[el_pints[0][num][2]] += h * beta * LocBondA[0][0];
         di[el_pints[0][num][5]] += h * beta * LocBondA[1][1];
         di[el_pints[0][num][8]] += h * beta * LocBondA[2][2];
         int Index;
         for (int i = 1; i < 3; i++)
         {
            for (int j = 0; j < i; j++)
            {
               GetIndex(el_pints[0][num][3 * i + 2], el_pints[0][num][3 * j + 2], Index);
               al[Index] += h * beta * LocBondA[i][j];
               au[Index] += h * beta * LocBondA[i][j];
            }
         }
         break;
      }
      case 3:
      {
         type beta = Beta();
         di[el_pints[0][num][6]] += h * beta * LocBondA[0][0];
         di[el_pints[0][num][7]] += h * beta * LocBondA[1][1];
         di[el_pints[0][num][8]] += h * beta * LocBondA[2][2];
         int Index;
         for (int i = 1; i < 3; i++)
         {
            for (int j = 0; j < i; j++)
            {
               GetIndex(el_pints[0][num][6 + i], el_pints[0][num][6 + j], Index);
               al[Index] += h * beta * LocBondA[i][j];
               au[Index] += h * beta * LocBondA[i][j];
            }
         }
         break;
      }
      }
   }
   // решение слау Los
   void SLAU()
   {
      Testing("Do reshenia");
      LUsq(di, al, au, ig, jg,q, B, (2 * n_r - 1) * (2 * n_z - 1));
      Testing("Posle reshenia");

   };
   void Testing(string str)
   {
      ofstream Test;
      Test.open("Test.txt", std::ios_base::app);

      Test << "\n"+str+"\nvar Di = new[] {";
      for (size_t i = 0; i < di.size(); i++)
      {
         Test << di[i] << ", ";
      }
      Test << "};\nvar Al = new[] {";
      for (size_t i = 0; i < al.size(); i++)
      {
         Test << al[i] << ", ";
      }
      Test << "};\nvar Au = new[] {";
      for (size_t i = 0; i < au.size(); i++)
      {
         Test << au[i] << ", ";
      }
      Test << "};\nvar Ig = new[] {";
      for (size_t i = 0; i < ig.size(); i++)
      {
         Test << ig[i] << ", ";
      }
      Test << "};\nvar Jg = new[] {";
      for (size_t i = 0; i < jg.size(); i++)
      {
         Test << jg[i] << ", ";
      }
      Test << "};\nvar B = new[] {";
      for (size_t i = 0; i < B.size(); i++)
      {
         Test << B[i] << ", ";
      }
      Test << "};\nvar q = new[] {";
      for (size_t i = 0; i < q.size(); i++)
      {
         Test << q[i] << ", ";
      }
      Test << "};\n";
      Test.close();
   }
   // вывод
   void answer()
   {
      std::cout << std::fixed << std::setprecision(15);
      type hr = (grid[el_pints[1][0][1]][0] - grid[el_pints[1][0][0]][0])/2.;
      type hz = (grid[el_pints[1][0][2]][1] - grid[el_pints[1][0][0]][1])/2.;
      for (int i = 0; i < (2 * n_r - 1) * (2 * n_z - 1); i++)
      {
         cout << setw(17) << q[i]<< setw(14) << "  " << RealF(grid[el_pints[1][0][0]][0] + hr * (i % (2 * n_r - 1)), grid[el_pints[1][0][0]][1] + hz * (i / (2 * n_z - 1)))
             << "  " << setw(17) << q[i] - RealF(grid[el_pints[1][0][0]][0] + hr * (i % (2 * n_r - 1)), grid[el_pints[1][0][0]][1] + hz * (i / (2 * n_z - 1))) << endl;
      }
   }
};

int main()
{
   MKA MX;
   MX.makingGrid();
   MX.profile();
   MX.global_build();
   MX.Testing("Do kraev");
   MX.two();
   MX.three();
   MX.one();
   MX.SLAU();
   MX.answer();
}
