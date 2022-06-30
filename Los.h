#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>

using namespace std;

typedef vector<double> Vector;
void LUsq(vector<double> di, vector<double> a_l, vector<double> a_u, vector<int> ig, vector<int> jg, vector<double>& q, vector<double> B, int N);
void LoS_precond(vector<double>& dif, vector<double>& a_l_f, vector<double>& a_u_f, vector<double>& di, vector<double>& a_l, vector<double>& a_u, vector<int>& ig, vector<int>& jg, vector<double>& q, vector<double> B, int N);
vector<double> Mult(vector<double>& v, vector<double> di, vector<double> a_l, vector<double> a_u, vector<int> ig, vector<int> jg);
vector<double> LUDirect(const  vector<double>& b, vector<double>& _dif, vector<double>& _alf, vector<int>& ig, vector<int>& jg);
vector<double> LUReverse(const  vector<double>& b, vector<double>& _dif, vector<double>& _alf, vector<int>& ig, vector<int>& jg, int N);
double scalar_prod(vector<double>& x, vector<double>& y);
void print_(vector<double>& di, vector<double>& ggu, vector<double>& ggl, vector<int>& ig, vector<int>& jg, vector<double>& B, int N);