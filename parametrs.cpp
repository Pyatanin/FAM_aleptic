#include "parmetrs.h"

double RealF(double x, double y)
{
   return x;
}

double Lambda(double x, double y)
{
   return 1.;
}

double Gamma(double x, double y)
{
   return 1.;
}

double Theta(double x, double y, int Border)
{
   switch (Border)
   {
   case 0:
   {
      return -1;
      break;
   }
   case 1:
   {
      return 0;
      break;
   }
   case 2:
   {
      return 1;
      break;
   }
   case 3:
   {
      return 0;
      break;
   }
   }
}

double Ubeta(double x, double y, int Border)
{
   switch (Border)
   {
   case 0:
   {
      return 0;
      break;
   }
   case 1:
   {
      return x;
      break;
   }
   case 2:
   {
      return 3;
      break;
   }
   case 3:
   {
      return x;
      break;
   }
   }
}

double Ug(double x, double y, int Border)
{
   switch (Border)
   {
   case 0:
   {
      return 1.;
      break;
   }
   case 1:
   {
      return x;
      break;
   }
   case 2:
   {
      return 2.;
      break;
   }
   case 3:
   {
      return x;
      break;
   }
   }
}

double Beta()
{
   return 1.0;
}

double F(double x, double y)
{
   return -1/x+x;
}