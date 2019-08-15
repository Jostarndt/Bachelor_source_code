#ifndef __jost_parameters_h
#define __jost_parameters_h

#include <deal.II/base/subscriptor.h>

namespace jost
{
struct Parameters : public dealii::Subscriptor
{
  double alpha0;
  double alpha1;
	double alpha;
	double K;
	double rho;
	double a;
	double b;
	double gamma;
	
	double f;
	double k;
};
}

#endif
