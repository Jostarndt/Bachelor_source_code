#include <amandus/apps.h>
#include <amandus/jost/explicit.h>
#include <amandus/jost/implicit.h>
#include <amandus/jost/matrix.h>
#include <deal.II/algorithms/newton.h>
#include <deal.II/algorithms/theta_timestepping.h>
#include <deal.II/base/function.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/dof_output_operator.h>
#include <deal.II/numerics/dof_output_operator.templates.h>

#include <boost/scoped_ptr.hpp>

#include <math.h>
#include <iostream>
#include <random>

template <int dim>
class Startup : public dealii::Function<dim>
{
public:
  Startup();
  virtual void vector_value_list(const std::vector<Point<dim>>& points,
                                 std::vector<Vector<double>>& values) const;
};

template <int dim>
Startup<dim>::Startup()
  : Function<dim>(2)
{
}

template <int dim>
void
Startup<dim>::vector_value_list(const std::vector<Point<dim>>& points,
                                std::vector<Vector<double>>& values) const
{
  AssertDimension(points.size(), values.size());


	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(-1, 1);
	for (unsigned int k = 0; k< points.size(); ++k)
	{
	values[k](0) = 10+dis(gen);
	values[k](1) = 9+dis(gen);
	}
}	

int
main(int argc, const char** argv)
{
  const unsigned int d = 2;

  std::ofstream logfile("deallog");
  deallog.attach(logfile);
  deallog.depth_console(10);

  AmandusParameters param;
  param.read(argc, argv);
  param.log_parameters(deallog);

  param.enter_subsection("Discretization");
  boost::scoped_ptr<const FiniteElement<d>> fe(FETools::get_fe_by_name<d, d>(param.get("FE")));

  Triangulation<d> tr(Triangulation<d>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr, 0., 1.);
  tr.refine_global(param.get_integer("Refinement"));
  param.leave_subsection();

  jost::Parameters parameters;
  parameters.alpha0 = 1.;
  parameters.alpha1 = 10;
	parameters.alpha =1.5; 
	parameters.a = 92;
	parameters.b = 64;
	parameters.K= 0.1;
	parameters.rho= 18.5;
	parameters.gamma = 100;
//  parameters.A = 3.4;
//  parameters.B = 1.;
  jost::Matrix<d> matrix_integrator(parameters);
  jost::ExplicitResidual<d> explicit_integrator(parameters);
  explicit_integrator.input_vector_names.push_back("Previous iterate");
  jost::ImplicitResidual<d> implicit_integrator(parameters);
  implicit_integrator.input_vector_names.push_back("Newton iterate");


  AmandusApplication<d> app(tr, *fe);
  AmandusResidual<d> expl(app, explicit_integrator);
  AmandusSolve<d> solver(app, matrix_integrator);
  AmandusResidual<d> residual(app, implicit_integrator);

  // Set up timestepping algorithm with embedded Newton solver

  param.enter_subsection("Output");
  Algorithms::DoFOutputOperator<Vector<double>, d> newout;
  newout.parse_parameters(param);
  newout.initialize(app.dofs());
  param.leave_subsection();

  Algorithms::Newton<Vector<double>> newton(residual, solver);
  newton.parse_parameters(param);

  Algorithms::ThetaTimestepping<Vector<double>> timestepping(expl, newton);
  timestepping.set_output(newout);
  timestepping.parse_parameters(param);

  // Now we prepare for the actual timestepping

  timestepping.notify(dealii::Algorithms::Events::remesh);
  dealii::Vector<double> solution;
  app.setup_system();
  app.setup_vector(solution);

  Startup<d> startup;
  VectorTools::interpolate(app.dofs(), startup, solution);

  dealii::AnyData indata;
  indata.add(&solution, "solution");
  dealii::AnyData outdata;
  timestepping(indata, outdata);
}
