
#ifndef __jost_matrix_h
#define __jost_matrix_h

#include <amandus/jost/parameters.h>
#include <amandus/integrator.h>
#include <deal.II/integrators/divergence.h>
#include <deal.II/integrators/l2.h>
#include <deal.II/integrators/laplace.h>
#include <deal.II/meshworker/integration_info.h>

using namespace dealii;
using namespace LocalIntegrators;


namespace jost
{

template <int dim>
class Matrix : public AmandusIntegrator<dim>
{
public:
  Matrix(const Parameters& par);

  virtual void cell(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const;
  virtual void boundary(MeshWorker::DoFInfo<dim>& dinfo,
                        MeshWorker::IntegrationInfo<dim>& info) const;
  virtual void face(MeshWorker::DoFInfo<dim>& dinfo1, MeshWorker::DoFInfo<dim>& dinfo2,
                    MeshWorker::IntegrationInfo<dim>& info1,
                    MeshWorker::IntegrationInfo<dim>& info2) const;

private:
  SmartPointer<const Parameters, class Matrix<dim>> parameters;
};

template <int dim>
Matrix<dim>::Matrix(const Parameters& par)
  : parameters(&par)
{
  this->use_boundary = false;
  this->use_face = true;
  this->input_vector_names.push_back("Newton iterate");
}

template <int dim>
void
Matrix<dim>::cell(MeshWorker::DoFInfo<dim>& dinfo, MeshWorker::IntegrationInfo<dim>& info) const
{
  AssertDimension(dinfo.n_matrices(), 4);
  // Assert (info.values.size() >0, ExcInternalError());
  const double factor = (this->timestep == 0.) ? 1. : this->timestep;
  if (this->timestep != 0.)
  {
    L2::mass_matrix(dinfo.matrix(0, false).matrix, info.fe_values(0));
    L2::mass_matrix(dinfo.matrix(3, false).matrix, info.fe_values(0));
  }

  Laplace::cell_matrix(
    dinfo.matrix(0, false).matrix, info.fe_values(0), parameters->alpha0 * factor);
  Laplace::cell_matrix(
    dinfo.matrix(3, false).matrix, info.fe_values(0), parameters->alpha1 * factor);
  if (info.values.size() > 0)
  {
    AssertDimension(info.values[0].size(), 2);
    AssertDimension(info.values[0][0].size(), info.fe_values(0).n_quadrature_points);
    AssertDimension(info.values[0][1].size(), info.fe_values(0).n_quadrature_points);
    std::vector<double> Du_ru(info.fe_values(0).n_quadrature_points);
    std::vector<double> Dv_ru(info.fe_values(0).n_quadrature_points);
    std::vector<double> Dv_rv(info.fe_values(0).n_quadrature_points);
    std::vector<double> Du_rv(info.fe_values(0).n_quadrature_points);
    for (unsigned int k = 0; k < Du_ru.size(); ++k)
    {
      const double u = info.values[0][0][k];
      const double v = info.values[0][1][k];

      Du_ru[k] =factor* (-v*v -u);
      Dv_ru[k] =factor* (-2*u*v);
      Dv_rv[k] =factor* (2*u*v -(parameters->f+parameters->k));
      Du_rv[k] =factor* (v*v);
    }
    L2::weighted_mass_matrix(dinfo.matrix(0, false).matrix, info.fe_values(0), Du_ru);
    L2::weighted_mass_matrix(dinfo.matrix(1, false).matrix, info.fe_values(0), Dv_ru);
    L2::weighted_mass_matrix(dinfo.matrix(3, false).matrix, info.fe_values(0), Dv_rv);
    L2::weighted_mass_matrix(dinfo.matrix(2, false).matrix, info.fe_values(0), Du_rv);
  }
}

template <int dim>
void
Matrix<dim>::boundary(MeshWorker::DoFInfo<dim>&, typename MeshWorker::IntegrationInfo<dim>&) const
{
}

template <int dim>
void
Matrix<dim>::face(MeshWorker::DoFInfo<dim>& dinfo1, MeshWorker::DoFInfo<dim>& dinfo2,
                  MeshWorker::IntegrationInfo<dim>& info1,
                  MeshWorker::IntegrationInfo<dim>& info2) const
{
  const unsigned int deg = info1.fe_values(0).get_fe().tensor_degree();
  const double factor = (this->timestep == 0.) ? 1. : this->timestep;
  Laplace::ip_matrix(dinfo1.matrix(0, false).matrix,
                     dinfo1.matrix(0, true).matrix,
                     dinfo2.matrix(0, true).matrix,
                     dinfo2.matrix(0, false).matrix,
                     info1.fe_values(0),
                     info2.fe_values(0),
                     2.*Laplace::compute_penalty(dinfo1, dinfo2, deg, deg),
                     parameters->alpha0 * factor);
  Laplace::ip_matrix(dinfo1.matrix(3, false).matrix,
                     dinfo1.matrix(3, true).matrix,
                     dinfo2.matrix(3, true).matrix,
                     dinfo2.matrix(3, false).matrix,
                     info1.fe_values(0),
                     info2.fe_values(0),
                     2.*Laplace::compute_penalty(dinfo1, dinfo2, deg, deg),
                     parameters->alpha1 * factor);
}
}

#endif
