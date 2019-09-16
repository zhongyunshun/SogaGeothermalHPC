#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <iostream>
#include <fstream>

#include <map>

using namespace dealii;

namespace Geothermal
{
using namespace dealii;

template <int dim>
class HeatEquation<dim>
{
public:
    HeatEquation();
    void run();

private:
    void gridinput();
    void setupsys();
    void solve_time_step();
    void output_results() const;
    Triangulation<dim> triangulation; //grid
    FE_Q<dim> fe;                     //element
    DoFHandler<dim> dof_handler;      //grid<->eleemnt

    ConstraintMatrix constraints; // hanging node

    SparsityPattern sparsity_pattern;    // sparsity
    SparseMatrix<double> mass_matrix;    // M
    SparseMatrix<double> laplace_matrix; //A
    SparseMatrix<double> system_matrix;  //M + k*theta*A

    Vector<double> solution;     // solution at n
    Vector<double> old_solution; //solution at n-1
    Vector<double> system_rhs;   //rhs

    double time;
    double time_step;
    unsigned int timestep_number;

    const double theta;
}

template <int dim>
class Righthandside : public Function<dim>
{
public:
    Righthandside() : Function<dim>(), period(0.5)
    {
    }
    virtual double value(const Point<dim> &p, const unsigned int component = 0) const;

private:
    const double period;
}

template <int dim>
class BoundaryValue : public Function<dim>
{
public:
}
} // namespace Geothermal

main()
{
    HeatEquation<3> heat_eq;
    heat_eq.run;
}