/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
#define ACADOS_OCP_NLP_OCP_NLP_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/types.h"
#include "acados/utils/external_function_generic.h"



/************************************************
* structures
************************************************/

/************************************************
* dims
************************************************/

/* ocp_nlp */

typedef struct
{
	void *cost_dims;
    int *nx;
    int *nu;
    int *nb;  // nbx + nbu
    int *nbx;
    int *nbu;
    int *ng;  // number of general linear constraints
    int *nh;  // number of path constraints - ONLY difference with ocp_qp_dims atm
    int *ns;  // number of soft constraints
    int *num_stages;
    int N;
	int memsize;
} ocp_nlp_dims;

/* ocp_nlp_cost_ls */

// TODO do ocp_nlp_cost instread ???

typedef struct
{
	int *nv; // number of variables
	int *ny; // number of outputs
    int N;
	int memsize;
} ocp_nlp_cost_ls_dims;



/************************************************
* cost
************************************************/

/* least squares */

typedef struct
{
	ocp_nlp_cost_ls_dims *dims;
	external_function_generic **nls_jac; // array of N+1 pointers; evaluation and jacobian of ls residuals
	struct blasfeo_dmat *Cyt;
	struct blasfeo_dmat *W;
    struct blasfeo_dvec *y_ref;
	int *nls_mask; // nonlinear least squares mask
} ocp_nlp_cost_ls;



/************************************************
* dynamics
************************************************/

/* ERK */

typedef struct
{
	ocp_nlp_dims *dims; // TODO dynamics dims ???
	external_function_generic **forw_vde;
	external_function_generic **adj_vde;
	external_function_generic **jac_ode;
} ocp_nlp_dynamics_erk;

/* IRK */

typedef struct
{
	ocp_nlp_dims *dims; // TODO dynamics dims ???
	external_function_generic **ode;
	external_function_generic **jac_x;
	external_function_generic **jac_xdot;
	external_function_generic **jac_u;
} ocp_nlp_dynamics_irk;

/* lifted IRK */

typedef struct
{
	ocp_nlp_dims *dims; // TODO dynamics dims ???
	external_function_generic **forw_vde;
	external_function_generic **jac_ode;
} ocp_nlp_dynamics_lifted_irk;



/************************************************
* constraints
************************************************/

typedef struct
{
	ocp_nlp_dims *dims; // TODO constraints dims ???
    int **idxb;
	struct blasfeo_dvec *d;
	struct blasfeo_dmat *DCt;
}
ocp_nlp_constraints;



/************************************************
* in
************************************************/

typedef struct
{
    ocp_nlp_dims *dims;

    // double **lh;
    // double **uh;
    // ocp_nlp_function *h;  // nonlinear path constraints

	// TODO array of structures or structures of arrays ???
	// void **cost; // ???
    void *cost;
	void *dynamics;
	void *constraints;

    // TODO(rien): what about invariants, e.g., algebraic constraints?

    bool freezeSens;  // TODO(dimitris): shouldn't this be in the integrator args?
} ocp_nlp_in;



/************************************************
* out
************************************************/

typedef struct
{
    ocp_nlp_dims *dims;
	struct blasfeo_dvec *ux;
	struct blasfeo_dvec *pi;
	struct blasfeo_dvec *lam;
	struct blasfeo_dvec *t;
	int memsize;
} ocp_nlp_out;



/************************************************
* memory
************************************************/

typedef struct
{
    ocp_nlp_dims *dims;
	struct blasfeo_dvec *cost_grad;
	struct blasfeo_dvec *dyn_fun;
	struct blasfeo_dvec *dyn_adj;
	struct blasfeo_dvec *ineq_fun;
	struct blasfeo_dvec *ineq_adj;
	int memsize;
} ocp_nlp_mem;



/************************************************
* residuals
************************************************/

typedef struct
{
    ocp_nlp_dims *dims;
	struct blasfeo_dvec *res_g; // stationarity
	struct blasfeo_dvec *res_b; // dynamics
	struct blasfeo_dvec *res_d; // inequality constraints
	struct blasfeo_dvec *res_m; // complementarity
	double inf_norm_res_g;
	double inf_norm_res_b;
	double inf_norm_res_d;
	double inf_norm_res_m;
	int memsize;
} ocp_nlp_res;



/************************************************
* function pointers
************************************************/

typedef struct
{
	// cost
	int (*cost_calculate_size) (ocp_nlp_cost_ls_dims *); // TODO ocp_nlp_cost_dims
	void *(*cost_assign) (ocp_nlp_cost_ls_dims *, void *); // TODO ocp_nlp_dims

	// dynamics
	int (*dynamics_calculate_size) (ocp_nlp_dims *); // calculate size
	void *(*dynamics_assign) (ocp_nlp_dims *, void *); // assign
	void (*dynamics_to_sim_in) (void *, sim_in **);

	// constraints
	int (*constraints_calculate_size) (ocp_nlp_dims *); // calculate size
	void *(*constraints_assign) (ocp_nlp_dims *, void *); // assign

	// all the others
    int (*fun)(ocp_nlp_in *qp_in, ocp_nlp_out *qp_out, void *args, void *mem, void *work);
    int (*calculate_args_size)(ocp_nlp_dims *dims, void *solver_);
    void *(*assign_args)(ocp_nlp_dims *dims, void *solver_, void *raw_memory);
    void *(*copy_args)(ocp_nlp_dims *dims, void *solver_, void *raw_memory, void *source_);
    void (*initialize_default_args)(void *args);
    int (*calculate_memory_size)(ocp_nlp_dims *dims, void *args);
    void *(*assign_memory)(ocp_nlp_dims *dims, void *args, void *raw_memory);
    int (*calculate_workspace_size)(ocp_nlp_dims *dims, void *args);
    ocp_qp_xcond_solver_config *qp_solver;
    sim_solver_config **sim_solvers;
} ocp_nlp_solver_config;





/************************************************
* headers
************************************************/

/************************************************
* dims
************************************************/

/* ocp_nlp */

//
int ocp_nlp_dims_calculate_size(int N);
//
ocp_nlp_dims *ocp_nlp_dims_assign(int N, void *raw_memory);
//
void ocp_nlp_dims_init(int *nx, int *nu, int *nbx, int *nbu, int *ng, int *nh, int *ns, void *ocp_nlp_cost_dims, ocp_nlp_dims *dims);

/* ocp_nlp_cost_ls */

//
int ocp_nlp_cost_ls_dims_calculate_size(int N);
//
ocp_nlp_cost_ls_dims *ocp_nlp_cost_ls_dims_assign(int N, void *raw_memory);
//
void ocp_nlp_cost_ls_dims_init(int *nv, int *ny, ocp_nlp_cost_ls_dims *dims);

/************************************************
* cost
************************************************/

/* least squares */

//
int ocp_nlp_cost_ls_calculate_size(ocp_nlp_cost_ls_dims *dims);
//
ocp_nlp_cost_ls *ocp_nlp_cost_ls_assign(ocp_nlp_cost_ls_dims *dims, void *raw_memory);

/************************************************
* dynamics
************************************************/

/* ERK */

//
int ocp_nlp_dynamics_erk_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_dynamics_erk *ocp_nlp_dynamics_erk_assign(ocp_nlp_dims *dims, void *raw_memory);
//
void ocp_nlp_dynamics_erk_to_sim_in(ocp_nlp_dynamics_erk *dynamics, sim_in **sim);

/* IRK */

//
int ocp_nlp_dynamics_irk_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_dynamics_irk *ocp_nlp_dynamics_irk_assign(ocp_nlp_dims *dims, void *raw_memory);
//
void ocp_nlp_dynamics_irk_to_sim_in(ocp_nlp_dynamics_irk *dynamics, sim_in **sim);

/* lifted IRK */

//
int ocp_nlp_dynamics_lifted_irk_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_dynamics_lifted_irk *ocp_nlp_dynamics_lifted_irk_assign(ocp_nlp_dims *dims, void *raw_memory);
//
void ocp_nlp_dynamics_lifted_irk_to_sim_in(ocp_nlp_dynamics_lifted_irk *dynamics, sim_in **sim);

/************************************************
* constraints
************************************************/

//
int ocp_nlp_constraints_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_constraints *ocp_nlp_constraints_assign(ocp_nlp_dims *dims, void *raw_memory);

/************************************************
* in
************************************************/

//
int ocp_nlp_in_calculate_size(ocp_nlp_dims *dims, ocp_nlp_solver_config *config);
//
ocp_nlp_in *assign_ocp_nlp_in(ocp_nlp_dims *dims, int num_stages, void *raw_memory, ocp_nlp_solver_config *config);

/************************************************
* out
************************************************/

//
int ocp_nlp_out_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_out *assign_ocp_nlp_out(ocp_nlp_dims *dims, void *raw_memory);

/************************************************
* memory
************************************************/

//
int ocp_nlp_mem_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_mem *ocp_nlp_mem_assign(ocp_nlp_dims *dims, void *raw_memory);

/************************************************
* residuals
************************************************/

//
int ocp_nlp_res_calculate_size(ocp_nlp_dims *dims);
//
ocp_nlp_res *ocp_nlp_res_assign(ocp_nlp_dims *dims, void *raw_memory);
//
void ocp_nlp_res_compute(ocp_nlp_in *in, ocp_nlp_out *out, ocp_nlp_res *res, ocp_nlp_mem *mem); //ocp_nlp_res_workspace *work);









/************************************************
* ?????
************************************************/

int number_of_primal_vars(ocp_nlp_dims *dims);

void cast_nlp_dims_to_qp_dims(ocp_qp_dims *qp_dims, ocp_nlp_dims *nlp_dims);

void cast_nlp_dims_to_sim_dims(sim_dims *sim_dims, ocp_nlp_dims *nlp_dims, int stage);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_OCP_NLP_OCP_NLP_COMMON_H_
