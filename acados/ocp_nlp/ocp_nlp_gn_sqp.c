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

#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"

// external
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
// blasfeo
#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_collocation_utils.h" // TODO remove ???
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/mem.h"



static int get_max_sim_workspace_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, ocp_nlp_gn_sqp_opts *args)
{
	/* ocp_qp_xcond_solver_config *qp_solver = config->qp_solver; */
	ocp_nlp_dynamics_config **dynamics = config->dynamics;

    int sim_work_size;

    int max_sim_work_size = 0;

    for (int ii = 0; ii < dims->N; ii++)
    {
        // sim_in_size = sim_in_calculate_size(dims->sim[ii]);
        // if (sim_in_size > *max_sim_in_size) *max_sim_in_size = sim_in_size;
        // sim_out_size = sim_out_calculate_size(dims->sim[ii]);
        // if (sim_out_size > *max_sim_out_size) *max_sim_out_size = sim_out_size;
        sim_work_size = dynamics[ii]->sim_solver->workspace_calculate_size(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim, args->sim_solvers_opts[ii]);
        if (sim_work_size > max_sim_work_size) max_sim_work_size = sim_work_size;
    }
    return max_sim_work_size;
}



/************************************************
* arguments
************************************************/



int ocp_nlp_gn_sqp_opts_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims)
{
	ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
//	sim_solver_config **sim_solvers = config->sim_solvers;
	ocp_nlp_dynamics_config **dynamics = config->dynamics;

    int size = 0;

    size += sizeof(ocp_nlp_gn_sqp_opts);

    size += qp_solver->opts_calculate_size(qp_solver, dims->qp_solver);

    size += dims->N*sizeof(void *);  //sim_solvers_opts

    for (int ii = 0; ii < dims->N; ii++)
    {
        size += dynamics[ii]->sim_solver->opts_calculate_size(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim);
    }

    return size;
}



ocp_nlp_gn_sqp_opts *ocp_nlp_gn_sqp_opts_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, void *raw_memory)
{
	ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
//	sim_solver_config **sim_solvers = config->sim_solvers;
	ocp_nlp_dynamics_config **dynamics = config->dynamics;

    ocp_nlp_gn_sqp_opts *opts;

    char *c_ptr = (char *) raw_memory;

    opts = (ocp_nlp_gn_sqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_gn_sqp_opts);

    opts->qp_solver_opts = qp_solver->opts_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += qp_solver->opts_calculate_size(qp_solver, dims->qp_solver);

    opts->sim_solvers_opts = (void **) c_ptr;
    c_ptr += dims->N*sizeof(void *);

    for (int ii = 0; ii < dims->N; ii++)
    {
        opts->sim_solvers_opts[ii] = dynamics[ii]->sim_solver->opts_assign(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim, c_ptr);
        c_ptr += dynamics[ii]->sim_solver->opts_calculate_size(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim);
    }

    assert((char*)raw_memory + ocp_nlp_gn_sqp_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_gn_sqp_opts_initialize_default(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, ocp_nlp_gn_sqp_opts *opts)
{

	ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
	ocp_nlp_dynamics_config **dynamics = config->dynamics;

	int ii;

	int N = dims->N;

	opts-> maxIter = 20;
	opts->min_res_g = 1e-12;
	opts->min_res_b = 1e-12;
	opts->min_res_d = 1e-12;
	opts->min_res_m = 1e-12;

	qp_solver->opts_initialize_default(qp_solver, opts->qp_solver_opts);

	for (ii=0; ii<N; ii++)
	{
		dynamics[ii]->sim_solver->opts_initialize_default(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim, opts->sim_solvers_opts[ii]);
	}

	return;

}



/************************************************
* memory
************************************************/



int ocp_nlp_gn_sqp_memory_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, ocp_nlp_gn_sqp_opts *args)
{
	ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
//	sim_solver_config **sim_solvers = config->sim_solvers;
	ocp_nlp_dynamics_config **dynamics = config->dynamics;

	// loop index
	int ii;

	// extract dims
    int N = dims->N;
	ocp_nlp_cost_dims **cost_dims = dims->cost;
	int ny;

    int size = 0;

    size += sizeof(ocp_nlp_gn_sqp_memory);

    size += qp_solver->memory_calculate_size(qp_solver, dims->qp_solver, args->qp_solver_opts);

    size += N*sizeof(void *);  // sim_solvers_mem

    for (int ii = 0; ii < N; ii++)
    {
        size += dynamics[ii]->sim_solver->memory_calculate_size(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim, args->sim_solvers_opts[ii]);
    }

	// nlp res
	size += ocp_nlp_res_calculate_size(dims);

	// nlp mem
	size += ocp_nlp_memory_calculate_size(config, dims);

	// blasfeo stuff
    size += 1*(N+1)*sizeof(struct blasfeo_dmat); // W_chol
    size += 1*(N+1)*sizeof(struct blasfeo_dvec); // ls_res

    for (ii = 0; ii < N+1; ii++)
    {
		ny = cost_dims[ii]->ny;
        size += 1*blasfeo_memsize_dvec(ny); // ls_res
        size += 1*blasfeo_memsize_dmat(ny, ny); // W_chol
    }

    size += 8; // initial align
    size += 8; // blasfeo_struct align
    size += 64; // blasfeo_mem align

//    make_int_multiple_of(64, &size);

    return size;
}



ocp_nlp_gn_sqp_memory *ocp_nlp_gn_sqp_memory_assign(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, ocp_nlp_gn_sqp_opts *args, void *raw_memory)
{
	ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
//	sim_solver_config **sim_solvers = config->sim_solvers;
	ocp_nlp_dynamics_config **dynamics = config->dynamics;

    char *c_ptr = (char *) raw_memory;

	// loop index
	int ii;

	// extract dims
    int N = dims->N;
	ocp_nlp_cost_dims **cost_dims = dims->cost;
	int ny;

	// initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)c_ptr;
    c_ptr += sizeof(ocp_nlp_gn_sqp_memory);

    // QP solver
    mem->qp_solver_mem = qp_solver->memory_assign(qp_solver, dims->qp_solver, args->qp_solver_opts, c_ptr);
    c_ptr += qp_solver->memory_calculate_size(qp_solver, dims->qp_solver, args->qp_solver_opts);

    mem->sim_solvers_mem = (void **) c_ptr;
    c_ptr += N*sizeof(void *);

    for (ii = 0; ii < N; ii++)
    {
//        mem->sim_solvers_mem[ii] = sim_solvers[ii]->memory_assign(sim_solvers[ii], dims->sim[ii], args->sim_solvers_opts[ii], c_ptr);
//        c_ptr += sim_solvers[ii]->memory_calculate_size(sim_solvers[ii], dims->sim[ii], args->sim_solvers_opts[ii]);
        mem->sim_solvers_mem[ii] = dynamics[ii]->sim_solver->memory_assign(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim, args->sim_solvers_opts[ii], c_ptr);
        c_ptr += dynamics[ii]->sim_solver->memory_calculate_size(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim, args->sim_solvers_opts[ii]);
    }

	// nlp res
	mem->nlp_res = ocp_nlp_res_assign(dims, c_ptr);
	c_ptr += mem->nlp_res->memsize;

	// nlp mem
	mem->nlp_mem = ocp_nlp_memory_assign(config, dims, c_ptr);
	c_ptr += ocp_nlp_memory_calculate_size(config, dims);

	// blasfeo_struct align
    align_char_to(8, &c_ptr);

    // set up local SQP data
    assign_blasfeo_dmat_structs(N+1, &mem->W_chol, &c_ptr);
    assign_blasfeo_dvec_structs(N+1, &mem->ls_res, &c_ptr);

	// blasfeo_mem align
    align_char_to(64, &c_ptr);

	// W_chol
    for (int ii = 0; ii <= N; ii++)
	{
		ny = cost_dims[ii]->ny;
        assign_blasfeo_dmat_mem(ny, ny, mem->W_chol+ii, &c_ptr);
	}
	// ls_res
    for (int ii = 0; ii <= N; ii++)
	{
		ny = cost_dims[ii]->ny;
        assign_blasfeo_dvec_mem(ny, mem->ls_res+ii, &c_ptr);
	}

	// dims
//    mem->dims = dims;

    assert((char *)raw_memory + ocp_nlp_gn_sqp_memory_calculate_size(config, dims, args) >= c_ptr);

    return mem;
}



/************************************************
* workspace
************************************************/



int ocp_nlp_gn_sqp_workspace_calculate_size(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, ocp_nlp_gn_sqp_opts *args)
{
	ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
//	sim_solver_config **sim_solvers = config->sim_solvers;
	ocp_nlp_dynamics_config **dynamics = config->dynamics;

	// loop index
	int ii;

	// extract dims
	int N = dims->N;
	int nb, ng;
	int nv, ny;

    int size = 0;

    size += sizeof(ocp_nlp_gn_sqp_work);

    size += ocp_qp_in_calculate_size(qp_solver, dims->qp_solver);
    size += ocp_qp_out_calculate_size(qp_solver, dims->qp_solver);
    size += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, args->qp_solver_opts);

    size += N*sizeof(sim_in *);
    size += N*sizeof(sim_out *);
    size += N*sizeof(void *);  // sim_work

    size += get_max_sim_workspace_size(config, dims, args);

    for (ii = 0; ii < N; ii++)
    {
        size += sim_in_calculate_size(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim);
        size += sim_out_calculate_size(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim);
    }


	// temporary stuff
    size += 2*(N+1)*sizeof(struct blasfeo_dvec); // tmp_ny, tmp_nbg
    size += 1*(N+1)*sizeof(struct blasfeo_dmat); // tmp_nv_ny
	size += 2*(N+1)*sizeof(double *); // ls_cost_in, ls_cost_jac_out

    for (ii = 0; ii < N+1; ii++)
    {
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
		nv = dims->cost[ii]->nx + dims->cost[ii]->nu;
		ny = dims->cost[ii]->ny;
        size += 1*blasfeo_memsize_dvec(ny); // tmp_ny
        size += 1*blasfeo_memsize_dvec(nb+ng); // tmp_nbg
        size += 1*blasfeo_memsize_dmat(nv, ny); // tmp_nv_ny
		size += 1*nv*sizeof(double); // ls_cost_in
		size += 1*(ny+ny*nv)*sizeof(double); // ls_cost_jac_out
    }

    size += 8;  // blasfeo_struct align
    size += 64;  // blasfeo_mem align

    return size;
}



static void ocp_nlp_gn_sqp_cast_workspace(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, ocp_nlp_gn_sqp_work *work, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_opts *args)
{
	ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
//	sim_solver_config **sim_solvers = config->sim_solvers;
	ocp_nlp_dynamics_config **dynamics = config->dynamics;

    char *c_ptr = (char *)work;
    c_ptr += sizeof(ocp_nlp_gn_sqp_work);

	// extract dims
    int N = dims->N;
	int nb, ng;
	int nv, ny;

	// blasfeo_struct align
    align_char_to(8, &c_ptr);

    // set up local SQP data
    assign_blasfeo_dvec_structs(N+1, &work->tmp_ny, &c_ptr);
    assign_blasfeo_dvec_structs(N+1, &work->tmp_nbg, &c_ptr);
    assign_blasfeo_dmat_structs(N+1, &work->tmp_nv_ny, &c_ptr);

	// ls_cost_in
	assign_double_ptrs(N+1, &work->ls_cost_in, &c_ptr);
	for (int ii=0; ii<=N; ii++)
	{
		nv = dims->cost[ii]->nx + dims->cost[ii]->nu;
		assign_double(nv, &work->ls_cost_in[ii], &c_ptr);
	}
	assign_double_ptrs(N+1, &work->ls_cost_jac_out, &c_ptr);
	for (int ii=0; ii<=N; ii++)
	{
		nv = dims->cost[ii]->nx + dims->cost[ii]->nu;
		ny = dims->cost[ii]->ny;
		assign_double(ny+ny+nv, &work->ls_cost_jac_out[ii], &c_ptr);
	}

	// blasfeo_mem align
    align_char_to(64, &c_ptr);

	// tmp_nv_ny
    for (int ii = 0; ii <= N; ii++)
	{
		nv = dims->cost[ii]->nx + dims->cost[ii]->nu;
		ny = dims->cost[ii]->ny;
        assign_blasfeo_dmat_mem(nv, ny, work->tmp_nv_ny+ii, &c_ptr);
	}
	// tmp_ny
    for (int ii = 0; ii <= N; ii++)
	{
		ny = dims->cost[ii]->ny;
        assign_blasfeo_dvec_mem(ny, work->tmp_ny+ii, &c_ptr);
	}
	// tmp_nbg
    for (int ii = 0; ii <= N; ii++)
	{
		nb = dims->constraints[ii]->nb;
		ng = dims->constraints[ii]->ng;
        assign_blasfeo_dvec_mem(nb+ng, work->tmp_nbg+ii, &c_ptr);
	}

    // set up QP solver
    work->qp_in = ocp_qp_in_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(qp_solver, dims->qp_solver);
    work->qp_out = ocp_qp_out_assign(qp_solver, dims->qp_solver, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(qp_solver, dims->qp_solver);

    work->qp_work = (void *)c_ptr;
    c_ptr += qp_solver->workspace_calculate_size(qp_solver, dims->qp_solver, args->qp_solver_opts);

    work->sim_in = (sim_in **) c_ptr;
    c_ptr += dims->N*sizeof(sim_in *);
    work->sim_out = (sim_out **) c_ptr;
    c_ptr += dims->N*sizeof(sim_out *);
    work->sim_solvers_work = (void **) c_ptr;
    c_ptr += dims->N*sizeof(void *);

    int max_sim_work_size = get_max_sim_workspace_size(config, dims, args);

    work->sim_solvers_work[0] = (void *)c_ptr;
    c_ptr += max_sim_work_size;

    for (int ii = 0; ii < dims->N; ii++)
    {
        work->sim_in[ii] = sim_in_assign(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim, c_ptr);
        c_ptr += sim_in_calculate_size(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim);
        work->sim_out[ii] = sim_out_assign(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim, c_ptr);
        c_ptr += sim_out_calculate_size(dynamics[ii]->sim_solver, dims->dynamics[ii]->sim);

        if (ii > 0) work->sim_solvers_work[ii] = work->sim_solvers_work[0];
    }

    assert((char *)work + ocp_nlp_gn_sqp_workspace_calculate_size(config, dims, args) >= c_ptr);
}



/************************************************
* solver
************************************************/



static void initialize_objective(ocp_nlp_dims *dims, ocp_nlp_in *nlp_in, ocp_nlp_gn_sqp_opts *args, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work)
{

    int N = nlp_in->dims->N;
	ocp_nlp_cost_dims **cost_dims = dims->cost;
	int nv;
	int ny;

    ocp_nlp_cost_ls_model **cost = (ocp_nlp_cost_ls_model**) nlp_in->cost;

	struct blasfeo_dmat *RSQrq = work->qp_in->RSQrq;

    for (int i = 0; i <= N; ++i)
	{

		nv = cost_dims[i]->nx + cost_dims[i]->nu;
		ny = cost_dims[i]->ny;
		// general Cyt

		// TODO recompute factorization only if W are re-tuned ???
		blasfeo_dpotrf_l(ny, &cost[i]->W, 0, 0, mem->W_chol+i, 0, 0);

		// linear ls
		// TODO avoid recomputing the Hessian if both W and Cyt do not change
		if (cost[i]->nls_mask==0)
		{
			blasfeo_dtrmm_rlnn(nv, ny, 1.0, mem->W_chol+i, 0, 0, &cost[i]->Cyt, 0, 0, work->tmp_nv_ny+i, 0, 0);
			blasfeo_dsyrk_ln(nv, ny, 1.0, work->tmp_nv_ny+i, 0, 0, work->tmp_nv_ny+i, 0, 0, 0.0, RSQrq+i, 0, 0, RSQrq+i, 0, 0);
		}

    }

	return;

}



static void linearize_update_qp_matrices(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_opts *args, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work)
{

	// loop index
	int i;

	// extract dims
    int N = nlp_in->dims->N;
	int nx, nu, nb, ng, nx1, nu1;
	int nv, ny;

    ocp_nlp_cost_ls_model **cost = (ocp_nlp_cost_ls_model **) nlp_in->cost;
	ocp_nlp_constraints_model **constraints = (ocp_nlp_constraints_model **) nlp_in->constraints;

    struct blasfeo_dmat *W_chol = mem->W_chol;
    struct blasfeo_dvec *ls_res = mem->ls_res;

    struct blasfeo_dvec *tmp_ny = work->tmp_ny;
    struct blasfeo_dvec *tmp_nbg = work->tmp_nbg;

    struct blasfeo_dmat *BAbt = work->qp_in->BAbt;
    struct blasfeo_dmat *RSQrq = work->qp_in->RSQrq;
	struct blasfeo_dmat *DCt = work->qp_in->DCt;

	ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    for (i = 0; i <= N; i++)
    {
		nv = dims->cost[i]->nx + dims->cost[i]->nu;
		ny = dims->cost[i]->ny;
		nx = dims->constraints[i]->nx;
		nu = dims->constraints[i]->nu;
		nb = dims->constraints[i]->nb;
		ng = dims->constraints[i]->ng;

		// dynamics

		// nlp mem: dyn_adj
		blasfeo_dvecse(nu+nx, 0.0, nlp_mem->dyn_adj+i, 0);

		if (i<N)
		{
			nx1 = dims->constraints[i+1]->nx;
			nu1 = dims->constraints[i+1]->nu;

			// pass state and control to integrator
			blasfeo_unpack_dvec(nu, nlp_out->ux+i, 0, work->sim_in[i]->u);
			blasfeo_unpack_dvec(nx, nlp_out->ux+i, nu, work->sim_in[i]->x);

			// call integrator
			config->dynamics[i]->sim_solver->evaluate(config->dynamics[i]->sim_solver, work->sim_in[i], work->sim_out[i], args->sim_solvers_opts[i],
				mem->sim_solvers_mem[i], work->sim_solvers_work[i]);

			// TODO(rien): transition functions for changing dimensions not yet implemented!

			// B
			blasfeo_pack_tran_dmat(nx1, nu, &work->sim_out[i]->S_forw[nx1*nx], nx1, &BAbt[i], 0, 0);
			// A
			blasfeo_pack_tran_dmat(nx1, nx, &work->sim_out[i]->S_forw[0], nx1, &BAbt[i], nu, 0);

			// nlp mem: dyn_fun
			blasfeo_pack_dvec(nx1, work->sim_out[i]->xn, nlp_mem->dyn_fun+i, 0);
			blasfeo_daxpy(nx1, -1.0, nlp_out->ux+i+1, nu1, nlp_mem->dyn_fun+i, 0, nlp_mem->dyn_fun+i, 0);
			// nlp mem: dyn_adj
			// TODO unless already computed in the simulation
			blasfeo_dgemv_n(nu+nx, nx1, -1.0, BAbt+i, 0, 0, nlp_out->pi+i, 0, 0.0, nlp_mem->dyn_adj+i, 0, nlp_mem->dyn_adj+i, 0);
		}

		// nlp mem: dyn_adj
		if(i>0)
			blasfeo_daxpy(nx, 1.0, nlp_out->pi+i-1, 0, nlp_mem->dyn_adj+i, nu, nlp_mem->dyn_adj+i, nu);



		// constraints
		// TODO merge dgemv_n and dgemv_t for general linear constraints

		// nlp_mem: ineq_fun
		blasfeo_dvecex_sp(nb, 1.0, constraints[i]->idxb, nlp_out->ux+i, 0, tmp_nbg+i, 0);
		blasfeo_dgemv_t(nu+nx, ng, 1.0, DCt+i, 0, 0, nlp_out->ux+i, 0, 0.0, tmp_nbg+i, nb, tmp_nbg+i, nb);
		blasfeo_daxpy(nb+ng, -1.0, tmp_nbg+i, 0, &constraints[i]->d, 0, nlp_mem->ineq_fun+i, 0);
		blasfeo_daxpy(nb+ng, -1.0, &constraints[i]->d, nb+ng, tmp_nbg+i, 0, nlp_mem->ineq_fun+i, nb+ng);

		// nlp_mem: ineq_adj
		blasfeo_dvecse(nu+nx, 0.0, nlp_mem->ineq_adj+i, 0);
		blasfeo_daxpy(nb+ng, -1.0, nlp_out->lam+i, nb+ng, nlp_out->lam+i, 0, tmp_nbg+i, 0);
		blasfeo_dvecad_sp(nb, 1.0, tmp_nbg+i, 0, constraints[i]->idxb, nlp_mem->ineq_adj+i, 0);
		blasfeo_dgemv_n(nu+nx, ng, 1.0, DCt+i, 0, 0, tmp_nbg+i, nb, 1.0, nlp_mem->ineq_adj+i, 0, nlp_mem->ineq_adj+i, 0);



        // cost
		// general Cyt

		if (cost[i]->nls_mask==0) // linear ls
		{

			blasfeo_dgemv_t(nv, ny, 1.0, &cost[i]->Cyt, 0, 0, nlp_out->ux+i, 0, -1.0, &cost[i]->y_ref, 0, ls_res+i, 0);

		}
		else // nonlinear ls
		{
			// unpack ls cost input
			blasfeo_unpack_dvec(nu, nlp_out->ux+i, 0, work->ls_cost_in[i]+nx);
			blasfeo_unpack_dvec(nx, nlp_out->ux+i, nu, work->ls_cost_in[i]);

			// evaluate external function (that assumes variables stacked as [x; u] )
			cost[i]->nls_jac->evaluate(cost[i]->nls_jac, work->ls_cost_in[i], work->ls_cost_jac_out[i]);

			// pack residuals into ls_res
			blasfeo_pack_dvec(ny, work->ls_cost_jac_out[i], ls_res+i, 0);
			// pack jacobian into Cyt
			blasfeo_pack_tran_dmat(ny, nx, work->ls_cost_jac_out[i]+ny, ny, &cost[i]->Cyt, nu, 0);
			blasfeo_pack_tran_dmat(ny, nu, work->ls_cost_jac_out[i]+ny+ny*nx, ny, &cost[i]->Cyt, 0, 0);

			blasfeo_daxpy(ny, -1.0, &cost[i]->y_ref, 0, ls_res+i, 0, ls_res+i, 0);

			blasfeo_dtrmm_rlnn(nv, ny, 1.0, W_chol+i, 0, 0, &cost[i]->Cyt, 0, 0, work->tmp_nv_ny+i, 0, 0);
			blasfeo_dsyrk_ln(nv, ny, 1.0, work->tmp_nv_ny+i, 0, 0, work->tmp_nv_ny+i, 0, 0, 0.0, RSQrq+i, 0, 0, RSQrq+i, 0, 0);
		}

		// nlp_mem: cost_grad

		// TODO use lower triangular chol of W to save n_y^2 flops
        blasfeo_dsymv_l(ny, ny, 1.0, &cost[i]->W, 0, 0, ls_res+i, 0, 0.0, tmp_ny+i, 0, tmp_ny+i, 0);
		blasfeo_dgemv_n(nv, ny, 1.0, &cost[i]->Cyt, 0, 0, tmp_ny+i, 0, 0.0, nlp_mem->cost_grad+i, 0, nlp_mem->cost_grad+i, 0);

		// TODO(rien) where should the update happen??? move to qp update ???
		if(i<N)
		{
			sim_rk_opts *opts = (sim_rk_opts*) args->sim_solvers_opts[i];
			if (opts->scheme != NULL && opts->scheme->type != exact)
			{
				for (int_t j = 0; j < nx; j++)
					DVECEL_LIBSTR(nlp_mem->cost_grad+i, nu+j) += work->sim_out[i]->grad[j];
				for (int_t j = 0; j < nu; j++)
					DVECEL_LIBSTR(nlp_mem->cost_grad+i, j) += work->sim_out[i]->grad[nx+j];
			}
		}

    }

	return;

}



// update QP rhs for SQP (step prim var, abs dual var)
static void sqp_update_qp_vectors(ocp_nlp_dims *dims, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_opts *args, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work)
{

	// loop index
	int i;

	// extract dims
    int N = nlp_in->dims->N;
	int nx, nu, nb, ng, nx1;

    struct blasfeo_dmat *RSQrq = work->qp_in->RSQrq;
    struct blasfeo_dvec *rq = work->qp_in->rq;
    struct blasfeo_dmat *BAbt = work->qp_in->BAbt;
    struct blasfeo_dvec *b = work->qp_in->b;
    struct blasfeo_dvec *d = work->qp_in->d;

	ocp_nlp_memory *nlp_mem = mem->nlp_mem;

	// g
	for (i=0; i<=N; i++)
	{

		nx = dims->constraints[i]->nx;
		nu = dims->constraints[i]->nu;
		blasfeo_dveccp(nu+nx, nlp_mem->cost_grad+i, 0, rq+i, 0);
        blasfeo_drowin(nu+nx, 1.0, rq+i, 0, RSQrq+i, nu+nx, 0); // XXX needed ???

	}

	// b
	for (i=0; i<N; i++)
	{
		nx = dims->constraints[i]->nx;
		nu = dims->constraints[i]->nu;
		nx1 = dims->constraints[i+1]->nx;
		blasfeo_dveccp(nx1, nlp_mem->dyn_fun+i, 0, b+i, 0);
		blasfeo_drowin(nx1, 1.0, b+i, 0, BAbt+i, nu+nx, 0); // XXX needed ???
	}

	// d
	for (i=0; i<=N; i++)
	{
		nb = dims->constraints[i]->nb;
		ng = dims->constraints[i]->ng;
		blasfeo_dveccp(2*nb+2*ng, nlp_mem->ineq_fun+i, 0, d+i, 0);
	}

	return;

}



static void update_variables(ocp_nlp_dims *dims, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_opts *args, ocp_nlp_gn_sqp_memory *mem, ocp_nlp_gn_sqp_work *work)
{

	// loop index
	int i, j;

	// extract dims
    int N = nlp_out->dims->N;
	int nx, nu, nb, ng, nx1;

    for (i = 0; i < N; i++)
    {
		nx1 = dims->constraints[i+1]->nx;
        for (j = 0; j < nx1; j++)
        {
            work->sim_in[i]->S_adj[j] = -DVECEL_LIBSTR(&work->qp_out->pi[i], j);
        }
    }

	// (full) step in primal variables
	for (i=0; i<=N; i++)
	{
		nx = dims->constraints[i]->nx;
		nu = dims->constraints[i]->nu;
		blasfeo_daxpy(nu+nx, 1.0, work->qp_out->ux+i, 0, nlp_out->ux+i, 0, nlp_out->ux+i, 0);
	}

	// absolute in dual variables
	for (i=0; i<N; i++)
	{
		nx1 = dims->constraints[i+1]->nx;
		blasfeo_dveccp(nx1, work->qp_out->pi+i, 0, nlp_out->pi+i, 0);
	}

	for (i=0; i<=N; i++)
	{
		nb = dims->constraints[i]->nb;
		ng = dims->constraints[i]->ng;
		blasfeo_dveccp(2*nb+2*ng, work->qp_out->lam+i, 0, nlp_out->lam+i, 0);
	}

	for (i=0; i<=N; i++)
	{
		nb = dims->constraints[i]->nb;
		ng = dims->constraints[i]->ng;
		blasfeo_dveccp(2*nb+2*ng, work->qp_out->t+i, 0, nlp_out->t+i, 0);
	}

	return;

}



// Simple fixed-step Gauss-Newton based SQP routine
int ocp_nlp_gn_sqp(ocp_nlp_solver_config *config, ocp_nlp_dims *dims, ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, ocp_nlp_gn_sqp_opts *opts, ocp_nlp_gn_sqp_memory *mem, void *work_)
{
	ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
	/* sim_solver_config **sim_solvers = config->sim_solvers; */

    ocp_nlp_gn_sqp_work *work = (ocp_nlp_gn_sqp_work*) work_;
    ocp_nlp_gn_sqp_cast_workspace(config, dims, work, mem, opts);

    int N = dims->N;
	int nx, nu, nx1;

    sim_rk_opts *sim_opts;

    // set up integrators

    for (int ii = 0; ii < N; ii++)
    {
		nx = dims->constraints[ii]->nx;
		nx1 = dims->constraints[ii+1]->nx;
		nu = dims->constraints[ii]->nu;

		ocp_nlp_dynamics_model *dynamics = nlp_in->dynamics[ii];
		work->sim_in[ii]->model = dynamics->sim_model;

        sim_opts = opts->sim_solvers_opts[ii];
        work->sim_in[ii]->T = nlp_in->Ts[ii];

        // TODO(dimitris): REVISE IF THIS IS CORRECT FOR VARYING DIMENSIONS!
        for (int jj = 0; jj < nx1 * (nx + nu); jj++)
            work->sim_in[ii]->S_forw[jj] = 0.0;
        for (int jj = 0; jj < nx1; jj++)
            work->sim_in[ii]->S_forw[jj * (nx + 1)] = 1.0;
        for (int jj = 0; jj < nx + nu; jj++)
            work->sim_in[ii]->S_adj[jj] = 0.0;
        // for (int jj = 0; jj < nlp_in->dims->num_stages[ii] * nx[ii+1]; jj++)
            // work->sim_in[ii]->grad_K[jj] = 0.0;
    }

    initialize_objective(dims, nlp_in, opts, mem, work);

	// initialize constraints
	for (int ii=0; ii<=N; ii++)
		config->constraints[ii]->initialize_qp(config->constraints[ii], dims->constraints[ii], nlp_in->constraints[ii], work->qp_in->idxb[ii], work->qp_in->DCt+ii, NULL, NULL);

	// start timer
    acados_timer timer;
    real_t total_time = 0;
    acados_tic(&timer);

	// main sqp loop
    int max_sqp_iterations =  opts->maxIter;
	int sqp_iter = 0;
    for ( ; sqp_iter < max_sqp_iterations; sqp_iter++)
    {

        linearize_update_qp_matrices(config, dims, nlp_in, nlp_out, opts, mem, work);

		// update QP rhs for SQP (step prim var, abs dual var)
        sqp_update_qp_vectors(dims, nlp_in, nlp_out, opts, mem, work);

		// compute nlp residuals
		ocp_nlp_res_compute(dims, nlp_in, nlp_out, mem->nlp_res, mem->nlp_mem);

		// TODO exit conditions on residuals
		if( (mem->nlp_res->inf_norm_res_g < opts->min_res_g) &
			(mem->nlp_res->inf_norm_res_b < opts->min_res_b) &
			(mem->nlp_res->inf_norm_res_d < opts->min_res_d) &
			(mem->nlp_res->inf_norm_res_m < opts->min_res_m) )
		{

			// save sqp iterations number
			mem->sqp_iter = sqp_iter;

			// stop timer
			total_time += acados_toc(&timer);

			return 0;

		}

//print_ocp_qp_in(work->qp_in);
//exit(1);

        int_t qp_status = qp_solver->evaluate(qp_solver, work->qp_in, work->qp_out,
            opts->qp_solver_opts, mem->qp_solver_mem, work->qp_work);

//print_ocp_qp_out(work->qp_out);
//exit(1);

        if (qp_status != 0)
        {
            printf("QP solver returned error status %d\n", qp_status);
            return -1;
        }

        update_variables(dims, nlp_out, opts, mem, work);

//ocp_nlp_dims_print(nlp_out->dims);
//ocp_nlp_out_print(nlp_out);
//exit(1);

        for (int_t i = 0; i < N; i++)
        {
            sim_rk_opts *rk_opts = (sim_rk_opts*) opts->sim_solvers_opts[i];
            if (rk_opts->scheme == NULL)
                continue;
            rk_opts->sens_adj = (rk_opts->scheme->type != exact);
            if (nlp_in->freezeSens) {
                // freeze inexact sensitivities after first SQP iteration !!
                rk_opts->scheme->freeze = true;
            }
        }

    }

	// stop timer
    total_time += acados_toc(&timer);

//	ocp_nlp_out_print(nlp_out);

	// save sqp iterations number
	mem->sqp_iter = sqp_iter;

	// maximum number of iterations reached
    return 1;

}
