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

// external
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// acados
// TODO(dimitris): remove most includes
#include "acados/sim/sim_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_irk_integrator.h"
#include "acados/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/external_function_generic.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// blasfeo
#include <blasfeo/include/blasfeo_target.h>
#include <blasfeo/include/blasfeo_common.h>
#include <blasfeo/include/blasfeo_d_aux.h>
#include <blasfeo/include/blasfeo_d_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_v_aux_ext_dep.h>
#include <blasfeo/include/blasfeo_d_blas.h>

// wt model
#include "examples/c/wt_model_nx6/wt_model.h"

// x0 and u for simulation
#include "examples/c/wt_model_nx6/u_x0.c"



int main()
{

	/************************************************
	* initialization
	************************************************/

    int ii, jj;

    int nx = 6;
    int nu = 2;
	int np = 1;

    int NF = nx + nu; // columns of forward seed

    double Ts = 0.2; // simulation time

	double *x_sim = malloc(sizeof(double)*nx*(nsim+1));

	for (ii=0; ii<nx; ii++)
		x_sim[ii] = x_ref[ii];

	/************************************************
	* external functions (explicit model)
	************************************************/

	// expl_ode_fun
	external_function_param_casadi expl_ode_fun;
	expl_ode_fun.casadi_fun = &casadi_expl_ode_fun;
	expl_ode_fun.casadi_work = &casadi_expl_ode_fun_work;
	expl_ode_fun.casadi_sparsity_in = &casadi_expl_ode_fun_sparsity_in;
	expl_ode_fun.casadi_sparsity_out = &casadi_expl_ode_fun_sparsity_out;
	expl_ode_fun.casadi_n_in = &casadi_expl_ode_fun_n_in;
	expl_ode_fun.casadi_n_out = &casadi_expl_ode_fun_n_out;
	external_function_param_casadi_create(&expl_ode_fun, np);

	// expl_vde_for
	external_function_param_casadi expl_vde_for;
	expl_vde_for.casadi_fun = &casadi_expl_vde_for;
	expl_vde_for.casadi_work = &casadi_expl_vde_for_work;
	expl_vde_for.casadi_sparsity_in = &casadi_expl_vde_for_sparsity_in;
	expl_vde_for.casadi_sparsity_out = &casadi_expl_vde_for_sparsity_out;
	expl_vde_for.casadi_n_in = &casadi_expl_vde_for_n_in;
	expl_vde_for.casadi_n_out = &casadi_expl_vde_for_n_out;
	external_function_param_casadi_create(&expl_vde_for, np);

	// expl_vde_adj
	external_function_param_casadi expl_vde_adj;
	expl_vde_adj.casadi_fun = &casadi_expl_vde_adj;
	expl_vde_adj.casadi_work = &casadi_expl_vde_adj_work;
	expl_vde_adj.casadi_sparsity_in = &casadi_expl_vde_adj_sparsity_in;
	expl_vde_adj.casadi_sparsity_out = &casadi_expl_vde_adj_sparsity_out;
	expl_vde_adj.casadi_n_in = &casadi_expl_vde_adj_n_in;
	expl_vde_adj.casadi_n_out = &casadi_expl_vde_adj_n_out;
	external_function_param_casadi_create(&expl_vde_adj, np);

	/************************************************
	* external functions (implicit model)
	************************************************/

	// impl_ode_fun
	external_function_param_casadi impl_ode_fun;
	impl_ode_fun.casadi_fun = &casadi_impl_ode_fun;
	impl_ode_fun.casadi_work = &casadi_impl_ode_fun_work;
	impl_ode_fun.casadi_sparsity_in = &casadi_impl_ode_fun_sparsity_in;
	impl_ode_fun.casadi_sparsity_out = &casadi_impl_ode_fun_sparsity_out;
	impl_ode_fun.casadi_n_in = &casadi_impl_ode_fun_n_in;
	impl_ode_fun.casadi_n_out = &casadi_impl_ode_fun_n_out;
	external_function_param_casadi_create(&impl_ode_fun, np);

	// impl_ode_jac_x
	external_function_param_casadi impl_ode_jac_x;
	impl_ode_jac_x.casadi_fun = &casadi_impl_ode_jac_x;
	impl_ode_jac_x.casadi_work = &casadi_impl_ode_jac_x_work;
	impl_ode_jac_x.casadi_sparsity_in = &casadi_impl_ode_jac_x_sparsity_in;
	impl_ode_jac_x.casadi_sparsity_out = &casadi_impl_ode_jac_x_sparsity_out;
	impl_ode_jac_x.casadi_n_in = &casadi_impl_ode_jac_x_n_in;
	impl_ode_jac_x.casadi_n_out = &casadi_impl_ode_jac_x_n_out;
	external_function_param_casadi_create(&impl_ode_jac_x, np);

	// impl_ode_jac_xdot
	external_function_param_casadi impl_ode_jac_xdot;
	impl_ode_jac_xdot.casadi_fun = &casadi_impl_ode_jac_xdot;
	impl_ode_jac_xdot.casadi_work = &casadi_impl_ode_jac_xdot_work;
	impl_ode_jac_xdot.casadi_sparsity_in = &casadi_impl_ode_jac_xdot_sparsity_in;
	impl_ode_jac_xdot.casadi_sparsity_out = &casadi_impl_ode_jac_xdot_sparsity_out;
	impl_ode_jac_xdot.casadi_n_in = &casadi_impl_ode_jac_xdot_n_in;
	impl_ode_jac_xdot.casadi_n_out = &casadi_impl_ode_jac_xdot_n_out;
	external_function_param_casadi_create(&impl_ode_jac_xdot, np);

	// impl_ode_jac_u
	external_function_param_casadi impl_ode_jac_u;
	impl_ode_jac_u.casadi_fun = &casadi_impl_ode_jac_u;
	impl_ode_jac_u.casadi_work = &casadi_impl_ode_jac_u_work;
	impl_ode_jac_u.casadi_sparsity_in = &casadi_impl_ode_jac_u_sparsity_in;
	impl_ode_jac_u.casadi_sparsity_out = &casadi_impl_ode_jac_u_sparsity_out;
	impl_ode_jac_u.casadi_n_in = &casadi_impl_ode_jac_u_n_in;
	impl_ode_jac_u.casadi_n_out = &casadi_impl_ode_jac_u_n_out;
	external_function_param_casadi_create(&impl_ode_jac_u, np);

	// impl_ode_fun_jac_x_xdot
	external_function_param_casadi impl_ode_fun_jac_x_xdot;
	impl_ode_fun_jac_x_xdot.casadi_fun = &casadi_impl_ode_fun_jac_x_xdot;
	impl_ode_fun_jac_x_xdot.casadi_work = &casadi_impl_ode_fun_jac_x_xdot_work;
	impl_ode_fun_jac_x_xdot.casadi_sparsity_in = &casadi_impl_ode_fun_jac_x_xdot_sparsity_in;
	impl_ode_fun_jac_x_xdot.casadi_sparsity_out = &casadi_impl_ode_fun_jac_x_xdot_sparsity_out;
	impl_ode_fun_jac_x_xdot.casadi_n_in = &casadi_impl_ode_fun_jac_x_xdot_n_in;
	impl_ode_fun_jac_x_xdot.casadi_n_out = &casadi_impl_ode_fun_jac_x_xdot_n_out;
	external_function_param_casadi_create(&impl_ode_fun_jac_x_xdot, np);

	// impl_ode_jac_x_xdot_u
	external_function_param_casadi impl_ode_jac_x_xdot_u;
	impl_ode_jac_x_xdot_u.casadi_fun = &casadi_impl_ode_jac_x_xdot_u;
	impl_ode_jac_x_xdot_u.casadi_work = &casadi_impl_ode_jac_x_xdot_u_work;
	impl_ode_jac_x_xdot_u.casadi_sparsity_in = &casadi_impl_ode_jac_x_xdot_u_sparsity_in;
	impl_ode_jac_x_xdot_u.casadi_sparsity_out = &casadi_impl_ode_jac_x_xdot_u_sparsity_out;
	impl_ode_jac_x_xdot_u.casadi_n_in = &casadi_impl_ode_jac_x_xdot_u_n_in;
	impl_ode_jac_x_xdot_u.casadi_n_out = &casadi_impl_ode_jac_x_xdot_u_n_out;
	external_function_param_casadi_create(&impl_ode_jac_x_xdot_u, np);

	// impl_ode_jac_x_u
	external_function_param_casadi impl_ode_jac_x_u;
	impl_ode_jac_x_u.casadi_fun = &casadi_impl_ode_jac_x_u;
	impl_ode_jac_x_u.casadi_work = &casadi_impl_ode_jac_x_u_work;
	impl_ode_jac_x_u.casadi_sparsity_in = &casadi_impl_ode_jac_x_u_sparsity_in;
	impl_ode_jac_x_u.casadi_sparsity_out = &casadi_impl_ode_jac_x_u_sparsity_out;
	impl_ode_jac_x_u.casadi_n_in = &casadi_impl_ode_jac_x_u_n_in;
	impl_ode_jac_x_u.casadi_n_out = &casadi_impl_ode_jac_x_u_n_out;
	external_function_param_casadi_create(&impl_ode_jac_x_u, np);




	int number_sim_solvers = 3;
	int nss;
	for (nss = 0; nss < number_sim_solvers; nss++)
	{
		/************************************************
		* sim plan & config
		************************************************/

		// choose plan
		sim_solver_plan plan;
		switch (nss)
		{

			case 0:
				plan.sim_solver = ERK;
				break;

			case 1:
			case 2:
				plan.sim_solver = IRK;
				break;

			default :
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);

		}


		// create correct config based on plan
		sim_solver_config *config = sim_config_create(plan);

		/************************************************
		* sim dims
		************************************************/

		void *dims = sim_dims_create(config);
		config->set_nx(dims, nx);
		config->set_nu(dims, nu);

		/************************************************
		* sim opts
		************************************************/

		sim_rk_opts *opts = sim_opts_create(config, dims);

	//		opts->ns = 4; // number of stages in rk integrator
	//		opts->num_steps = 5; // number of integration steps
		opts->sens_adj = true;
		opts->sens_forw = true;


		switch (nss)
		{

			case 0:
				opts->ns = 4; // number of stages in rk integrator
				opts->num_steps = 10; // number of integration steps
				break;

			case 1:
				opts->ns = 2; // number of stages in rk integrator
				opts->num_steps = 6; // number of integration steps
				break;

			case 2:
				opts->ns = 8; // number of stages in rk integrator
				opts->num_steps = 3; // number of integration steps
				break;

			default :
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);

		}

		/************************************************
		* sim in / out
		************************************************/

		sim_in *in = sim_in_create(config, dims);
		sim_out *out = sim_out_create(config, dims);

		in->T = Ts;

		// external functions
		switch (nss)
		{
			case 0:
			{
				sim_set_model(config, in, "expl_ode_fun", &expl_ode_fun);
				sim_set_model(config, in, "expl_vde_for", &expl_vde_for);
				sim_set_model(config, in, "expl_vde_adj", &expl_vde_adj);
				break;
			}
			case 1:
			case 2:
			{
				sim_set_model(config, in, "impl_ode_fun", &impl_ode_fun);
				sim_set_model(config, in, "impl_ode_fun_jac_x_xdot", &impl_ode_fun_jac_x_xdot);
				sim_set_model(config, in, "impl_ode_jac_x_xdot_u", &impl_ode_jac_x_xdot_u);
				sim_set_model(config, in, "impl_ode_jac_x_u", &impl_ode_jac_x_u);
				break;
			}
			default :
			{
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);
			}
		}

		// seeds forw
		for (ii = 0; ii < nx * NF; ii++)
			in->S_forw[ii] = 0.0;
		for (ii = 0; ii < nx; ii++)
			in->S_forw[ii * (nx + 1)] = 1.0;

		// seeds adj
		for (ii = 0; ii < nx; ii++)
			in->S_adj[ii] = 1.0;

		/************************************************
		* sim solver
		************************************************/

		// print solver info
		switch (nss)
		{

			case 0:
				printf("\n\nsim solver: ERK, ns=%d, num_steps=%d\n", opts->ns, opts->num_steps);
				plan.sim_solver = ERK;
				break;

			case 1:
			case 2:
				printf("\n\nsim solver: IRK, ns=%d, num_steps=%d\n", opts->ns, opts->num_steps);
				plan.sim_solver = IRK;
				break;

			default :
				printf("\nnot enough sim solvers implemented!\n");
				exit(1);

		}

		sim_solver *sim_solver = sim_create(config, dims, opts);

		int acados_return;

		acados_timer timer;
		acados_tic(&timer);

		int nsim0 = nsim;

		double cpu_time = 0.0;
		double la_time = 0.0;
		double ad_time = 0.0;

		// to avoid unstable behavior introduce a small pi-controller for rotor speed tracking
		double uctrl = 0.0;
		double uctrlI = 0.0;
		double kI = 1e-1;
		double kP = 10;
		double tmp, ctrlErr;

		for (ii=0; ii<nsim; ii++)
		{
			// update initial state
			for (jj = 0; jj < nx; jj++)
				in->x[jj] = x_sim[ii*nx+jj];

			// compute inputs
			for (jj = 0; jj < nu; jj++)
				in->u[jj] = u_sim[ii*nu+jj];
			tmp = in->u[1] - uctrl;
			in->u[1] = tmp>0.0 ? tmp : 0.0;

			// update parameters
			switch (nss)
			{
				case 0:
				{
					expl_ode_fun.set_param(&expl_ode_fun, p_sim+ii*np);
					expl_vde_for.set_param(&expl_vde_for, p_sim+ii*np);
					expl_vde_for.set_param(&expl_vde_adj, p_sim+ii*np);
					break;
				}
				case 1:
				case 2:
				{
					impl_ode_fun.set_param(&impl_ode_fun, p_sim+ii*np);
					impl_ode_fun_jac_x_xdot.set_param(&impl_ode_fun_jac_x_xdot, p_sim+ii*np);
					impl_ode_jac_x_xdot_u.set_param(&impl_ode_jac_x_xdot_u, p_sim+ii*np);
					impl_ode_jac_x_u.set_param(&impl_ode_jac_x_u, p_sim+ii*np);
					break;
				}
				default :
				{
					printf("\nnot enough sim solvers implemented!\n");
					exit(1);
				}
			}



			// d_print_mat(1, nx, in->x, 1);
			// d_print_mat(1, nu, in->u, 1);

			// execute simulation step with current input and state
			acados_return = sim_solve(sim_solver, in, out);
			if (acados_return != 0)
			{
				printf("error in sim solver\n");
				return ACADOS_FAILURE;
			}

			cpu_time += out->info->CPUtime;
			la_time += out->info->LAtime;
			ad_time += out->info->ADtime;

			// d_print_mat(1, nx, out->xn, 1);
			// d_print_mat(1, nx, x_ref+ii*nx, 1);

			// extract state at next time step
			for (jj = 0; jj < nx; jj++)
				x_sim[(ii+1)*nx+jj] = out->xn[jj];

			// update PI-controller
			ctrlErr = x_ref[nx*(ii+1)] - x_sim[nx*(ii+1)];
			uctrlI = uctrlI + kI*ctrlErr*Ts;
			uctrl = kP*ctrlErr + uctrlI;

			// if (ii < nsim-1)
			// 	printf("\nii = %d, sim error = %e\n", ii, ctrlErr);
		}
		double total_cpu_time = acados_toc(&timer);

		/************************************************
		* printing
		************************************************/

		printf("\nxn: \n");
		// for (ii=0; ii<nx; ii++)
		// 	printf("%8.5f ", x_sim[nsim0*nx+ii]);
		// printf("\n");
		d_print_e_mat(1, nx, &x_sim[nsim0*nx], 1);

		double *S_forw_out;
		S_forw_out = NULL;
		if(opts->sens_forw){
			S_forw_out = out->S_forw;
			printf("\nS_forw_out: \n");
			d_print_e_mat(nx, NF, S_forw_out, nx);
			// for (ii=0;ii<nx;ii++){
			// 	for (jj=0;jj<NF;jj++)
			// 		printf("%8.5f ", S_forw_out[jj*nx+ii]);
			// 	printf("\n");
			// }
		}


		if(opts->sens_adj){
			double *S_adj_out = out->S_adj;
			printf("\nS_adj_out: \n");
			d_print_e_mat(1, nx+nu, S_adj_out, 1);
		}

if(opts->sens_forw){		// debug adjoints
      struct blasfeo_dmat S_forw_result;
      struct blasfeo_dvec adjoint_seed;
      struct blasfeo_dvec forw_times_seed;

      int Sf_mem_size = blasfeo_memsize_dmat(nx, nx+nu);
      int adj_s_mem_size = blasfeo_memsize_dvec(nx);
      int check_mem_size = blasfeo_memsize_dvec(nx+nu);

      void *Sf_mem; v_zeros_align(&Sf_mem, Sf_mem_size);
      void *seed_mem; v_zeros_align(&seed_mem, adj_s_mem_size);
      void *check_mem; v_zeros_align(&check_mem, check_mem_size);

      blasfeo_create_dmat(nx, nu+nx, &S_forw_result, Sf_mem);
      blasfeo_create_dvec(nx, &adjoint_seed, seed_mem);
      blasfeo_create_dvec(nu+nx, &forw_times_seed, check_mem);

      blasfeo_pack_dmat(nx, nx+nu, S_forw_out, nx, &S_forw_result, 0, 0);
      blasfeo_pack_dvec(nx, in->S_adj, &adjoint_seed, 0);

      blasfeo_dgemv_t(nx, nx+nu, 1.0, &S_forw_result, 0, 0, &adjoint_seed, 0, 0.0, &forw_times_seed, 0, &forw_times_seed, 0);
      printf("S_forw^T * adj_seed = \n");
      blasfeo_print_exp_tran_dvec(nx+nu, &forw_times_seed, 0);

      v_free_align(Sf_mem);
      v_free_align(seed_mem);
      v_free_align(check_mem);
		}

    #if 0
		printf("\n");
		printf("cpt: %8.4f [ms]\n", 1000*out->info->CPUtime);
		printf("AD cpt: %8.4f [ms]\n", 1000*out->info->ADtime);

	#endif

		// printf("time split: %f ms CPU, %f ms LA, %f ms AD\n\n", cpu_time, la_time, ad_time);
		printf("\n\ntime for %d simulation steps: %f ms (AD time: %f ms (%5.2f%%))\n\n", nsim, 1e3*total_cpu_time, 1e3*ad_time, 1e2*ad_time/cpu_time);

		/************************************************
		* free memory
		************************************************/

		free(sim_solver);
		free(in);
		free(out);

		free(opts);
		free(config);

	}

	free(x_sim);

	// explicit model
	external_function_param_casadi_free(&expl_ode_fun);
	external_function_param_casadi_free(&expl_vde_for);
	external_function_param_casadi_free(&expl_vde_adj);
	// implicit model
	external_function_param_casadi_free(&impl_ode_fun);
	external_function_param_casadi_free(&impl_ode_jac_x);
	external_function_param_casadi_free(&impl_ode_jac_xdot);
	external_function_param_casadi_free(&impl_ode_jac_u);
	external_function_param_casadi_free(&impl_ode_fun_jac_x_xdot);
	external_function_param_casadi_free(&impl_ode_jac_x_xdot_u);
	external_function_param_casadi_free(&impl_ode_jac_x_u);

	printf("\nsuccess!\n\n");

    return 0;
}