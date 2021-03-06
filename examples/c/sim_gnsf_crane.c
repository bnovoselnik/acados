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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// acados
#include <acados/sim/sim_common.h>
#include <acados/sim/sim_gnsf.h>

#include "acados/utils/external_function_generic.h"
#include "acados/utils/math.h"

#include "acados_c/external_function_interface.h"
#include "acados_c/sim_interface.h"

// model
#include "examples/c/gnsf_crane_model/gnsf_crane_model.h"

int main() {

/************************************************
*   external functions
************************************************/

    // phi_fun
    external_function_casadi phi_fun;
    phi_fun.casadi_fun            = &casadi_phi_fun;
    phi_fun.casadi_work           = &casadi_phi_fun_work;
    phi_fun.casadi_sparsity_in    = &casadi_phi_fun_sparsity_in;
    phi_fun.casadi_sparsity_out   = &casadi_phi_fun_sparsity_out;
    phi_fun.casadi_n_in           = &casadi_phi_fun_n_in;
    phi_fun.casadi_n_out          = &casadi_phi_fun_n_out;
	external_function_casadi_create(&phi_fun);

    // phi_fun_jac_y
    external_function_casadi phi_fun_jac_y;
    phi_fun_jac_y.casadi_fun            = &casadi_phi_fun_jac_y;
    phi_fun_jac_y.casadi_work           = &casadi_phi_fun_jac_y_work;
    phi_fun_jac_y.casadi_sparsity_in    = &casadi_phi_fun_jac_y_sparsity_in;
    phi_fun_jac_y.casadi_sparsity_out   = &casadi_phi_fun_jac_y_sparsity_out;
    phi_fun_jac_y.casadi_n_in           = &casadi_phi_fun_jac_y_n_in;
    phi_fun_jac_y.casadi_n_out          = &casadi_phi_fun_jac_y_n_out;
	external_function_casadi_create(&phi_fun_jac_y);

    // phi_jac_y_uhat
    external_function_casadi phi_jac_y_uhat;
    phi_jac_y_uhat.casadi_fun                = &casadi_phi_jac_y_uhat;
    phi_jac_y_uhat.casadi_work               = &casadi_phi_jac_y_uhat_work;
    phi_jac_y_uhat.casadi_sparsity_in        = &casadi_phi_jac_y_uhat_sparsity_in;
    phi_jac_y_uhat.casadi_sparsity_out       = &casadi_phi_jac_y_uhat_sparsity_out;
    phi_jac_y_uhat.casadi_n_in               = &casadi_phi_jac_y_uhat_n_in;
    phi_jac_y_uhat.casadi_n_out              = &casadi_phi_jac_y_uhat_n_out;

	external_function_casadi_create(&phi_jac_y_uhat);

    // f_lo_fun_jac_x1k1uz
    external_function_casadi f_lo_fun_jac_x1k1uz;
    f_lo_fun_jac_x1k1uz.casadi_fun            = &casadi_f_lo_fun_jac_x1k1uz;
    f_lo_fun_jac_x1k1uz.casadi_work           = &casadi_f_lo_fun_jac_x1k1uz_work;
    f_lo_fun_jac_x1k1uz.casadi_sparsity_in    = &casadi_f_lo_fun_jac_x1k1uz_sparsity_in;
    f_lo_fun_jac_x1k1uz.casadi_sparsity_out   = &casadi_f_lo_fun_jac_x1k1uz_sparsity_out;
    f_lo_fun_jac_x1k1uz.casadi_n_in           = &casadi_f_lo_fun_jac_x1k1uz_n_in;
    f_lo_fun_jac_x1k1uz.casadi_n_out          = &casadi_f_lo_fun_jac_x1k1uz_n_out;
	external_function_casadi_create(&f_lo_fun_jac_x1k1uz);

    // get_matrices_fun
    external_function_casadi get_matrices_fun;
    get_matrices_fun.casadi_fun            = &casadi_get_matrices_fun;
    get_matrices_fun.casadi_work           = &casadi_get_matrices_fun_work;
    get_matrices_fun.casadi_sparsity_in    = &casadi_get_matrices_fun_sparsity_in;
    get_matrices_fun.casadi_sparsity_out   = &casadi_get_matrices_fun_sparsity_out;
    get_matrices_fun.casadi_n_in           = &casadi_get_matrices_fun_n_in;
    get_matrices_fun.casadi_n_out          = &casadi_get_matrices_fun_n_out;
	external_function_casadi_create(&get_matrices_fun);

/************************************************
* Set up sim_gnsf structs
************************************************/

    sim_solver_plan plan;

    plan.sim_solver = GNSF;
    sim_solver_config *config = sim_config_create(plan);

    void *dims = sim_dims_create(config);
    sim_gnsf_dims *gnsf_dim = (sim_gnsf_dims *) dims;
    int nx = 9;
    int nu = 2;
    int num_stages = 3;
    int num_steps = 1;

    gnsf_dim->nx = nx;
    gnsf_dim->nx1 = 8;
    gnsf_dim->nx2 = 1;
    gnsf_dim->nu = nu;
    gnsf_dim->n_out = 2;
    gnsf_dim->nz = 1;
    gnsf_dim->ny = 4;
    gnsf_dim->nuhat = 1;


    // set up sim_dims
    // sim_dims *dims = (sim_dims *) gnsf_dim; // typecasting works as gnsf_dims has entries of sim_dims at the beginning

    // set up gnsf_opts
    int opts_size = config->opts_calculate_size(config, dims);
	void *opts_mem = malloc(opts_size);
    sim_rk_opts *opts = config->opts_assign(config, dims, opts_mem);
    config->opts_initialize_default(config, dims, opts);
    opts->sens_adj = true;

    // opts->A_mat[0] = 8.696121837661337e-02;
    // opts->A_mat[1] = 1.881133909077607e-01;
    // opts->A_mat[2] = 1.671849744614559e-01;
    // opts->A_mat[3] = 1.774751870097904e-01;
    // opts->A_mat[4] = -2.660256900501954e-02;
    // opts->A_mat[5] = 1.630292692796755e-01;
    // opts->A_mat[6] = 3.539526857259798e-01;
    // opts->A_mat[7] = 3.134484623618029e-01;
    // opts->A_mat[8] = 1.262633154623132e-02;
    // opts->A_mat[9] = -2.787789181794523e-02;
    // opts->A_mat[10] = 1.630455246283593e-01;
    // opts->A_mat[11] = 3.526773629130557e-01;
    // opts->A_mat[12] = -3.554980917825141e-03;
    // opts->A_mat[13] = 6.735231630508876e-03;
    // opts->A_mat[14] = -1.419318481579561e-02;
    // opts->A_mat[15] = 8.695898771535188e-02;

    // opts->b_vec[0] =     1.739190668368515e-01;
    // opts->b_vec[1] =     3.260788016404259e-01;
    // opts->b_vec[2] =     3.260666571808232e-01;
    // opts->b_vec[3] =     1.739354743418989e-01;
    
    opts->ns = num_stages;
    opts->num_steps = num_steps;
    
    sim_in *in = sim_in_create(config, dims);

    for (int ii = 0; ii < nx *(nx +nu); ii++) {
        in->S_forw[ii] = 0.0;
    }
    for (int ii = 0; ii < nx; ii++) {
        in->S_forw[ii+ ii*nx] = 1.0;
        in->x[ii] = 0.0;
    }
    for (int ii = 0; ii < nx; ii++) {
        in->S_adj[ii] = 1.0;
    }
    in->x[2] = 0.8;
    in->u[0] = 40.108149413030752;
    in->u[1] = -50.446662212534974;
    in->T = 0.1;

    // set up gnsf_model
    gnsf_model *model = in->model;
    // set external functions
    model->f_lo_fun_jac_x1_x1dot_u_z = (external_function_generic *) &f_lo_fun_jac_x1k1uz;
    model->phi_fun = (external_function_generic *) &phi_fun;
    model->phi_fun_jac_y = (external_function_generic *) &phi_fun_jac_y;
    model->phi_jac_y_uhat = (external_function_generic *) &phi_jac_y_uhat;
    external_function_generic *get_model_matrices = (external_function_generic *) &get_matrices_fun;
    sim_gnsf_import_matrices(gnsf_dim, model, get_model_matrices);

    sim_solver *sim_solver = sim_create(config, dims, opts);

    sim_gnsf_precompute(config, gnsf_dim, model, opts, sim_solver->mem, sim_solver->work, in->T);

    // gnsf_import_precomputed(gnsf_dim, model, But_KK_YY_ZZ_LO_fun);

    // set up sim_out
    sim_out *out = sim_out_create(config, dims);

    int NREP = 100;
    double casadi_times[NREP];
    double gnsf_times[NREP];
    for (int i = 0; i < NREP; i++) {
        int acados_return = sim_solve(sim_solver, in, out);
        if (acados_return != 0)
        {
            printf("error in sim solver\n");
            return ACADOS_FAILURE;
        }
        casadi_times[i] = out->info->ADtime;
        gnsf_times[i] = out->info->CPUtime;
    }

    double casadi_time = minimum_of_doubles(casadi_times, NREP);
    double gnsf_time = minimum_of_doubles(gnsf_times, NREP);

    // PRINTING
    printf("Newton_iter = %d,\t num_steps = %d \n", opts->newton_iter, opts->num_steps);
    printf("xf =\n");
    d_print_e_mat(1, nx, out->xn, 1);
    printf("forw_Sensitivities = \n");
    d_print_e_mat(nx, nx + nu, out->S_forw, nx);
    printf("adj Sensitivities =\n");
    d_print_e_mat(1, nx + nu, out->S_adj, 1);
    
    printf("gnsf_time  =  %f [ms]\n", gnsf_time*1000);
    printf("casadi_time =  %f [ms]\t minimum of %d executions \n", casadi_time*1000, NREP);

    free(config);
    free(gnsf_dim);
    free(in);
    free(out);
    free(opts_mem);
    free(sim_solver);

	external_function_casadi_free(&phi_fun);
	external_function_casadi_free(&phi_fun_jac_y);
	external_function_casadi_free(&phi_jac_y_uhat);
    external_function_casadi_free(&f_lo_fun_jac_x1k1uz);
	external_function_casadi_free(&get_matrices_fun);

    return 0;
}
