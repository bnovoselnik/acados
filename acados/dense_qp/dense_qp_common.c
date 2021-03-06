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
#include <assert.h>
#include <stddef.h>
// hpipm
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_ipm.h"
#include "hpipm/include/hpipm_d_dense_qp_kkt.h"
#include "hpipm/include/hpipm_d_dense_qp_res.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/dense_qp/dense_qp_common.h"
#include "acados/utils/types.h"

/************************************************
 * config
 ************************************************/

int dense_qp_solver_config_calculate_size()
{
    int size = 0;

    size += sizeof(qp_solver_config);

    return size;
}

qp_solver_config *dense_qp_solver_config_assign(void *raw_memory)
{
    char *c_ptr = raw_memory;

    qp_solver_config *config = (qp_solver_config *) c_ptr;
    c_ptr += sizeof(qp_solver_config);

    return config;
}

/************************************************
 * dims
 ************************************************/

int dense_qp_dims_calculate_size()
{
    int size = sizeof(dense_qp_dims);

    size += d_memsize_dense_qp_dim();

    return size;
}

dense_qp_dims *dense_qp_dims_assign(void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_dims *dims = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    d_create_dense_qp_dim(dims, c_ptr);
    c_ptr += d_memsize_dense_qp_dim();

    assert((char *) raw_memory + dense_qp_dims_calculate_size() == c_ptr);

    return dims;
}

/************************************************
 * in
 ************************************************/

int dense_qp_in_calculate_size(void *config, dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_in);
    size += sizeof(dense_qp_dims);
    size += d_memsize_dense_qp(dims);
    return size;
}

dense_qp_in *dense_qp_in_assign(void *config, dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_in *qp_in = (dense_qp_in *) c_ptr;
    c_ptr += sizeof(dense_qp_in);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp(dims, qp_in, c_ptr);
    c_ptr += d_memsize_dense_qp(dims);

    qp_in->dim = (dense_qp_dims *) c_ptr;
    c_ptr += sizeof(dense_qp_dims);

    qp_in->dim->nv = dims->nv;
    qp_in->dim->ne = dims->ne;
    qp_in->dim->nb = dims->nb;
    qp_in->dim->ng = dims->ng;
    qp_in->dim->ns = dims->ns;

    assert((char *) raw_memory + dense_qp_in_calculate_size(config, dims) == c_ptr);

    return qp_in;
}

/************************************************
 * out
 ************************************************/

int dense_qp_out_calculate_size(void *config, dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_out);
    size += d_memsize_dense_qp_sol(dims);
    size += sizeof(dense_qp_info);
    return size;
}

dense_qp_out *dense_qp_out_assign(void *config, dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_out *qp_out = (dense_qp_out *) c_ptr;
    c_ptr += sizeof(dense_qp_out);

    assert((size_t) c_ptr % 8 == 0 && "memory not 8-byte aligned!");

    d_create_dense_qp_sol(dims, qp_out, c_ptr);
    c_ptr += d_memsize_dense_qp_sol(dims);

    qp_out->misc = (void *) c_ptr;
    c_ptr += sizeof(dense_qp_info);

    assert((char *) raw_memory + dense_qp_out_calculate_size(config, dims) == c_ptr);

    return qp_out;
}

/************************************************
 * res
 ************************************************/

int dense_qp_res_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_res);
    size += d_memsize_dense_qp_res(dims);
    return size;
}

dense_qp_res *dense_qp_res_assign(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_res *qp_res = (dense_qp_res *) c_ptr;
    c_ptr += sizeof(dense_qp_res);

    d_create_dense_qp_res(dims, qp_res, c_ptr);
    c_ptr += d_memsize_dense_qp_res(dims);

    assert((char *) raw_memory + dense_qp_res_calculate_size(dims) == c_ptr);

    return qp_res;
}

int dense_qp_res_workspace_calculate_size(dense_qp_dims *dims)
{
    int size = sizeof(dense_qp_res_ws);
    size += d_memsize_dense_qp_res_workspace(dims);
    return size;
}

dense_qp_res_ws *dense_qp_res_workspace_assign(dense_qp_dims *dims, void *raw_memory)
{
    char *c_ptr = (char *) raw_memory;

    dense_qp_res_ws *res_ws = (dense_qp_res_ws *) c_ptr;
    c_ptr += sizeof(dense_qp_res_ws);

    d_create_dense_qp_res_workspace(dims, res_ws, c_ptr);
    c_ptr += d_memsize_dense_qp_res_workspace(dims);

    assert((char *) raw_memory + dense_qp_res_workspace_calculate_size(dims) == c_ptr);

    return res_ws;
}

void dense_qp_res_compute(dense_qp_in *qp_in, dense_qp_out *qp_out, dense_qp_res *qp_res,
                          dense_qp_res_ws *res_ws)
{
    int nvd = qp_in->dim->nv;
    // int ned = qp_in->dim->ne;
    int ngd = qp_in->dim->ng;
    int nbd = qp_in->dim->nb;
    int nsd = qp_in->dim->ns;

    int *idxb = qp_in->idxb;
    int *idxs = qp_in->idxs;

    struct blasfeo_dvec *tmp_nbg = res_ws->tmp_nbg;

    // compute slacks for general constraints
    blasfeo_dgemv_t(nvd, ngd, 1.0, qp_in->Ct, 0, 0, qp_out->v, 0, -1.0, qp_in->d, nbd, qp_out->t,
                    nbd);
    blasfeo_dgemv_t(nvd, ngd, -1.0, qp_in->Ct, 0, 0, qp_out->v, 0, -1.0, qp_in->d, 2 * nbd + ngd,
                    qp_out->t, 2 * nbd + ngd);

    // compute slacks for bounds
    blasfeo_dvecex_sp(nbd, 1.0, idxb, qp_out->v, 0, tmp_nbg + 0, 0);
    blasfeo_daxpby(nbd, 1.0, tmp_nbg + 0, 0, -1.0, qp_in->d, 0, qp_out->t, 0);
    blasfeo_daxpby(nbd, -1.0, tmp_nbg + 0, 0, -1.0, qp_in->d, nbd + ngd, qp_out->t, nbd + ngd);

    // soft
    blasfeo_dvecad_sp(nsd, 1.0, qp_out->v, nvd, idxs, qp_out->t, 0);
    blasfeo_dvecad_sp(nsd, 1.0, qp_out->v, nvd+nsd, idxs, qp_out->t, nbd+ngd);
//    blasfeo_daxpy(2*nsd, -1.0, qp_out->v, nvd, qp_out->t, 2*nbd+2*ngd, qp_out->t, 2*nbd+2*ngd);
    blasfeo_dvecse(2*nsd, 0.0, qp_out->t, 2*nbd+2*ngd);

    // compute residuals
    d_compute_res_dense_qp(qp_in, qp_out, qp_res, res_ws);
}



void dense_qp_res_compute_nrm_inf(dense_qp_res *qp_res, double res[4])
{
    int nv = qp_res->dim->nv;
    int nb = qp_res->dim->nb;
    int ne = qp_res->dim->ne;
    int ng = qp_res->dim->ng;
    int ns = qp_res->dim->ns;

    blasfeo_dvecnrm_inf(nv + 2 * ns, qp_res->res_g, 0, &res[0]);
    blasfeo_dvecnrm_inf(ne, qp_res->res_b, 0, &res[1]);
    blasfeo_dvecnrm_inf(2 * nb + 2 * ng + 2 * ns, qp_res->res_d, 0, &res[2]);
    blasfeo_dvecnrm_inf(2 * nb + 2 * ng + 2 * ns, qp_res->res_m, 0, &res[3]);
}


