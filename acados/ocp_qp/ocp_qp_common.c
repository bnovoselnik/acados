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

#include "acados/ocp_qp/ocp_qp_common.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <assert.h>

#include "acados/ocp_qp/ocp_qp_condensing_hpipm.h"
#include "acados/ocp_qp/ocp_qp_condensing_qpoases.h"
#include "acados/ocp_qp/ocp_qp_hpipm.h"
#ifdef ACADOS_WITH_HPMPC
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#endif
#ifdef OOQP
#include "acados/ocp_qp/ocp_qp_ooqp.h"
#endif
#include "acados/ocp_qp/ocp_qp_qpdunes.h"
#include "acados/utils/types.h"
#include "hpipm/include/hpipm_d_cond.h"
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"

int_t ocp_qp_in_calculate_size(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                               const int_t *nc) {

    int_t bytes = sizeof(ocp_qp_in);

    bytes += 5*(N+1)*sizeof(int_t);  // nx, nu, nb, nc, ns
    bytes += 3*N*sizeof(real_t *);  // A, B, b
    bytes += 11*(N+1)*sizeof(real_t *);  // ...
    bytes += 1*(N+1)*sizeof(int_t *);  // idxb

    for (int_t k = 0; k < N+1; k++) {

        if (k < N) {
            bytes += nx[k+1]*nx[k]*sizeof(real_t);  // A
            bytes += nx[k+1]*nu[k]*sizeof(real_t);  // B
            bytes += nx[k+1]*sizeof(real_t);  // b
        }

        bytes += nx[k]*nx[k]*sizeof(real_t);  // Q
        bytes += nu[k]*nx[k]*sizeof(real_t);  // S
        bytes += nu[k]*nu[k]*sizeof(real_t);  // R
        bytes += nx[k]*sizeof(real_t);  // q
        bytes += nu[k]*sizeof(real_t);  // r
        bytes += nb[k]*sizeof(int_t);  // idxb
        bytes += 2*nb[k]*sizeof(real_t);  // lb, ub
        bytes += nc[k]*nx[k]*sizeof(real_t);  // Cx
        bytes += nc[k]*nu[k]*sizeof(real_t);  // Cu
        bytes += 2*nc[k]*sizeof(real_t);  // lc, uc
    }

    bytes = (bytes+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    bytes += ALIGNMENT;

    return bytes;
}


char *assign_ocp_qp_in(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                       const int_t *nc, ocp_qp_in **qp_in, void *ptr) {

    // pointer to initialize QP data to zero
    char *c_ptr_QPdata;

    // char pointer
    char *c_ptr = (char *) ptr;

    *qp_in = (ocp_qp_in *) c_ptr;
    c_ptr += sizeof(ocp_qp_in);

    // copy dimensions to workspace
    (*qp_in)->N = N;

    (*qp_in)->nx = (int_t *) c_ptr;
    memcpy(c_ptr, nx, (N+1)*sizeof(int_t));
    c_ptr += (N+1)*sizeof(int_t);

    (*qp_in)->nu = (int_t *) c_ptr;
    memcpy(c_ptr, nu, (N+1)*sizeof(int_t));
    c_ptr += (N+1)*sizeof(int_t);

    (*qp_in)->nb = (int_t *) c_ptr;
    memcpy(c_ptr, nb, (N+1)*sizeof(int_t));
    c_ptr += (N+1)*sizeof(int_t);

    (*qp_in)->nc = (int_t *) c_ptr;
    memcpy(c_ptr, nc, (N+1)*sizeof(int_t));
    c_ptr += (N+1)*sizeof(int_t);

    // not used atm, set to zero
    (*qp_in)->ns = (int_t *) c_ptr;
    memset(c_ptr, 0, (N+1)*sizeof(int_t));
    c_ptr += (N+1)*sizeof(int_t);

    // assign double pointers
    (*qp_in)->A = (const real_t **) c_ptr;
    c_ptr += N*sizeof(real_t *);

    (*qp_in)->B = (const real_t **) c_ptr;
    c_ptr += N*sizeof(real_t *);

    (*qp_in)->b = (const real_t **) c_ptr;
    c_ptr += N*sizeof(real_t *);

    (*qp_in)->Q = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->S = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->R = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->q = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->r = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->idxb = (const int_t **) c_ptr;
    c_ptr += (N+1)*sizeof(int_t *);

    (*qp_in)->lb = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->ub = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->Cx = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->Cu = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->lc = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_in)->uc = (const real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    // assign pointers to ints
    for (int_t k = 0; k < N+1; k++) {
        (*qp_in)->idxb[k] = (int_t *) c_ptr;
        c_ptr += nb[k]*sizeof(int_t);
    }

    // align data
    size_t l_ptr = (size_t) c_ptr;
    l_ptr = (l_ptr+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    c_ptr = (char *) l_ptr;

    // assign pointers to doubles
    c_ptr_QPdata = c_ptr;

    for (int_t k = 0; k < N+1; k++) {
        // printf("%zu MODULO %d = %zu\n", (size_t)c_ptr, 8, (size_t)c_ptr % 8);
        assert((size_t)c_ptr % 8 == 0);

        if (k < N) {
            (*qp_in)->A[k] = (real_t *) c_ptr;
            c_ptr += nx[k+1]*nx[k]*sizeof(real_t);

            (*qp_in)->B[k] = (real_t *) c_ptr;
            c_ptr += nx[k+1]*nu[k]*sizeof(real_t);

            (*qp_in)->b[k] = (real_t *) c_ptr;
            c_ptr += nx[k+1]*sizeof(real_t);
        }

        (*qp_in)->Q[k] = (real_t *) c_ptr;
        c_ptr += nx[k]*nx[k]*sizeof(real_t);

        (*qp_in)->S[k] = (real_t *) c_ptr;
        c_ptr += nu[k]*nx[k]*sizeof(real_t);

        (*qp_in)->R[k] = (real_t *) c_ptr;
        c_ptr += nu[k]*nu[k]*sizeof(real_t);

        (*qp_in)->q[k] = (real_t *) c_ptr;
        c_ptr += nx[k]*sizeof(real_t);

        (*qp_in)->r[k] = (real_t *) c_ptr;
        c_ptr += nu[k]*sizeof(real_t);

        (*qp_in)->lb[k] = (real_t *) c_ptr;
        c_ptr += nb[k]*sizeof(real_t);

        (*qp_in)->ub[k] = (real_t *) c_ptr;
        c_ptr += nb[k]*sizeof(real_t);

        (*qp_in)->Cx[k] = (real_t *) c_ptr;
        c_ptr += nc[k]*nx[k]*sizeof(real_t);

        (*qp_in)->Cu[k] = (real_t *) c_ptr;
        c_ptr += nc[k]*nu[k]*sizeof(real_t);

        (*qp_in)->lc[k] = (real_t *) c_ptr;
        c_ptr += nc[k]*sizeof(real_t);

        (*qp_in)->uc[k] = (real_t *) c_ptr;
        c_ptr += nc[k]*sizeof(real_t);
    }

    // set QP data to zero (mainly for valgrind)
    for (char *idx = c_ptr_QPdata; idx < c_ptr; idx++)
        *idx = 0;

    return c_ptr;
}


ocp_qp_in *create_ocp_qp_in(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                            const int_t *nc) {

    ocp_qp_in *qp_in;

    int_t bytes = ocp_qp_in_calculate_size(N, nx, nu, nb, nc);

    // TODO(dimitris): replace with acados_malloc to replace malloc at one place if not supported
    void *ptr = malloc(bytes);

    // // set a value for debugging
    // char *c_ptr = (char *) ptr;
    // for (int_t i = 0; i < bytes; i++) c_ptr[i] = 13;

    char *ptr_end = assign_ocp_qp_in(N, nx, nu, nb, nc, &qp_in, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    // for (int_t i = 0; i < bytes; i++) printf("%d - ", c_ptr[i]);
    // exit(1);

    return qp_in;
}


int_t ocp_qp_out_calculate_size(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                                const int_t *nc) {

    int_t bytes = sizeof(ocp_qp_out);

    bytes += 4*(N+1)*sizeof(real_t *);  // u, x, lam_b, lam_c
    bytes += N*sizeof(real_t *);  // pi

    for (int_t k = 0; k < N+1; k++) {
        bytes += (nx[k] + nu[k])*sizeof(real_t);  // u, x
        if (k < N)
            bytes += (nx[k+1])*sizeof(real_t);  // pi
        bytes += (nb[k] + nc[k])*sizeof(real_t);  // lam_b, lam_c
    }

    bytes = (bytes+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    bytes += ALIGNMENT;

    return bytes;
}


char *assign_ocp_qp_out(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                        const int_t *nc, ocp_qp_out **qp_out,
    void *ptr) {

    // char pointer
    char *c_ptr = (char *) ptr;

    *qp_out = (ocp_qp_out *) c_ptr;
    c_ptr += sizeof(ocp_qp_out);

    // assign double pointers
    (*qp_out)->x = (real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_out)->u = (real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_out)->pi = (real_t **) c_ptr;
    c_ptr += N*sizeof(real_t *);

    (*qp_out)->lam_b = (real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    (*qp_out)->lam_c = (real_t **) c_ptr;
    c_ptr += (N+1)*sizeof(real_t *);

    // align data
    size_t l_ptr = (size_t) c_ptr;
    l_ptr = (l_ptr+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    c_ptr = (char *) l_ptr;

    // NOTE(dimitris): splitted the loops below to be able to print primal/dual solution at once

    // assign pointers to QP solution
    for (int_t k = 0; k < N+1; k++) {
        assert((size_t)c_ptr % 8 == 0);

        (*qp_out)->x[k] = (real_t *) c_ptr;
        c_ptr += nx[k]*sizeof(real_t);

        (*qp_out)->u[k] = (real_t *) c_ptr;
        c_ptr += nu[k]*sizeof(real_t);
    }

    for (int_t k = 0; k < N; k++) {
        assert((size_t)c_ptr % 8 == 0);
        (*qp_out)->pi[k] = (real_t *) c_ptr;
        c_ptr += nx[k+1]*sizeof(real_t);
    }

    for (int_t k = 0; k < N+1; k++) {
        assert((size_t)c_ptr % 8 == 0);
        (*qp_out)->lam_b[k] = (real_t *) c_ptr;
        c_ptr += nb[k]*sizeof(real_t);

        (*qp_out)->lam_c[k] = (real_t *) c_ptr;
        c_ptr += nc[k]*sizeof(real_t);
    }
    return c_ptr;
}

ocp_qp_out *create_ocp_qp_out(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                              const int_t *nc) {

    ocp_qp_out *qp_out;

    int_t bytes = ocp_qp_out_calculate_size(N, nx, nu, nb, nc);
    void *ptr = malloc(bytes);
    char *ptr_end = assign_ocp_qp_out(N, nx, nu, nb, nc, &qp_out, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    return qp_out;
}

int_t ocp_qp_res_calculate_size(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                                const int_t *nc) {
    int_t bytes = sizeof(ocp_qp_res);

    bytes += 10 * (N + 1) * sizeof(real_t *);  // res_r, res_q,
                                               // res_d_lb, res_d_ub,
                                               // res_d_lg, res_d_ug,
                                               // res_m_lb, res_m_ub,
                                               // res_m_lg, res_m_ug
    bytes += 1 * N * sizeof(real_t *);  // res_b

    for (int_t k = 0; k < N + 1; k++) {
        bytes += (nx[k] + nu[k]) * sizeof(real_t);         // res_r, res_q
        if (k < N)
            bytes += (nx[k + 1]) * sizeof(real_t);  // res_b
        bytes += 4 * (nb[k] + nc[k]) * sizeof(real_t);  // res_d_lb, res_d_ub,
                                                        // res_d_lg, res_d_ug,
                                                        // res_m_lb, res_m_ub,
                                                        // res_m_lg, res_m_ug
    }

    bytes = (bytes + ALIGNMENT - 1) / ALIGNMENT * ALIGNMENT;
    bytes += ALIGNMENT;

    return bytes;
}

char *assign_ocp_qp_res(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                        const int_t *nc, ocp_qp_res **qp_res, void *ptr) {
    // char pointer
    char *c_ptr = (char *)ptr;

    *qp_res = (ocp_qp_res *)c_ptr;
    c_ptr += sizeof(ocp_qp_res);

    // assign double pointers
    (*qp_res)->res_r = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*qp_res)->res_q = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*qp_res)->res_b = (real_t **)c_ptr;
    c_ptr += N * sizeof(real_t *);

    (*qp_res)->res_d_lb = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*qp_res)->res_d_ub = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*qp_res)->res_d_lg = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*qp_res)->res_d_ug = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*qp_res)->res_m_lb = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*qp_res)->res_m_ub = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*qp_res)->res_m_lg = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    (*qp_res)->res_m_ug = (real_t **)c_ptr;
    c_ptr += (N + 1) * sizeof(real_t *);

    // align data
    size_t l_ptr = (size_t)c_ptr;
    l_ptr = (l_ptr + ALIGNMENT - 1) / ALIGNMENT * ALIGNMENT;
    c_ptr = (char *)l_ptr;

    // assign pointers to QP residuals
    for (int_t k = 0; k < N + 1; k++) {
        assert((size_t)c_ptr % 8 == 0);

        (*qp_res)->res_r[k] = (real_t *)c_ptr;
        c_ptr += nu[k] * sizeof(real_t);

        (*qp_res)->res_q[k] = (real_t *)c_ptr;
        c_ptr += nx[k] * sizeof(real_t);

        (*qp_res)->res_d_lb[k] = (real_t *)c_ptr;
        c_ptr += nb[k] * sizeof(real_t);

        (*qp_res)->res_d_ub[k] = (real_t *)c_ptr;
        c_ptr += nb[k] * sizeof(real_t);

        (*qp_res)->res_d_lg[k] = (real_t *)c_ptr;
        c_ptr += nc[k] * sizeof(real_t);

        (*qp_res)->res_d_ug[k] = (real_t *)c_ptr;
        c_ptr += nc[k] * sizeof(real_t);

        (*qp_res)->res_m_lb[k] = (real_t *)c_ptr;
        c_ptr += nb[k] * sizeof(real_t);

        (*qp_res)->res_m_ub[k] = (real_t *)c_ptr;
        c_ptr += nb[k] * sizeof(real_t);

        (*qp_res)->res_m_lg[k] = (real_t *)c_ptr;
        c_ptr += nc[k] * sizeof(real_t);

        (*qp_res)->res_m_ug[k] = (real_t *)c_ptr;
        c_ptr += nc[k] * sizeof(real_t);
    }

    for (int_t k = 0; k < N; k++) {
        assert((size_t)c_ptr % 8 == 0);

        (*qp_res)->res_b[k] = (real_t *)c_ptr;
        c_ptr += nx[k + 1] * sizeof(real_t);
    }

    return c_ptr;
}

ocp_qp_res *create_ocp_qp_res(const int_t N, const int_t *nx, const int_t *nu, const int_t *nb,
                              const int_t *nc) {
    ocp_qp_res *qp_res;

    int_t bytes = ocp_qp_res_calculate_size(N, nx, nu, nb, nc);
    void *ptr = malloc(bytes);
    char *ptr_end = assign_ocp_qp_res(N, nx, nu, nb, nc, &qp_res, ptr);
    assert((char *)ptr + bytes >= ptr_end);
    (void)ptr_end;

    return qp_res;
}

int_t ocp_qp_res_memory_calculate_size(const ocp_qp_in *qp_in) {
    int N = qp_in->N;
    int *nx = (int *)qp_in->nx;
    int *nu = (int *)qp_in->nu;
    int *nb = (int *)qp_in->nb;
    int *ng = (int *)qp_in->nc;
    int *ns = (int *)qp_in->ns;

    struct d_ocp_qp qp;
    qp.N = N;
    qp.nx = nx;
    qp.nu = nu;
    qp.nb = nb;
    qp.ng = ng;
    qp.ns = ns;

    int size = 0;

    size += sizeof(ocp_qp_res_memory);

    size += 1 * sizeof(struct d_ocp_qp);      // qp
    size += 1 * sizeof(struct d_ocp_qp_sol);  // qp_sol
    size += 1 * sizeof(struct d_ocp_qp_res);  // res_workspace

    size += d_memsize_ocp_qp(N, nx, nu, nb, ng, ns);
    size += d_memsize_ocp_qp_sol(N, nx, nu, nb, ng, ns);
    size += d_memsize_ocp_qp_res(&qp);
    size += 4 * (N + 1) * sizeof(double *);  // lam_lb lam_ub lam_lg lam_ug
    size += 1 * (N + 1) * sizeof(int *);     // hidxb_rev
    for (int_t ii = 0; ii <= N; ii++) {
        size += nb[ii] * sizeof(int);  // hidxb_rev
        size += 2 * (nb[ii] + ng[ii]) * sizeof(double);  // lam_lb lam_ub lam_lg lam_ug
    }

    size = (size + 63) / 64 * 64;  // make multipl of typical cache line size
    size += 1 * 64;                // align once to typical cache line size

    return size;
}

char *assign_ocp_qp_res_memory(const ocp_qp_in *qp_in, ocp_qp_res_memory **qp_res_mem, void *ptr) {
    // extract problem size
    int N = qp_in->N;
    int *nx = (int *)qp_in->nx;
    int *nu = (int *)qp_in->nu;
    int *nb = (int *)qp_in->nb;
    int *ng = (int *)qp_in->nc;
    int *ns = (int *)qp_in->ns;

    // char pointer
    char *c_ptr = (char *) ptr;

    *qp_res_mem = (ocp_qp_res_memory *)c_ptr;
    c_ptr += sizeof(ocp_qp_res_memory);

    //
    (*qp_res_mem)->qp = (struct d_ocp_qp *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp);
    //
    (*qp_res_mem)->qp_sol = (struct d_ocp_qp_sol *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp_sol);
    //
    (*qp_res_mem)->res_workspace = (struct d_ocp_qp_res *)c_ptr;
    c_ptr += 1 * sizeof(struct d_ocp_qp_res);
    //
    (*qp_res_mem)->hlam_lb = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*qp_res_mem)->hlam_ub = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*qp_res_mem)->hlam_lg = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*qp_res_mem)->hlam_ug = (double **)c_ptr;
    c_ptr += (N + 1) * sizeof(double *);
    //
    (*qp_res_mem)->hidxb_rev = (int **)c_ptr;
    c_ptr += (N + 1) * sizeof(int *);

    //
    struct d_ocp_qp *qp = (*qp_res_mem)->qp;
    //
    struct d_ocp_qp_sol *qp_sol = (*qp_res_mem)->qp_sol;
    //
    struct d_ocp_qp_res *res_workspace = (*qp_res_mem)->res_workspace;

    // align memory to typical cache line size
    size_t s_ptr = (size_t)c_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    c_ptr = (char *)s_ptr;

    // ocp qp structure
    d_create_ocp_qp(N, nx, nu, nb, ng, ns, qp, c_ptr);
    c_ptr += qp->memsize;
    // ocp qp sol structure
    d_create_ocp_qp_sol(N, nx, nu, nb, ng, ns, qp_sol, c_ptr);
    c_ptr += qp_sol->memsize;
    // ipm workspace structure
    d_create_ocp_qp_res(qp, res_workspace, c_ptr);
    c_ptr += res_workspace->memsize;

    // assign multipliers
    for (int_t ii = 0; ii <= N; ii++) {
        (*qp_res_mem)->hlam_lb[ii] = (double *)c_ptr;
        c_ptr += nb[ii] * sizeof(double);

        (*qp_res_mem)->hlam_lg[ii] = (double *)c_ptr;
        c_ptr += ng[ii] * sizeof(double);

        (*qp_res_mem)->hlam_ub[ii] = (double *)c_ptr;
        c_ptr += nb[ii] * sizeof(double);

        (*qp_res_mem)->hlam_ug[ii] = (double *)c_ptr;
        c_ptr += ng[ii] * sizeof(double);
    }

    //
    for (int_t ii = 0; ii <= N; ii++) {
        (*qp_res_mem)->hidxb_rev[ii] = (int *)c_ptr;
        c_ptr += nb[ii] * sizeof(int);
    }

    return c_ptr;
}

ocp_qp_res_memory *create_ocp_qp_res_memory(const ocp_qp_in *qp_in) {
    ocp_qp_res_memory *mem;
    int_t memory_size = ocp_qp_res_memory_calculate_size(qp_in);
    void *raw_memory = malloc(memory_size);
    char *ptr_end = assign_ocp_qp_res_memory(qp_in, &mem, raw_memory);
    assert((char *)raw_memory + memory_size >= ptr_end);
    (void)ptr_end;

    return mem;
}

void ocp_qp_calculate_res(const ocp_qp_in *qp_in, const ocp_qp_out *qp_out, ocp_qp_res *qp_res,
                          ocp_qp_res_memory *mem) {
    // loop indices
    int ii, jj, kk;

    // extract ocp problem size
    int N = qp_in->N;
    int *nx = (int *) qp_in->nx;
    int *nu = (int *) qp_in->nu;
    int *nb = (int *) qp_in->nb;
    int *ng = (int *) qp_in->nc;
    int *ns = (int *) qp_in->ns;

    // extract memory
    double **hlam_lb = mem->hlam_lb;
    double **hlam_ub = mem->hlam_ub;
    double **hlam_lg = mem->hlam_lg;
    double **hlam_ug = mem->hlam_ug;
    struct d_ocp_qp *qp = mem->qp;
    struct d_ocp_qp_sol *qp_sol = mem->qp_sol;
    struct d_ocp_qp_res *res_workspace = mem->res_workspace;
    int **hidxb_rev = (int **) mem->hidxb_rev;

    // extract input data
    double **hA = (double **)qp_in->A;
    double **hB = (double **)qp_in->B;
    double **hb = (double **)qp_in->b;
    double **hQ = (double **)qp_in->Q;
    double **hS = (double **)qp_in->S;
    double **hR = (double **)qp_in->R;
    double **hq = (double **)qp_in->q;
    double **hr = (double **)qp_in->r;
    double **hd_lb = (double **)qp_in->lb;
    double **hd_ub = (double **)qp_in->ub;
    double **hC = (double **)qp_in->Cx;
    double **hD = (double **)qp_in->Cu;
    double **hd_lg = (double **)qp_in->lc;
    double **hd_ug = (double **)qp_in->uc;
    int **hidxb = (int **)qp_in->idxb;

    // extract solution struct members
    double **hx = qp_out->x;
    double **hu = qp_out->u;
    double **hpi = qp_out->pi;
    double **hlam_b = qp_out->lam_b;
    double **hlam_c = qp_out->lam_c;

    // extract output struct members
    double **hres_r = qp_res->res_r;
    double **hres_q = qp_res->res_q;
    double **hres_b = qp_res->res_b;
    double **hres_d_lb = qp_res->res_d_lb;
    double **hres_d_ub = qp_res->res_d_ub;
    double **hres_d_lg = qp_res->res_d_lg;
    double **hres_d_ug = qp_res->res_d_ug;
    double **hres_m_lb = qp_res->res_m_lb;
    double **hres_m_ub = qp_res->res_m_ub;
    double **hres_m_lg = qp_res->res_m_lg;
    double **hres_m_ug = qp_res->res_m_ug;

    // expand multipliers
    for (kk = 0; kk <= N; kk++) {
        // expand multipliers for lb and ub
        for (ii = 0; ii < nb[kk]; ii++) {
            if (hlam_b[kk][ii] >= 0.0)
                hlam_lb[kk][ii] = hlam_b[kk][ii];
            else
                hlam_ub[kk][ii] = hlam_b[kk][ii];
        }

        // expand multipliers for lg and ug
        for (ii = 0; ii < ng[kk]; ii++) {
            if (hlam_c[kk][ii] >= 0.0)
                hlam_lg[kk][ii] = hlam_c[kk][ii];
            else
                hlam_ug[kk][ii] = hlam_c[kk][ii];
        }
    }

    // compute bounds indices in order [u; x]
    for (ii = 0; ii <= N; ii++) {
        for (jj = 0; jj < nb[ii]; jj++) {
            if (hidxb[ii][jj] < nx[ii]) {  // state constraint
                hidxb_rev[ii][jj] = hidxb[ii][jj]+nu[ii];
            } else {  // input constraint
                hidxb_rev[ii][jj] = hidxb[ii][jj]-nx[ii];
            }
        }
    }

    // convert to ocp qp structure
    d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb_rev, hd_lb, hd_ub,
                           hC, hD, hd_lg, hd_ug, NULL, NULL, NULL, NULL, NULL, qp);

    // convert to ocp qp sol structure
    d_cvt_colmaj_to_ocp_qp_sol(qp, hu, hx, NULL, NULL, hpi, hlam_lb, hlam_ub, hlam_lg, hlam_ug,
                               NULL, NULL, qp_sol);

    // compute residuals
    d_compute_res_ocp_qp(qp, qp_sol, res_workspace);

    // convert residuals to column major
    d_cvt_ocp_qp_res_to_colmaj(qp, res_workspace, hres_r, hres_q, NULL, NULL, hres_b, hres_d_lb,
                               hres_d_ub, hres_d_lg, hres_d_ug, NULL, NULL, hres_m_lb, hres_m_ub,
                               hres_m_lg, hres_m_ug, NULL, NULL);

    return;
}

void ocp_qp_in_copy_dynamics(const real_t *A, const real_t *B, const real_t *b,
                             ocp_qp_in *qp_in, int_t stage) {

    real_t **hA = (real_t **) qp_in->A;
    real_t **hB = (real_t **) qp_in->B;
    real_t **hb = (real_t **) qp_in->b;

    memcpy(hA[stage], A, qp_in->nx[stage+1]*qp_in->nx[stage]*sizeof(real_t));
    memcpy(hB[stage], B, qp_in->nx[stage+1]*qp_in->nu[stage]*sizeof(real_t));
    memcpy(hb[stage], b, qp_in->nx[stage+1]*sizeof(real_t));
}

void ocp_qp_in_copy_objective(const real_t *Q, const real_t *S, const real_t *R, const real_t *q,
                              const real_t *r, ocp_qp_in *qp_in, int_t stage) {

        real_t **hQ = (real_t **) qp_in->Q;
        real_t **hS = (real_t **) qp_in->S;
        real_t **hR = (real_t **) qp_in->R;
        real_t **hq = (real_t **) qp_in->q;
        real_t **hr = (real_t **) qp_in->r;

        memcpy(hQ[stage], Q, qp_in->nx[stage]*qp_in->nx[stage]*sizeof(real_t));
        memcpy(hS[stage], S, qp_in->nu[stage]*qp_in->nx[stage]*sizeof(real_t));
        memcpy(hR[stage], R, qp_in->nu[stage]*qp_in->nu[stage]*sizeof(real_t));
        memcpy(hq[stage], q, qp_in->nx[stage]*sizeof(real_t));
        memcpy(hr[stage], r, qp_in->nu[stage]*sizeof(real_t));
}

ocp_qp_solver *create_ocp_qp_solver(const ocp_qp_in *qp_in, const char *solver_name,
                                    void *solver_options) {
    ocp_qp_solver *qp_solver = (ocp_qp_solver *) malloc(sizeof(ocp_qp_solver));

    qp_solver->qp_in = (ocp_qp_in *) qp_in;
    qp_solver->qp_out = create_ocp_qp_out(qp_in->N, qp_in->nx, qp_in->nu, qp_in->nb, qp_in->nc);
    qp_solver->args = solver_options;

    if (!strcmp(solver_name, "qpdunes")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_qpdunes_create_arguments(QPDUNES_NONLINEAR_MPC);
        qp_solver->fun = &ocp_qp_qpdunes;
        qp_solver->initialize = &ocp_qp_qpdunes_initialize;
        qp_solver->destroy = &ocp_qp_qpdunes_destroy;
#ifdef OOQP
    } else if (!strcmp(solver_name, "ooqp")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_ooqp_create_arguments();
        qp_solver->fun = &ocp_qp_ooqp;
        qp_solver->initialize = &ocp_qp_ooqp_initialize;
        qp_solver->destroy = &ocp_qp_ooqp_destroy;
#endif
    } else if (!strcmp(solver_name, "condensing_qpoases")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_condensing_qpoases_create_arguments(qp_in);
        qp_solver->fun = &ocp_qp_condensing_qpoases;
        qp_solver->initialize = &ocp_qp_condensing_qpoases_initialize;
        qp_solver->destroy = &ocp_qp_condensing_qpoases_destroy;
#ifdef ACADOS_WITH_HPMPC
    } else if (!strcmp(solver_name, "hpmpc")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_hpmpc_create_arguments(qp_in, HPMPC_DEFAULT_ARGUMENTS);
        qp_solver->fun = &ocp_qp_hpmpc;
        qp_solver->initialize = &ocp_qp_hpmpc_initialize;
        qp_solver->destroy = &ocp_qp_hpmpc_destroy;
#endif
    } else if (!strcmp(solver_name, "condensing_hpipm")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_condensing_hpipm_create_arguments(qp_in);
        qp_solver->fun = &ocp_qp_condensing_hpipm;
        qp_solver->initialize = &ocp_qp_condensing_hpipm_initialize;
        qp_solver->destroy = &ocp_qp_condensing_hpipm_destroy;
    } else if (!strcmp(solver_name, "hpipm")) {
        if (qp_solver->args == NULL)
            qp_solver->args = ocp_qp_hpipm_create_arguments(qp_in);
        qp_solver->fun = &ocp_qp_hpipm;
        qp_solver->initialize = &ocp_qp_hpipm_initialize;
        qp_solver->destroy = &ocp_qp_hpipm_destroy;
    } else {
        printf("Chosen QP solver not available\n");
        exit(1);
    }
    qp_solver->initialize(qp_solver->qp_in, qp_solver->args, &qp_solver->mem, &qp_solver->work);

    return qp_solver;
}
