/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#ifndef FMPRB_H
#define FMPRB_H

#include "fmpr.h"

typedef struct
{
    fmpr_struct mid;
    fmpr_struct rad;
}
fmprb_struct;

typedef fmprb_struct fmprb_t[1];

#define fmprb_midref(x) (&(x)->mid)
#define fmprb_radref(x) (&(x)->rad)


#define FMPRB_RAD_PREC 30
#define FMPRB_MIN_ADJUSTED_PREC 

static __inline__ void
fmprb_init(fmprb_t x)
{
    fmpr_init(fmprb_midref(x));
    fmpr_init(fmprb_radref(x));
}

static __inline__ void
fmprb_clear(fmprb_t x)
{
    fmpr_clear(fmprb_midref(x));
    fmpr_clear(fmprb_radref(x));
}

static __inline__ fmprb_struct *
_fmprb_vec_init(long n)
{
    long i;
    fmprb_struct * v = flint_malloc(sizeof(fmprb_struct) * n);

    for (i = 0; i < n; i++)
        fmprb_init(v + i);

    return v;
}

static __inline__ void
_fmprb_vec_clear(fmprb_struct * v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        fmprb_clear(v + i);
    flint_free(v);
}

static __inline__ int
fmprb_is_exact(const fmprb_t x)
{
    return fmpr_is_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_zero(fmprb_t x)
{
    fmpr_zero(fmprb_midref(x));
    fmpr_zero(fmprb_radref(x));
}

static __inline__ int
fmprb_is_zero(const fmprb_t x)
{
    return fmpr_is_zero(fmprb_midref(x)) && fmpr_is_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_set(fmprb_t x, const fmprb_t y)
{
    fmpr_set(fmprb_midref(x), fmprb_midref(y));
    fmpr_set(fmprb_radref(x), fmprb_radref(y));
}

void fmprb_set_round(fmprb_t z, const fmprb_t x, long prec);

static __inline__ void
fmprb_neg(fmprb_t x, const fmprb_t y)
{
    fmpr_neg(fmprb_midref(x), fmprb_midref(y));
    fmpr_set(fmprb_radref(x), fmprb_radref(y));
}

static __inline__ void
fmprb_abs(fmprb_t x, const fmprb_t y)
{
    fmpr_abs(fmprb_midref(x), fmprb_midref(y));
    fmpr_set(fmprb_radref(x), fmprb_radref(y));
}


static __inline__ void
fmprb_set_si(fmprb_t x, long y)
{
    fmpr_set_si(fmprb_midref(x), y);
    fmpr_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_set_ui(fmprb_t x, ulong y)
{
    fmpr_set_ui(fmprb_midref(x), y);
    fmpr_zero(fmprb_radref(x));
}

static __inline__ void
fmprb_set_fmpz(fmprb_t x, const fmpz_t y)
{
    fmpr_set_fmpz(fmprb_midref(x), y);
    fmpr_zero(fmprb_radref(x));
}

static __inline__ int
fmprb_is_one(const fmprb_t f)
{
    return fmpr_is_one(fmprb_midref(f)) && fmpr_is_zero(fmprb_radref(f));
}

static __inline__ void
fmprb_one(fmprb_t f)
{
    fmprb_set_ui(f, 1UL);
}

static __inline__ void
fmprb_adjust(fmprb_t x)
{
    if (!fmprb_is_exact(x))
    {
        /* reduce precision here */
    }
}

void fmprb_add(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_add_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_add_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_add_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);

void fmprb_addmul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_addmul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_addmul_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_addmul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);

void fmprb_div(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_div_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_div_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_div_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);
void fmprb_fmpz_div_fmpz(fmprb_t y, const fmpz_t num, const fmpz_t den, long prec);
void fmprb_ui_div(fmprb_t z, ulong x, const fmprb_t y, long prec);

void fmprb_div_2expm1_ui(fmprb_t y, const fmprb_t x, ulong n, long prec);

void fmprb_mul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_mul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_mul_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_mul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);

void fmprb_sqrt(fmprb_t z, const fmprb_t x, long prec);
void fmprb_sqrt_ui(fmprb_t z, ulong x, long prec);
void fmprb_sqrt_fmpz(fmprb_t z, const fmpz_t x, long prec);

void fmprb_sub(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_sub_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_sub_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_sub_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);

void fmprb_submul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec);
void fmprb_submul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec);
void fmprb_submul_si(fmprb_t z, const fmprb_t x, long y, long prec);
void fmprb_submul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec);

void fmprb_pow_ui(fmprb_t y, const fmprb_t b, ulong e, long prec);
void fmprb_ui_pow_ui(fmprb_t y, ulong b, ulong e, long prec);
void fmprb_si_pow_ui(fmprb_t y, long b, ulong e, long prec);

void fmprb_log(fmprb_t z, const fmprb_t x, long prec);
void fmprb_log_ui(fmprb_t z, ulong x, long prec);
void fmprb_log_fmpz(fmprb_t z, const fmpz_t x, long prec);

void fmprb_exp(fmprb_t z, const fmprb_t x, long prec);

void fmprb_fac_ui(fmprb_t x, ulong n, long prec);

void fmprb_bin_ui(fmprb_t x, const fmprb_t n, ulong k, long prec);
void fmprb_bin_uiui(fmprb_t x, ulong n, ulong k, long prec);

void fmprb_fib_fmpz(fmprb_t f, const fmpz_t n, long prec);
void fmprb_fib_ui(fmprb_t f, ulong n, long prec);

void fmprb_const_pi_chudnovsky(fmprb_t x, long prec);
void fmprb_const_pi(fmprb_t x, long prec);

void fmprb_const_euler_brent_mcmillan(fmprb_t res, long prec);
void fmprb_const_zeta3_bsplit(fmprb_t x, long prec);

void fmprb_zeta_ui_bsplit(fmprb_t x, ulong s, long prec);
void fmprb_zeta_ui_euler_product(fmprb_t z, ulong s, long prec);
void fmprb_zeta_ui_bernoulli(fmprb_t x, ulong n, long prec);
void fmprb_zeta_ui_vec_borwein(fmprb_struct * z, ulong start, long num, ulong step, long prec);
void fmprb_zeta_ui(fmprb_t x, ulong n, long prec);

void fmprb_zeta_ui_vec_even(fmprb_struct * x, ulong start, long num, long prec);
void fmprb_zeta_ui_vec_odd(fmprb_struct * x, ulong start, long num, long prec);
void fmprb_zeta_ui_vec(fmprb_struct * x, ulong start, long num, long prec);

static __inline__ void
fmprb_add_error_2exp_si(fmprb_t x, long err)
{
    fmpr_t t;
    fmpr_init(t);
    fmpz_set_ui(fmpr_manref(t), 1);
    fmpz_set_si(fmpr_expref(t), err);
    fmpr_add(fmprb_radref(x), fmprb_radref(x), t, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_clear(t);
}

static __inline__ void
fmprb_printd(const fmprb_t x, long digits)
{
    fmpr_printd(fmprb_midref(x), digits);
    printf(" +/- ");
    fmpr_printd(fmprb_radref(x), 5);
}

static __inline__ void
fmprb_mul_2exp_si(fmprb_t y, const fmprb_t x, long e)
{
    fmpr_mul_2exp_si(fmprb_midref(y), fmprb_midref(x), e);
    fmpr_mul_2exp_si(fmprb_radref(y), fmprb_radref(x), e);
}

static __inline__ void
fmprb_set_fmpq(fmprb_t y, const fmpq_t x, long prec)
{
    fmprb_fmpz_div_fmpz(y, fmpq_numref(x), fmpq_denref(x), prec);
}

int fmprb_contains_fmpr(const fmprb_t x, const fmpr_t y);
int fmprb_contains_fmpq(const fmprb_t x, const fmpq_t y);
int fmprb_contains_fmpz(const fmprb_t x, const fmpz_t y);
int fmprb_contains_mpfr(const fmprb_t x, const mpfr_t y);
int fmprb_contains_zero(const fmprb_t x);

static __inline__ long
fmprb_rel_error_bits(const fmprb_t x)
{
    fmpz_t midmag, radmag;
    long result;

    if (fmpr_is_zero(fmprb_radref(x)))
        return -FMPR_PREC_EXACT;
    if (fmpr_is_special(fmprb_midref(x)) || fmpr_is_special(fmprb_radref(x)))
        return FMPR_PREC_EXACT;

    fmpz_init(midmag);
    fmpz_init(radmag);

    fmpz_add_ui(midmag, fmpr_expref(fmprb_midref(x)),
        fmpz_bits(fmpr_manref(fmprb_midref(x))));
    fmpz_add_ui(radmag, fmpr_expref(fmprb_radref(x)),
        fmpz_bits(fmpr_manref(fmprb_radref(x))));
    fmpz_add_ui(radmag, radmag, 1);

    result = _fmpz_sub_small(radmag, midmag);

    fmpz_clear(midmag);
    fmpz_clear(radmag);

    return result;
}

static __inline__ long
fmprb_rel_accuracy_bits(const fmprb_t x)
{
    return -fmprb_rel_error_bits(x);
}

#endif

