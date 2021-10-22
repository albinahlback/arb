/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_poly.h"
#include "flint/arith.h"
#include "acb.h"
#include "acb_poly.h"
#include "arb_fmpz_poly.h"

static
void
_lerch_transcendent_s_negative_integer(
        acb_t res,
        const acb_t z,
        slong s,
        const acb_t a,
        slong prec
        )
{
    /* Calculate the Lerch transcendent via
     *
     *      \phi(x, a, -n) = -B_{n + 1}(a, x) / (n + 1)
     *
     * where B_{n + 1}(a, x) is some sort of a generalized Bernoulli polynomial
     * found in Apostol's paper "On the Lerch zeta function". Here one finds
     * that this generalized Bernoulli polynomial coincides with
     *
     *      B_{n + 1}(a, x) = B_{n + 1}(a + x).
     */
    acb_t tmp;

    acb_init(tmp);

    acb_add(tmp, z, a, prec);
    acb_bernoulli_poly_ui(res, 1 - s, tmp, prec);
    acb_div_si(res, res, s - 1, prec);

    acb_clear(tmp);
}

static
void
_lerch_transcendent_directsum(
        acb_t res,
        const acb_t z,
        const acb_t s,
        const acb_t a,
        slong prec
        )
{
    slong ix, bound;
    arb_t logabsz, klowbnd, invabsz, aprec;
    acb_t zpow, negs, tmp;

    arb_init(logabsz);
    arb_init(klowbnd);
    arb_init(invabsz);
    arb_init(aprec);
    acb_init(zpow);
    acb_init(negs);
    acb_init(tmp);

    if (arb_is_negative(acb_realref(a)))
    {
        
    }
    else
    {
        bound = prec^
    }

    /* logabsz = log |z| */
    arb_mul(logabsz, acb_realref(z), acb_realref(z), prec);
    arb_addmul(logabsz, acb_imagref(z), acb_imagref(z), prec);
    arb_log(logabsz, logabsz, prec);
    arb_div_si(logabsz, logabsz, 2, prec);

    /* klowbnd = max(-Re(a), Re(s) / ln|z| - Re(a)) */

    arb_div(klowbnd, acb_realref(s), logabsz, prec);
    arb_sub(klowbnd, klowbnd, acb_realref(a), prec);
    arb_max(klowbnd, klowbnd, acb_realref(a), prec);

    acb_abs(invabsz, z, prec);
    arb_inv(invabsz, invabsz, prec);
    arb_set_si(aprec, prec);
    arb_pow(aprec, aprec, invabsz, prec);

    acb_neg(negs, s);

    /* First set res = 1 / a^s */
    acb_pow(res, a, negs, prec);

    /* Add z / (1 + a)^s to res */
    acb_add_si(tmp, a, 1, prec);
    acb_pow(tmp, tmp, negs, prec);
    acb_addmul(res, z, tmp, prec);

    /* Add z^ix / (ix + a)^s to res */
    acb_set(zpow, z);
    for (ix = 2; ix <= bound; ix++)
    {
        acb_mul(zpow, zpow, z, prec);
        acb_add_si(tmp, a, ix, prec);
        acb_pow(tmp, tmp, negs, prec);
        acb_addmul(res, zpow, tmp, prec);
    }

    arb_clear(lninvz);
    arb_clear(klowbnd);
    arb_clear(invabsz);
    arb_clear(aprec);
    acb_clear(zpow);
    acb_clear(negs);
    acb_clear(tmp);
}

void
lerch_transcendent(
        acb_t res,
        const acb_t z,
        const acb_t s,
        const acb_t a,
        slong prec
        )
{
    arb_t tmp;

    /* If z = 0, return 1 / a^s */
    if (acb_is_zero(z))
    {
        acb_neg(res, s);
        acb_pow(res, a, res, prec);
        return;
    }

    /* If z = 1, return \zeta(s, a) */
    if (acb_is_one(z))
    {
        acb_hurwitz_zeta(res, s, a, prec);
        return;
    }

    /* If s = 0, return \frac{1}{1 - z} */
    if (acb_is_zero(s))
    {
        acb_sub_si(res, z, 1, prec);
        acb_neg(res, res);
        acb_inv(res, res, prec);
        return;
    }

    /* If s is a negative integer */
    if (acb_is_int(s) && arb_is_negative(acb_realref(s)))
    {
        slong k = arf_get_si(arb_midref(acb_realref(s)), ARF_RND_NEAR);
        _lerch_transcendent_s_negative_integer(res, z, k, a, prec);
        return;
    }

    arb_init(tmp);

    arb_mul(tmp, acb_realref(z), acb_realref(z), prec);
    arb_addmul(tmp, acb_imagref(z), acb_imagref(z), prec);
    arb_mul_ui(tmp, tmp, 10, prec);
    arb_sub_ui(tmp, tmp, 7, prec);

    /* If |z| < 0.7 < 1 / \sqrt{2}, calculate via direct summation. The number
     * of iteration required for n correct bits should be around n^2
     * iterations. */
    if (arb_is_negative(tmp))
    {
        _lerch_transcendent_directsum(res, z, s, a, prec);
    }

    arb_clear(tmp);
}
