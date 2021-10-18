/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_poly.h"
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
    /* Calculate
     *
     *      \sum_{n = 0}^\infty z^n / (n + a)^s
     *  =
     *      (1 - z)^{s - 1} p_k(z)
     */

    acb_t tmp;
    fmpz_poly_t poly;

    acb_init(tmp);

    /* tmp = (1 - z)^{-s - 1} */
    acb_sub_ui(tmp, z, 1, prec);
    acb_neg(tmp, tmp);
    acb_pow_si(tmp, tmp, s - 1, prec);

    /* If a is one, use Eulerian polynomial via
     *      z A_n(z) (1 - z)^{-s - 1}.          */
    if (acb_is_one(a))
    {
        fmpz_poly_eulerian_polynomial(poly, -s);
        arb_fmpz_poly_evaluate_acb(res, poly, z, prec);
        acb_mul(res, res, tmp, prec);
        acb_mul(res, res, z, prec);
    }

    /* Else use the generalized Bernoulli polynomial B_n(a, x) via the identity
     *      \phi(x, a, -n) = - B_{n + 1}(a, x) / (n + 1)
     * found in Apostol's paper "On the Lerch zeta function". */

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
}
