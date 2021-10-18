/*
    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "flint/fmpz.h"

void
arb_bin_ui(arb_t x, const arb_t n, ulong k, slong prec)
{
    if (k == 0)
    {
        arb_one(x);
    }
    else if (k == 1)
    {
        arb_set_round(x, n, prec);
    }
    else
    {
        arb_t t, u;

        arb_init(t);
        arb_init(u);

        arb_sub_ui(t, n, k - 1, prec);
        arb_rising_ui(t, t, k, prec);
        arb_fac_ui(u, k, prec);
        arb_div(x, t, u, prec);

        arb_clear(t);
        arb_clear(u);
    }
}

static
void
_arb_bin_uiui_small(arb_t x, ulong n, ulong k, slong prec)
{
    fmpz_t b;
    fmpz_init(b);
    fmpz_bin_uiui(b, n, k);
    arb_set_round_fmpz(x, b, prec);
    fmpz_clear(b);
}

void
arb_bin_uiui(arb_t x, ulong n, ulong k, slong prec)
{
    arb_t t;

    if (n < 3000)
        _arb_bin_uiui_small(x, n, k, prec);
    else
    {
        arb_init(t);
        arb_set_ui(t, n);
        arb_bin_ui(x, t, k, prec);
        arb_clear(t);
    }
}

