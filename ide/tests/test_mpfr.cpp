//
// Created by thallock on 10/24/21.
//


#include "catch2/catch.hpp"

#include <mpfr.h>
#include <iostream>

#include <random>


#define PRECISION 512


TEST_CASE("test mpfr", "[mpfr]") {
    mpfr_t s, t, u;

    mpfr_init2(s, PRECISION);
    mpfr_set_d(s, 3.14, MPFR_RNDD);
    mpfr_init2(t, PRECISION);
    mpfr_set_d(t, 2.4, MPFR_RNDD);
    mpfr_init2(u, PRECISION);

    mpfr_div(u, s, t, MPFR_RNDD);

    fprintf(stdout, "\ns=");
    mpfr_out_str(stdout, 10, 0, s, MPFR_RNDD);
    fprintf(stdout, "\nt=");
    mpfr_out_str(stdout, 10, 0, t, MPFR_RNDD); fflush(stdout);
    fprintf(stdout, "\nu=");
    mpfr_out_str(stdout, 10, 0, u, MPFR_RNDD); fflush(stdout);
    fprintf(stdout, "\n");
}
