#include <iostream>
#include "Constant.h"
#include "Header.h"
#include <random>

using namespace std;

int main(int argc, char **argv) {

    // NOTE set seed to make it easier to check with python code
    mt19937 e2(0);

    uint32_t *a1 = new uint32_t[pad_N];
    uint32_t *a2 = new uint32_t[pad_N];

    //NOTE random and deepcopy
    for (int i = 0; i < N; i++)
    {
        a1[i] = e2() % q;
        a2[i] = e2() % q;
    }

    // NOTE transform a1, a2 to b1, b2
    uint32_t **b1 = new uint32_t*[Y_coef];
    uint32_t **b2 = new uint32_t*[Y_coef];
    for(int i=0; i<Z_coef; i++) {
        b1[i] = new uint32_t[Z_coef];
        b2[i] = new uint32_t[Z_coef];
    }

    transform_1d_to_2d(a1, b1);
    transform_1d_to_2d(a2, b2);

    // NOTE NTT in Z(512) axis

    for (int i = 0; i < Y_coef; i++)
    {
        radix2_ntt_mr(b1[i], Z_coef, false);
        radix2_ntt_mr(b2[i], Z_coef, false);
    }

    // PART NTT in Y(3) axis
    uint32_t **radix3_ntt = new uint32_t*[3];
    uint32_t **radix3_intt = new uint32_t*[3];
    for(int i=0; i<3; i++) {
        radix3_ntt[i] = new uint32_t[3];
        radix3_intt[i] = new uint32_t[3];
    }

    load_radix3_variable(1, radix3_ntt);
    load_radix3_variable(2, radix3_intt);

    //PART NTT in Y(3) axis
    
    radix3_mul_mont(radix3_ntt, b1, Z_coef);
    radix3_mul_mont(radix3_ntt, b2, Z_coef);

    //PART pointwise multiplication
    uint32_t **b3 = pointwise_mul_2d_mr(Y_coef, Z_coef, b1, b2);

    // PART INTT in Y(3) axis
    radix3_mul_mont(radix3_intt, b3, Z_coef);

    // PART INTT in Z(512) axis
    for (int i = 0; i < Y_coef; i++)
    {
        radix2_ntt_mr(b3[i], Z_coef, true);
    }

    // NOTE transform b3 to a3
    uint32_t *a3 = transform_2d_to_1d(b3);
    display_1d_mtx(a3, pad_N);    

    // NOTE delete all the pointer
    delete [] a1;
    delete [] a2;
    for(int i=0; i<Y_coef; i++) {
        delete [] b1[i];
        delete [] b2[i];
        delete [] b3[i];
        delete [] radix3_ntt[i];
        delete [] radix3_intt[i];
    }
    delete [] b1;
    delete [] b2;
    delete [] b3;
    delete [] a3;
    delete [] radix3_ntt;
    delete [] radix3_intt;

    cout << "pass" << endl;
    return 0;
}

