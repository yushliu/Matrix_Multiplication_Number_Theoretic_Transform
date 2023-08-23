#include <iostream>
#include "Header.h"
#include "Constant.h"
#include <bitset>
#include <string>
#include <random>
using namespace std;

int main()
{
    //NOTE set seed to make it easier to check with python code
    mt19937 e2(0);
    
    //random a1, a2
    Mtx1d *a1 = new Mtx1d(pad_N);
    Mtx1d *a2 = new Mtx1d(pad_N);
    Mtx1d *school_book_a1 = new Mtx1d(pad_N);
    Mtx1d *school_book_a2 = new Mtx1d(pad_N);

    //NOTE random and deepcopy
    for (int i = 0; i < N; i++) {
        a1->data[i] = e2()%q;
        a2->data[i] = e2()%q;
    }

    for(int i=0; i<pad_N; i++) {
        school_book_a1->data[i] = a1->data[i];
        school_book_a2->data[i] = a2->data[i];
    }

    //NOTE transform a1, a2 to b1, b2
    Mtx2d *b1 = new Mtx2d(Y_coef, Z_coef);
    Mtx2d *b2 = new Mtx2d(Y_coef, Z_coef);

    transform_1d_to_2d(a1, b1);
    transform_1d_to_2d(a2, b2);

    //NOTE same pointer mtx1 and mtx2
    //cause 6 "Mtx1d constructor error"
    Mtx1d *b1_1d_ptr = new Mtx1d[Y_coef];
    Mtx1d *b2_1d_ptr = new Mtx1d[Y_coef];

    for(int i=0; i<Y_coef; i++) {
        b1_1d_ptr[i].set_same_ptr(Z_coef, b1->data[i]);
        b2_1d_ptr[i].set_same_ptr(Z_coef, b2->data[i]);
    }

    //NOTE NTT in Z(512) axis

    for(int i=0; i<Y_coef; i++) {
        recursive_radix2_ntt(&b1_1d_ptr[i], false);
        recursive_radix2_ntt(&b2_1d_ptr[i], false);
    }

    //PART NTT in Y(3) axis
    Mtx2d *radix3_ntt = new Mtx2d(3, 3);//don't use Y_coef, since it should always be 3
    Mtx2d *radix3_intt = new Mtx2d(3, 3);
    
    load_radix3_variable(1, radix3_ntt);
    load_radix3_variable(2, radix3_intt);

    //NOTE how to do it column wise
    for(int i=0; i<Z_coef; i++) {
        Mtx1d *temp_b1 = mul_2d_with_3_points(radix3_ntt,
        b1->data[0][i], b1->data[1][i], b1->data[2][i]);
        Mtx1d *temp_b2 = mul_2d_with_3_points(radix3_ntt,
        b2->data[0][i], b2->data[1][i], b2->data[2][i]);

        b1->data[0][i] = temp_b1->data[0];
        b1->data[1][i] = temp_b1->data[1];
        b1->data[2][i] = temp_b1->data[2];
        b2->data[0][i] = temp_b2->data[0];
        b2->data[1][i] = temp_b2->data[1];
        b2->data[2][i] = temp_b2->data[2];

        delete temp_b1;
        delete temp_b2;
    }

    //PART pointwise multiplication
    //works fine

    Mtx2d *b3 = pointwise_mul_2d(Y_coef, Z_coef, b1, b2);

    //PART inverse NTT in Y(3) axis
    //works fine

    for (int i = 0; i < Z_coef; i++) {
        Mtx1d *temp_b3 = mul_2d_with_3_points(radix3_intt,
        b3->data[0][i], b3->data[1][i], b3->data[2][i]);

        b3->data[0][i] = temp_b3->data[0];
        b3->data[1][i] = temp_b3->data[1];
        b3->data[2][i] = temp_b3->data[2];

        delete temp_b3;
    }

    //NOTE inverse NTT in Z(512) axis

    Mtx1d *b3_1d_ptr = new Mtx1d[Y_coef];

    for (int i = 0; i < Y_coef; i++) {
        b3_1d_ptr[i].set_same_ptr(Z_coef, b3->data[i]);
    }

    for (int i = 0; i < Y_coef; i++) {
        recursive_radix2_ntt(&b3_1d_ptr[i], true);
    }

    //NOTE transform b3 to a3
    for(int i=0; i<Y_coef; i++) {
        mod_Mtx1d_with_value(&b3_1d_ptr[i], q);
    }

    Mtx1d *a3 = transform_2d_to_1d(b3);
    cout << a3 << endl;

    delete a1;
    delete a2;
    delete school_book_a1;
    delete school_book_a2;
    delete b1;
    delete b2;
    delete[] b1_1d_ptr;
    delete[] b2_1d_ptr;
    delete radix3_ntt;
    delete radix3_intt;
    delete b3;
    delete[] b3_1d_ptr;
    delete a3;

    return 0;
}