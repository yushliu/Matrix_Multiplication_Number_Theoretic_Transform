#ifndef HEADER_H
#define HEADER_H

#include "Constant.h"
#include <iostream>
#include <bitset>
#include "Montgomery_reduction.h"
using namespace std;

class IntegerModRing {
private:
    long long range;
    long long value;

public:

IntegerModRing(){
    this->range = prime_q;
    this->value = 0;    
};

IntegerModRing(long long value){
    this->range = prime_q;
    this->value = value;
};

long long getRange(){
    return this->range;
};

long long getValue(){
    return this->value;
};

void setValue(long long value){
    this->value = value%this->range;
};

void setValue(IntegerModRing other){
    this->value = other.value;
    if(this->range != other.range) {
        cout << "range not equal when setValue" << endl;
    }
};

IntegerModRing operator+(IntegerModRing other){
    IntegerModRing result;
    result.value = (this->value + other.value) % this->range;
    return result;
};

IntegerModRing operator-(IntegerModRing other){
    IntegerModRing result;
    result.value = (this->value - other.value) % this->range;
    
    if(result.value < 0) {
        result.value += this->range;
    }
    
    return result;
};

IntegerModRing operator*(IntegerModRing other){
    IntegerModRing result;

    result.value = mon_red(this->value, other.value, this->range);
    return result;
};

IntegerModRing operator=(IntegerModRing other){
    this->value = other.value % this->range;
    if(this->range != other.range) {
        cout << "range not equal when = " << endl;
    }
    return *this;
};

friend ostream &operator<<(ostream &os, const IntegerModRing &obj) {
    os << obj.value;
    return os;
}

};

class Mtx1d {
private:
int length;
bool if_activate_destructor;

public:

IntegerModRing *data;

Mtx1d() {
    this->length = 0;
    this->data = new IntegerModRing[1];
    this->if_activate_destructor = false;
};

Mtx1d(int length) {
    this->length = length;
    this->data = new IntegerModRing[length];
    this->if_activate_destructor = false;
}

void set_same_ptr(int length, IntegerModRing *ptr) {
    this->length = length;
    delete[] this->data;
    this->data = ptr;
    this->if_activate_destructor = true;
}

int get_length() {
    return this->length;
}

~Mtx1d() {
    if(this->if_activate_destructor == false) {
        delete[] this->data; 
    }
}

void set_value(int index, long long value) {
    this->data[index].setValue(value);
}

//operator
Mtx1d operator=(Mtx1d other) {
    
    if(this->length != other.length) {
        cout << "Mtx1d operator= error, length not the same" << endl;
        return *this;
    }

    //NOTE deepcopy
    for(int i=0; i<this->length; i++) {
        this->data[i] = other.data[i];
    }

    return *this;
}

friend ostream &operator<<(ostream &os, const Mtx1d *obj) {
    os << "Mtx1d (" << obj->length << ")" << endl;
    for(int i = 0; i < obj->length; i++) {
        os << obj->data[i] << " ";
    }
    os << endl;
    return os;
}

};

class Mtx2d {
private:

int row;
int col;

public:

IntegerModRing **data;

Mtx2d()
{
    cout << "Mtx2d constructor error" << endl;
    this->row = 0;
    this->col = 0;
    this->data = NULL;
};

Mtx2d(int row, int column) {
    this->row = row;
    this->col = column;
    this->data = new IntegerModRing*[row];
    for(int i = 0; i < row; i++) {
        this->data[i] = new IntegerModRing[column];
    }    
}

~Mtx2d() {
    for(int i = 0; i < this->row; i++) {
        delete[] this->data[i];
    }
    delete[] this->data;
}

void set_value(int row, int column, long long value) {
    this->data[row][column].setValue(value);
}

int get_row() {
    return this->row;
}

int get_col() {
    return this->col;
}

friend ostream &operator<<(ostream &os, const Mtx2d *obj) {
    os << "Mtx2d(" << obj->row << ", " << obj->col << ")" << endl;
    for(int i = 0; i < obj->row; i++) {
        for(int j = 0; j < obj->col; j++) {
            os << obj->data[i][j] << " ";
        }
        os << endl;
    }
    return os;
}
};

Mtx1d* mul_2d_with_1d(Mtx2d *A2d, Mtx1d *A1d) {

    if (A2d->get_col() != A1d->get_length()) {
        cout << "Mtx2d * Mtx1d error" << endl;
        return new Mtx1d();
    }

    Mtx1d *result = new Mtx1d(A2d->get_row());

    for(int i=0; i<result->get_length(); i++) {

        IntegerModRing sum = IntegerModRing(0);
        
        for(int j=0; j<A2d->get_col(); j++) {
            sum = sum + (A2d->data[i][j] * A1d->data[j]);
        }

        result->data[i].setValue(sum);
    }

    return result;
}

Mtx1d* mul_2d_with_3_points(Mtx2d *A2d, IntegerModRing p0, IntegerModRing p1, IntegerModRing p2) {
    
    int coef = 3;

    if(A2d->get_col() != coef) {
        cout << "Mtx2d * mul_2d_with_3_points error with column not equals to 3" << endl;
        return new Mtx1d();
    }

    Mtx1d *result = new Mtx1d(coef);

    for(int i=0; i<coef; i++) {
        
        IntegerModRing sum = IntegerModRing(0);

        sum = sum + (A2d->data[i][0] * p0);
        sum = sum + (A2d->data[i][1] * p1);
        sum = sum + (A2d->data[i][2] * p2);

        result->data[i].setValue(sum);
    }

    return result;
}

//transform from a1 to b1
void transform_1d_to_2d(Mtx1d *a, Mtx2d *b) {
    for(int i=0; i<Y_coef*Z_coef; i++) {
        b->data[i%Y_coef][i%Z_coef].setValue(a->data[i]);
    }
}

Mtx1d* transform_2d_to_1d(Mtx2d *a) {

    if (Y_coef != 3 || Z_coef != 512) {
        cout << "transform_1d_to_2d does not support current Y_coef and Z_coef" << endl;
    }

    Mtx1d *result = new Mtx1d(Y_coef*Z_coef);
    
    for(int i=0; i<Y_coef; i++) {
        for(int j=0; j<Z_coef; j++) {
            result->data[(1024*i+513*j)%pad_N] = a->data[i][j];
        }
    }

    return result;
}

Mtx1d* pointwise_mul_1d(Mtx1d *A, Mtx1d *B) {
    
    if(A->get_length() != B->get_length()) {
        cout << "pointwise_mul_1d error with length not match" << endl;
        return new Mtx1d();
    }

    int length = A->get_length();
    Mtx1d *result = new Mtx1d(length);

    for(int i=0; i<length; i++) {
        result->data[i] = A->data[i] * B->data[i];
    }

    return result;
}

Mtx2d* pointwise_mul_2d(int row, int col, Mtx2d *A, Mtx2d *B) {
    Mtx2d *result = new Mtx2d(row, col);

    for(int i=0; i<row; i++) {
        for(int j=0; j<col; j++) {
            result->data[i][j] = A->data[i][j] * B->data[i][j];
        }
    }

    return result;
}

// choice 1: NTT
// choice 2: INTT
void load_radix3_variable(int choice, Mtx2d* a) {


    if(a->get_col()!=3 || a->get_row()!=3) {
        cout << "radix col or row error" << endl;
        return;
    }

    if(choice==1) {

        a->data[0][0].setValue(1);
        a->data[0][1].setValue(1025640636);
        a->data[0][2].setValue(394137924);

        a->data[1][0].setValue(1);
        a->data[1][1].setValue(394137924);
        a->data[1][2].setValue(1025640636);

        a->data[2][0].setValue(1);
        a->data[2][1].setValue(1);
        a->data[2][2].setValue(1);

    } else if(choice==2) {

        a->data[0][0].setValue(946519041);
        a->data[0][1].setValue(946519041);
        a->data[0][2].setValue(946519041);

        a->data[1][0].setValue(131379308);
        a->data[1][1].setValue(341880212);
        a->data[1][2].setValue(946519041);

        a->data[2][0].setValue(341880212);
        a->data[2][1].setValue(131379308);
        a->data[2][2].setValue(946519041);

    } else {
        cout << "choice not found" << endl;
    }
}

//TODO check if it's prime
//TODO did not check if len(seq) is 2^k, k is int
//inverse = false means NTT
//inverse = true means INTT
void recursive_radix2_ntt(Mtx1d *seq, bool inverse) {
    
    //bit reverse the seq
    for(unsigned int i=1; i<seq->get_length(); i++) {
        string binary_j = std::bitset<power_of_2>(i).to_string();
        reverse(binary_j.begin(), binary_j.end());
        unsigned int j = stoi(binary_j, 0, 2);
        if(i < j) {
            IntegerModRing temp = seq->data[i];
            seq->data[i] = seq->data[j];
            seq->data[j] = temp;
        }
    }

    IntegerModRing root(root_test);
    if(inverse) {
        root = rt_inverse;
    }

    Mtx1d *roots_seq = new Mtx1d(seq->get_length()/2);
    roots_seq->data[0].setValue(1);
    for(int i=1; i<roots_seq->get_length(); i++) {
        roots_seq->data[i] = roots_seq->data[i-1] * root;
    }

    //Main radix2 starts here
    int h = 2;
    while(h <= seq->get_length()) {
        for(int i=0; i<seq->get_length(); i+=h) {
            for(int j=0; j<h/2; j++) {

                IntegerModRing u = seq->data[i+j];
                IntegerModRing v = seq->data[i+j+h/2] * roots_seq->data[(seq->get_length()/h)*j];

                seq->data[i+j] = u + v;
                seq->data[i+j+h/2] = u - v;
            }
        }
        h *= 2;
    }

    if(inverse) {
        IntegerModRing rv(rv_inverse);

        for(int i=0; i<seq->get_length(); i++) {
            seq->data[i] = seq->data[i] * rv;
        }
    }
}

void mod_Mtx1d_with_value(Mtx1d *A, int value) {
    for(int i=0; i<A->get_length(); i++) {
        A->data[i].setValue(A->data[i].getValue()%value);
    }
}

#endif