#include <iostream>
using namespace std;

#ifndef MonRed
#define MonRed

bool RelativelyPrime(int a, int b);
long long modInverse(long long A, long long M);

long long mon_red(long long a, long long b, long long n) {

    // PART 1: Choose r
    long long r = n+1;
    while(RelativelyPrime(r, n) == false) {
        r++;
    }

    //cout << "r " << r << endl << "n " << n << endl;

    //PART 2 : Compute k
    long long r_inv = modInverse(r, n);
    long long k  = (r_inv*r - 1)/n;

    //PART 3: Compute a, b 
    a = (a*r)%n;
    b = (b*r)%n;

    //PART 4: Compute x
    long long x = a*b;

    //PART 5: Compute s
    long long s = (x*k)%r;

    //PART 6, 7: Compute u`
    long long u = (x + s*n)/r;

    //PART 8: Compute c
    long long c = u<n ? u : (u-n);
    c = (c*r_inv)%n;

    return c;
}

bool RelativelyPrime(int a, int b)
{ // Assumes a, b > 0
    return (a < b) ? RelativelyPrime(b, a) : !(a % b) ? (b == 1)
                : RelativelyPrime(b, a % b);
}

long long modInverse(long long A, long long M)
{
    for (long long X = 1; X < M; X++)
        if (((A % M) * (X % M)) % M == 1)
            return X;

    return -1;
}

#endif