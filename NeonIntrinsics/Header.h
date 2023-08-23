#include <iostream>
#include "Constant.h"
using namespace std;

#ifndef Header_no_class_h
#define Header_no_class_h

const uint32_t Mprime = 2248402433;
const uint32_t r2_mod_m_B = 422821617;
const uint32_t r2_mod_m_Bprime = 4187178225;

void display_1d_mtx(uint32_t *a, int len)
{
    cout << "display 1d mtx:" << endl;
    for (int i = 0; i < len; i++)
    {
        cout << a[i] << " ";
    }
    cout << endl;
}

void display_2d_mtx(uint32_t **a, int len1, int len2)
{
    cout << "display 2d mtx:" << endl;
    for (int i = 0; i < len1; i++)
    {
        for (int j = 0; j < len2; j++)
        {
            cout << a[i][j] << " ";
        }
        cout << endl << endl;
    }
}

// transform from a1 to b1
void transform_1d_to_2d(uint32_t *a, uint32_t **b)
{
    if(Z_coef != 512) {
        cout << "Z coef not equals to 512" << endl;
    }
    for(int i=0; i<Y_coef*Z_coef; i++) {
        b[i%Y_coef][i%Z_coef] = a[i];
    }
}

uint32x4_t mul_high_32(uint32x4_t a, uint32x4_t b)
{
    uint32x4_t getlow = vdupq_n_u32(0x0000ffff);
    uint32x4_t a_low = vandq_u32(a, getlow);
    uint32x4_t b_low = vandq_u32(b, getlow);
    uint32x4_t a_high = vshrq_n_u32(a, 16);
    uint32x4_t b_high = vshrq_n_u32(b, 16);

    // step1
    uint32x4_t temp1 = vmulq_u32(a_low, b_low);
    temp1 = vshrq_n_u32(temp1, 16);

    // step2
    temp1 = vmlaq_u32(temp1, a_low, b_high);
    temp1 = vmlaq_u32(temp1, a_high, b_low);
    temp1 = vshrq_n_u32(temp1, 16);

    // step3
    temp1 = vmlaq_u32(temp1, a_high, b_high);
    return temp1;
}

uint32x4_t mont(uint32x4_t a, uint32x4_t b, uint32_t n)
{
    uint32x4_t r2_mod_m_B_mtx = vdupq_n_u32(r2_mod_m_B);
    uint32x4_t r1 = mul_high_32(b, r2_mod_m_B_mtx);

    uint32x4_t r2_mod_m_Bprime_mtx = vdupq_n_u32(r2_mod_m_Bprime);
    uint32x4_t r2i = vmulq_u32(b, r2_mod_m_Bprime_mtx);

    uint32x4_t M_mtx = vdupq_n_u32(prime_q);
    uint32x4_t r2o = mul_high_32(M_mtx, r2i);

    uint32x4_t B = vsubq_u32(r1, r2o);
    if (r1[0] < r2o[0])
        B[0] += prime_q;
    if (r1[1] < r2o[1])
        B[1] += prime_q;
    if (r1[2] < r2o[2])
        B[2] += prime_q;
    if (r1[3] < r2o[3])
        B[3] += prime_q;
    

    uint32x4_t Mprime_mtx = vdupq_n_u32(Mprime);
    uint32x4_t Bprime = vmulq_u32(B, Mprime_mtx);

    r1 = mul_high_32(a, B);
    r2i = vmulq_u32(a, Bprime);
    r2o = mul_high_32(M_mtx, r2i);

    uint32x4_t result = vsubq_u32(r1, r2o);
    if (r1[0] < r2o[0])
        result[0] += prime_q;
    if (r1[1] < r2o[1])
        result[1] += prime_q;
    if (r1[2] < r2o[2])
        result[2] += prime_q;
    if (r1[3] < r2o[3])
        result[3] += prime_q;

    return result;
}

uint32_t mont_single(uint32_t a, uint32_t b, uint32_t n)
{
    uint64_t temp1, temp2;
    temp1 = b;
    temp2 = r2_mod_m_B;
    uint32_t r1 = (temp1*temp2) >> 32;

    uint32_t r2i = b*r2_mod_m_Bprime;

    temp1 = prime_q;
    temp2 = r2i;
    uint32_t r2o = (temp1*temp2) >> 32;

    uint32_t B = r1 - r2o;
    if (r1 < r2o) B += prime_q;

    uint32_t Bprime = B*Mprime;

    temp1 = a;
    temp2 = B;
    r1 = (temp1*temp2) >> 32;

    r2i = a*Bprime;

    temp1 = prime_q;
    temp2 = r2i;
    r2o = (temp1*temp2) >> 32;

    uint32_t result = r1 - r2o;
    if (r1 < r2o) result += prime_q;

    return result;
}

void swap(uint32_t &a, uint32_t &b) {
    uint32_t temp = a;
    a = b;
    b = temp;
}

void swap_all(uint32_t *seq);


// inverse = false means NTT
// inverse = true means INTT
void radix2_ntt_mr(uint32_t *seq, int len, bool inverse)
{
    swap_all(seq);

    uint32_t *roots_seq;

    roots_seq = (inverse)? root_seq_negative : root_seq_positive;

    int h = 2;
    while (h <= len)
    {
        for (int i = 0; i < len; i += h)
        {
            for (int j = 0; j < h / 2; j++)
            {
                if (h >= 8)
                {
                    uint32x4_t u = vld1q_u32(&seq[i + j]);
                    uint32x4_t temp = vld1q_u32(&seq[i + j + h / 2]);
                    uint32x4_t roots_4;
                    roots_4[0] = roots_seq[(len / h) * j];
                    roots_4[1] = roots_seq[(len / h) * (j + 1)];
                    roots_4[2] = roots_seq[(len / h) * (j + 2)];
                    roots_4[3] = roots_seq[(len / h) * (j + 3)];

                    uint32x4_t v = mont(temp, roots_4, prime_q);
                    temp = vaddq_u32(u, v);
                    if (temp[0] >= prime_q)
                        temp[0] -= prime_q;
                    if (temp[1] >= prime_q)
                        temp[1] -= prime_q;
                    if (temp[2] >= prime_q)
                        temp[2] -= prime_q;
                    if (temp[3] >= prime_q)
                        temp[3] -= prime_q;

                    vst1q_u32(&seq[i + j], temp);

                    temp = vsubq_u32(u, v);
                    if (u[0] < v[0])
                        temp[0] += prime_q;
                    if (u[1] < v[1])
                        temp[1] += prime_q;
                    if (u[2] < v[2])
                        temp[2] += prime_q;
                    if (u[3] < v[3])
                        temp[3] += prime_q;

                    vst1q_u32(&seq[i + j + h / 2], temp);

                    j += 3;
                }
                else
                {

                    uint32_t u = seq[i + j];
                    uint32_t v = mont_single(seq[i + j + h / 2], roots_seq[(len / h) * j], prime_q);

                    seq[i + j] = (u + v);
                    if (seq[i + j] >= prime_q)
                        seq[i + j] -= prime_q;

                    seq[i + j + h / 2] = (u - v);
                    if (u < v)
                        seq[i + j + h / 2] += prime_q;
                }
            }
        }
        h *= 2;
    }

    if (inverse)
    {
        for (int i = 0; i < len; i += 4)
        {
            uint32x4_t seq_vec32 = vld1q_u32(&seq[i]);
            uint32x4_t rv_inverse_mtx = vdupq_n_u32(rv_inverse);
            seq_vec32 = mont(seq_vec32, rv_inverse_mtx, prime_q);
            vst1q_u32(&seq[i], seq_vec32);
        }
    }
}

// choice 1: NTT
// choice 2: INTT
void load_radix3_variable(int choice, uint32_t **a)
{
    if (choice == 1)
    {
        a[0][0] = 1;
        a[0][1] = 1025640636;
        a[0][2] = 394137924;

        a[1][0] = 1;
        a[1][1] = 394137924;
        a[1][2] = 1025640636;

        a[2][0] = 1;
        a[2][1] = 1;
        a[2][2] = 1;
    }
    else if (choice == 2)
    {
        a[0][0] = 946519041;
        a[0][1] = 946519041;
        a[0][2] = 946519041;

        a[1][0] = 131379308;
        a[1][1] = 341880212;
        a[1][2] = 946519041;

        a[2][0] = 341880212;
        a[2][1] = 131379308;
        a[2][2] = 946519041;
    }
    else
    {
        cout << "choice not found" << endl;
    }
}

uint32_t *transform_2d_to_1d(uint32_t **a)
{
    if (Y_coef != 3 || Z_coef != 512)
    {
        cout << "transform_1d_to_2d does not support current Y_coef and Z_coef" << endl;
    }

    uint32_t *result = new uint32_t[Y_coef * Z_coef];

    for (int i = 0; i < Y_coef; i++)
    {
        for (int j = 0; j < Z_coef; j++)
        {
            result[(1024 * i + 513 * j) % pad_N] = a[i][j]%q;//q = 2048
        }
    }

    return result;
}

uint32_t **pointwise_mul_2d_mr(int row, int col, uint32_t **A, uint32_t **B)
{
    uint32_t **result = new uint32_t *[row];
    for (int i = 0; i < row; i++)
    {
        result[i] = new uint32_t[col];
    }

    if (row != 3 || col != 512)
    {
        cout << "pointwise_mul_2d_mont_red error with row or col not match" << endl;
        return result;
    }

    uint32x4_t a_vec32;
    uint32x4_t b_vec32;
    uint32x4_t c_vec32;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j += 4)
        {
            a_vec32 = vld1q_u32(&A[i][j]);
            b_vec32 = vld1q_u32(&B[i][j]);
            c_vec32 = mont(a_vec32, b_vec32, prime_q);

            vst1q_u32(&result[i][j], c_vec32);
        }
    }
    return result;
}

void radix3_mul_mont(uint32_t **var, uint32_t **A, int len)
{
    uint32x4_t row0 = vld1q_u32(var[0]);
    uint32x4_t row1 = vld1q_u32(var[1]);
    uint32x4_t row2 = vld1q_u32(var[2]);
    uint32x4_t col;
    row0[3] = 0;
    row1[3] = 0;
    row2[3] = 0;
    col[3] = 0;

    for (int i = 0; i < len; i++)
    {
        col[0] = A[0][i];
        col[1] = A[1][i];
        col[2] = A[2][i];

        A[0][i] = vaddvq_u32(mont(row0, col, prime_q));
        A[1][i] = vaddvq_u32(mont(row1, col, prime_q));
        A[2][i] = vaddvq_u32(mont(row2, col, prime_q));

        if (A[0][i] >= prime_q)
            A[0][i] -= prime_q;

        if (A[1][i] >= prime_q)
            A[1][i] -= prime_q;
    
        if (A[2][i] >= prime_q)
            A[2][i] -= prime_q;
    }
}

void swap_all(uint32_t *seq) {
    swap(seq[1], seq[256]);
    swap(seq[2], seq[128]);
    swap(seq[3], seq[384]);
    swap(seq[4], seq[64]);
    swap(seq[5], seq[320]);
    swap(seq[6], seq[192]);
    swap(seq[7], seq[448]);
    swap(seq[8], seq[32]);
    swap(seq[9], seq[288]);
    swap(seq[10], seq[160]);
    swap(seq[11], seq[416]);
    swap(seq[12], seq[96]);
    swap(seq[13], seq[352]);
    swap(seq[14], seq[224]);
    swap(seq[15], seq[480]);
    swap(seq[17], seq[272]);
    swap(seq[18], seq[144]);
    swap(seq[19], seq[400]);
    swap(seq[20], seq[80]);
    swap(seq[21], seq[336]);
    swap(seq[22], seq[208]);
    swap(seq[23], seq[464]);
    swap(seq[24], seq[48]);
    swap(seq[25], seq[304]);
    swap(seq[26], seq[176]);
    swap(seq[27], seq[432]);
    swap(seq[28], seq[112]);
    swap(seq[29], seq[368]);
    swap(seq[30], seq[240]);
    swap(seq[31], seq[496]);
    swap(seq[33], seq[264]);
    swap(seq[34], seq[136]);
    swap(seq[35], seq[392]);
    swap(seq[36], seq[72]);
    swap(seq[37], seq[328]);
    swap(seq[38], seq[200]);
    swap(seq[39], seq[456]);
    swap(seq[41], seq[296]);
    swap(seq[42], seq[168]);
    swap(seq[43], seq[424]);
    swap(seq[44], seq[104]);
    swap(seq[45], seq[360]);
    swap(seq[46], seq[232]);
    swap(seq[47], seq[488]);
    swap(seq[49], seq[280]);
    swap(seq[50], seq[152]);
    swap(seq[51], seq[408]);
    swap(seq[52], seq[88]);
    swap(seq[53], seq[344]);
    swap(seq[54], seq[216]);
    swap(seq[55], seq[472]);
    swap(seq[57], seq[312]);
    swap(seq[58], seq[184]);
    swap(seq[59], seq[440]);
    swap(seq[60], seq[120]);
    swap(seq[61], seq[376]);
    swap(seq[62], seq[248]);
    swap(seq[63], seq[504]);
    swap(seq[65], seq[260]);
    swap(seq[66], seq[132]);
    swap(seq[67], seq[388]);
    swap(seq[69], seq[324]);
    swap(seq[70], seq[196]);
    swap(seq[71], seq[452]);
    swap(seq[73], seq[292]);
    swap(seq[74], seq[164]);
    swap(seq[75], seq[420]);
    swap(seq[76], seq[100]);
    swap(seq[77], seq[356]);
    swap(seq[78], seq[228]);
    swap(seq[79], seq[484]);
    swap(seq[81], seq[276]);
    swap(seq[82], seq[148]);
    swap(seq[83], seq[404]);
    swap(seq[85], seq[340]);
    swap(seq[86], seq[212]);
    swap(seq[87], seq[468]);
    swap(seq[89], seq[308]);
    swap(seq[90], seq[180]);
    swap(seq[91], seq[436]);
    swap(seq[92], seq[116]);
    swap(seq[93], seq[372]);
    swap(seq[94], seq[244]);
    swap(seq[95], seq[500]);
    swap(seq[97], seq[268]);
    swap(seq[98], seq[140]);
    swap(seq[99], seq[396]);
    swap(seq[101], seq[332]);
    swap(seq[102], seq[204]);
    swap(seq[103], seq[460]);
    swap(seq[105], seq[300]);
    swap(seq[106], seq[172]);
    swap(seq[107], seq[428]);
    swap(seq[109], seq[364]);
    swap(seq[110], seq[236]);
    swap(seq[111], seq[492]);
    swap(seq[113], seq[284]);
    swap(seq[114], seq[156]);
    swap(seq[115], seq[412]);
    swap(seq[117], seq[348]);
    swap(seq[118], seq[220]);
    swap(seq[119], seq[476]);
    swap(seq[121], seq[316]);
    swap(seq[122], seq[188]);
    swap(seq[123], seq[444]);
    swap(seq[125], seq[380]);
    swap(seq[126], seq[252]);
    swap(seq[127], seq[508]);
    swap(seq[129], seq[258]);
    swap(seq[131], seq[386]);
    swap(seq[133], seq[322]);
    swap(seq[134], seq[194]);
    swap(seq[135], seq[450]);
    swap(seq[137], seq[290]);
    swap(seq[138], seq[162]);
    swap(seq[139], seq[418]);
    swap(seq[141], seq[354]);
    swap(seq[142], seq[226]);
    swap(seq[143], seq[482]);
    swap(seq[145], seq[274]);
    swap(seq[147], seq[402]);
    swap(seq[149], seq[338]);
    swap(seq[150], seq[210]);
    swap(seq[151], seq[466]);
    swap(seq[153], seq[306]);
    swap(seq[154], seq[178]);
    swap(seq[155], seq[434]);
    swap(seq[157], seq[370]);
    swap(seq[158], seq[242]);
    swap(seq[159], seq[498]);
    swap(seq[161], seq[266]);
    swap(seq[163], seq[394]);
    swap(seq[165], seq[330]);
    swap(seq[166], seq[202]);
    swap(seq[167], seq[458]);
    swap(seq[169], seq[298]);
    swap(seq[171], seq[426]);
    swap(seq[173], seq[362]);
    swap(seq[174], seq[234]);
    swap(seq[175], seq[490]);
    swap(seq[177], seq[282]);
    swap(seq[179], seq[410]);
    swap(seq[181], seq[346]);
    swap(seq[182], seq[218]);
    swap(seq[183], seq[474]);
    swap(seq[185], seq[314]);
    swap(seq[187], seq[442]);
    swap(seq[189], seq[378]);
    swap(seq[190], seq[250]);
    swap(seq[191], seq[506]);
    swap(seq[193], seq[262]);
    swap(seq[195], seq[390]);
    swap(seq[197], seq[326]);
    swap(seq[199], seq[454]);
    swap(seq[201], seq[294]);
    swap(seq[203], seq[422]);
    swap(seq[205], seq[358]);
    swap(seq[206], seq[230]);
    swap(seq[207], seq[486]);
    swap(seq[209], seq[278]);
    swap(seq[211], seq[406]);
    swap(seq[213], seq[342]);
    swap(seq[215], seq[470]);
    swap(seq[217], seq[310]);
    swap(seq[219], seq[438]);
    swap(seq[221], seq[374]);
    swap(seq[222], seq[246]);
    swap(seq[223], seq[502]);
    swap(seq[225], seq[270]);
    swap(seq[227], seq[398]);
    swap(seq[229], seq[334]);
    swap(seq[231], seq[462]);
    swap(seq[233], seq[302]);
    swap(seq[235], seq[430]);
    swap(seq[237], seq[366]);
    swap(seq[239], seq[494]);
    swap(seq[241], seq[286]);
    swap(seq[243], seq[414]);
    swap(seq[245], seq[350]);
    swap(seq[247], seq[478]);
    swap(seq[249], seq[318]);
    swap(seq[251], seq[446]);
    swap(seq[253], seq[382]);
    swap(seq[255], seq[510]);
    swap(seq[259], seq[385]);
    swap(seq[261], seq[321]);
    swap(seq[263], seq[449]);
    swap(seq[265], seq[289]);
    swap(seq[267], seq[417]);
    swap(seq[269], seq[353]);
    swap(seq[271], seq[481]);
    swap(seq[275], seq[401]);
    swap(seq[277], seq[337]);
    swap(seq[279], seq[465]);
    swap(seq[281], seq[305]);
    swap(seq[283], seq[433]);
    swap(seq[285], seq[369]);
    swap(seq[287], seq[497]);
    swap(seq[291], seq[393]);
    swap(seq[293], seq[329]);
    swap(seq[295], seq[457]);
    swap(seq[299], seq[425]);
    swap(seq[301], seq[361]);
    swap(seq[303], seq[489]);
    swap(seq[307], seq[409]);
    swap(seq[309], seq[345]);
    swap(seq[311], seq[473]);
    swap(seq[315], seq[441]);
    swap(seq[317], seq[377]);
    swap(seq[319], seq[505]);
    swap(seq[323], seq[389]);
    swap(seq[327], seq[453]);
    swap(seq[331], seq[421]);
    swap(seq[333], seq[357]);
    swap(seq[335], seq[485]);
    swap(seq[339], seq[405]);
    swap(seq[343], seq[469]);
    swap(seq[347], seq[437]);
    swap(seq[349], seq[373]);
    swap(seq[351], seq[501]);
    swap(seq[355], seq[397]);
    swap(seq[359], seq[461]);
    swap(seq[363], seq[429]);
    swap(seq[367], seq[493]);
    swap(seq[371], seq[413]);
    swap(seq[375], seq[477]);
    swap(seq[379], seq[445]);
    swap(seq[383], seq[509]);
    swap(seq[391], seq[451]);
    swap(seq[395], seq[419]);
    swap(seq[399], seq[483]);
    swap(seq[407], seq[467]);
    swap(seq[411], seq[435]);
    swap(seq[415], seq[499]);
    swap(seq[423], seq[459]);
    swap(seq[431], seq[491]);
    swap(seq[439], seq[475]);
    swap(seq[447], seq[507]);
    swap(seq[463], seq[487]);
    swap(seq[479], seq[503]);
}

#endif
