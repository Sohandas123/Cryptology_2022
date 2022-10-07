#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>

int **IP, **InvIP;           //IP : Initial Permutation matrix & InvIP : Matrix of Inverse of Initial Permutation
int **E;                     // Expansion matrix
unsigned char **L, **R;      // for L0,L1,...,L16 and R0,R1,... ,R16
unsigned char **ER;          // after applying Expansion() of R0.R1,....,R16
unsigned char **K;           // K[0] will be the 64-bit key but rest of 16s will be round keys
unsigned char *C, *D;        // C : Left 28-bit of round key & D : Right 28-bit of round key 
unsigned char *ciphertxt;

void Apply_Per_PI(unsigned char *, unsigned char *, int **, int, int);  /* (before permutation, after permutation, permutation matrix, row-size, column-size ) */
void Expansion(unsigned char *, unsigned char *, int **);    /* (before expansion, after expansion, expansion matrix ) */
void makeInputS(unsigned char *, unsigned char *);   /*  */
void Apply_Sbox(unsigned char, unsigned char *, int (*S)[16]);
void modify_Sout(unsigned char *, unsigned char *);
void Apply_P_in_f(unsigned char *, unsigned char *);
void XOR(unsigned char *, unsigned char *, unsigned char *, int); /* (input1, input2, output, size) */
void trans56to28(unsigned char *, unsigned char *, unsigned char *);
void leftshift(unsigned char*,unsigned char*);
void join_CD(unsigned char *, unsigned char *, unsigned char *);

void round_Key(int);

void Feistel(int);

void DES(unsigned char plaintxt[8], unsigned char key[8], unsigned char *ciphertxt);



/* PC_1 permutation for key scheduling */
int **PC1; /* due to some function implementation we use it same as PC_1 */
int PC_1[7][8] = {
    {57, 49, 41, 33, 25, 17, 9, 1},
    {58, 50, 42, 34, 26, 18, 10, 2},
    {59, 51, 43, 35, 27, 19, 11, 3},
    {60, 52, 44, 36, 63, 55, 47, 39},
    {31, 23, 15, 7, 62, 54, 46, 38},
    {30, 22, 14, 6, 61, 53, 45, 37},
    {29, 21, 13, 5, 28, 20, 12, 4}};

/* PC_2 permutation for key scheduling */
int **PC2; /* due to some function implementation we use it same as PC_1 */
int PC_2[6][8] = {
    {14, 17, 11, 24, 1, 5, 3, 28},
    {15, 6, 21, 10, 23, 19, 12, 4},
    {26, 8, 16, 7, 27, 20, 13, 2},
    {41, 52, 31, 37, 47, 55, 30, 40},
    {51, 45, 33, 48, 44, 49, 39, 56},
    {34, 53, 46, 42, 50, 36, 29, 32}};

/* S_box_1 */
int S1[4][16] = {
    {14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7},
    {0, 15, 7, 4, 14, 2, 13, 1, 10, 6, 12, 11, 9, 5, 3, 8},
    {4, 1, 14, 8, 13, 6, 2, 11, 15, 12, 9, 7, 3, 10, 5, 0},
    {15, 12, 8, 2, 4, 9, 1, 7, 5, 11, 3, 14, 10, 0, 6, 13}};
/* S_box_2 */
int S2[4][16] = {
    {15, 1, 8, 14, 6, 11, 3, 4, 9, 7, 2, 13, 12, 0, 5, 10},
    {3, 13, 4, 7, 15, 2, 8, 14, 12, 0, 1, 10, 6, 9, 11, 5},
    {0, 14, 7, 11, 10, 4, 13, 1, 5, 8, 12, 6, 9, 3, 2, 15},
    {13, 8, 10, 1, 3, 15, 4, 2, 11, 6, 7, 12, 0, 5, 14, 9}};
/* S_box_3 */
int S3[4][16] = {
    {10, 0, 9, 14, 6, 3, 15, 5, 1, 13, 12, 7, 11, 4, 2, 8},
    {13, 7, 0, 9, 3, 4, 6, 10, 2, 8, 5, 14, 12, 11, 15, 1},
    {13, 6, 4, 9, 8, 15, 3, 0, 11, 1, 2, 12, 5, 10, 14, 7},
    {1, 10, 13, 0, 6, 9, 8, 7, 4, 15, 14, 3, 11, 5, 2, 12}};

/* S_box_4 */
int S4[4][16] = {
    {7, 13, 14, 3, 0, 6, 9, 10, 1, 2, 8, 5, 11, 12, 4, 15},
    {13, 8, 11, 5, 6, 15, 0, 3, 4, 7, 2, 12, 1, 10, 14, 9},
    {10, 6, 9, 0, 12, 11, 7, 13, 15, 1, 3, 14, 5, 2, 8, 4},
    {3, 15, 0, 6, 10, 1, 13, 8, 9, 4, 5, 11, 12, 7, 2, 14}};
/* S_box_5 */
int S5[4][16] = {
    {2, 12, 4, 1, 7, 10, 11, 6, 8, 5, 3, 15, 13, 0, 14, 9},
    {14, 11, 2, 12, 4, 7, 13, 1, 5, 0, 15, 10, 3, 9, 8, 6},
    {4, 2, 1, 11, 10, 13, 7, 8, 15, 9, 12, 5, 6, 3, 0, 14},
    {11, 8, 12, 7, 1, 14, 2, 13, 6, 15, 0, 9, 10, 4, 5, 3}};
/* S_box_6 */
int S6[4][16] = {
    {12, 1, 10, 15, 9, 2, 6, 8, 0, 13, 3, 4, 14, 7, 5, 11},
    {10, 15, 4, 2, 7, 12, 9, 5, 6, 1, 13, 14, 0, 11, 3, 8},
    {9, 14, 15, 5, 2, 8, 12, 3, 7, 0, 4, 10, 1, 13, 11, 6},
    {4, 3, 2, 12, 9, 5, 15, 10, 11, 14, 1, 7, 6, 0, 8, 13}};
/* S_box_7 */
int S7[4][16] = {
    {4, 11, 2, 14, 15, 0, 8, 13, 3, 12, 9, 7, 5, 10, 6, 1},
    {13, 0, 11, 7, 4, 9, 1, 10, 14, 3, 5, 12, 2, 15, 8, 6},
    {1, 4, 11, 13, 12, 3, 7, 14, 10, 15, 6, 8, 0, 5, 9, 2},
    {6, 11, 13, 8, 1, 4, 10, 7, 9, 5, 0, 15, 14, 2, 3, 12}};
/* S_box_8 */
int S8[4][16] = {
    {13, 2, 8, 4, 6, 15, 11, 1, 10, 9, 3, 14, 5, 0, 12, 7},
    {1, 15, 13, 8, 10, 3, 7, 4, 12, 5, 6, 11, 0, 14, 9, 2},
    {7, 11, 4, 1, 9, 12, 14, 2, 0, 6, 10, 13, 15, 3, 5, 8},
    {2, 1, 14, 7, 4, 10, 8, 13, 15, 12, 9, 0, 3, 5, 6, 11}};

/*Permutation P, used in Feistel function */
int P_in_f[4][8] = {
    {16, 7, 20, 21, 29, 12, 28, 17},
    {1, 15, 23, 26, 5, 18, 31, 10},
    {2, 8, 24, 14, 32, 27, 3, 9},
    {19, 13, 30, 6, 22, 11, 4, 25}};





int main()
{   
    int i,j, m = 0;
    //  unsigned char plaintxt[8] = "SohanDas";
    unsigned char plaintxt[9];  /*    = {'S', 'o', 'h', 'a', 'n', 'D', 'a', 's'}; */
    unsigned char key[8], encmsg[8], IV[8]  = "hellohii";
    ciphertxt = (unsigned char *)calloc(8, sizeof(unsigned char));
    for(i=0; i<8; i++){
        key[i] = 0x69;
        ciphertxt[i] = 0x00;
        // printf("IV : %x | key : %x | ciphertxt ; %x\n",IV[i],key[i],ciphertxt[i]);
    }

    int fd = open("msg.txt",O_RDONLY);
    int fd_enc = open("EncryptedMsg.txt",O_WRONLY | O_CREAT | O_TRUNC, 0666);
    while(1) {
        for(i=0; i<8; i++){
            plaintxt[i] = 0x00;
            ciphertxt[i] = 0x00;
        }
        read(fd, plaintxt, 8*sizeof(unsigned char));
        if(strcmp(plaintxt,ciphertxt)==0)     /* while loop breaking condition */
            break;
        printf("\n Plaintxt : ");
        for(i=0; i<8; i++)
            printf(" %c",plaintxt[i]);
        printf("\n Ciphertxt : ");
        DES(IV, key, ciphertxt);
        XOR(plaintxt, ciphertxt, encmsg, 8);
        write(fd_enc, encmsg, 8*sizeof(unsigned char));
        for(i=0; i<8; i++)
            printf(" %c",encmsg[i]);
        printf("\n");
        for(i=0; i<8; i++) {
            IV[i] = ciphertxt[i];
        }
    }
    close(fd);
    close(fd_enc);


    // /* ---------------------- Decoration Part---------------------- */
    // printf("\n   Plaintext   :");
    // for (i = 0; i<8; i += 1){
    //     printf("  %c ",plaintxt[i]);
    // }                                    /* Plaintext printing */
    // printf("\n (hexadecimal) :");
    // for (i = 0; i<8; i += 1){
    //     printf(" %x ",plaintxt[i]);
    // }
    // /* -----------------------------------------------*/
    // printf("\n Private Key\n (hexadecimal) :");
    // for (i = 0; i<8; i += 1){
    //     printf("  %x ",key[i]);           /* Private Key Printing */
    // }
    // /* ----------------------------------------------- */
    // printf("\n  Ciphertext\n (hexadecimal) :");
    // for (i = 0; i<8; i++){
    //     printf(" %x ",ciphertxt[i]);    /* Ciphertext Printing */
    // }
    // printf("\n");

    return 0;
}






/* -----------------------------------Function Definitions -------------------------------------------*/


void XOR(unsigned char *inp_1, unsigned char *inp_2, unsigned char *opt, int size) /* (input1, input2, output, size) */
{
    int i;
    for (i = 0; i < size; i++)
    {
        opt[i] = inp_1[i] ^ inp_2[i];
    }
}

void Apply_Per_PI(unsigned char *pl, unsigned char *P_pl, int **PI, int row, int column)
{
    int i, j, k, m, n ;
    unsigned char b;
    for (i = 0; i < row; i++)
    {
        // printf("%d. \n",i);
        P_pl[i] = 0x00;
        for (j = 0; j < column; j++)
        {
            k = PI[i][j] - 1; 
            m = k / 8;     
            n = 7 - k % 8;
            b = (pl[m] >> n) & 1; 
            P_pl[i] = P_pl[i] | (b << (7 - j));
        }
        // printf(" \n");
    }
}

void Expansion(unsigned char *R, unsigned char *MR, int **M)
{
    int i, j, k, m, n, total = -1;
    unsigned char b;
    for (i = 0; i < 8; i++)
    {
        for (j = 0; j < 6; j++)
        {
            k = M[i][j] - 1; /* Similar to Apply_Per_P() */
            m = k / 8;
            n = 7 - k % 8; /* Working Correctly */
            b = (R[m] >> n) & 1;

            total += 1;
            m = total / 8;
            n = 7 - total % 8;
            MR[m] = MR[m] | (b << n);
        }
    }
}

void makeInputS(unsigned char *E_R, unsigned char *inpt_S) // working fine
{
    inpt_S[0] = 0x00 | ((E_R[0] & 0xfc) >> 2);      // 0xfc = 1111 1100
    inpt_S[1] = 0x00 | ((E_R[0] & 0x03) << 4);      // 0x03 = 0000 0011
    inpt_S[1] = inpt_S[1] | ((E_R[1] & 0xf0) >> 4); // 0xf0 = 1111 0000
    inpt_S[2] = 0x00 | ((E_R[1] & 0x0f) << 2);      // 0x0f = 0000 1111
    inpt_S[2] = inpt_S[2] | ((E_R[2] & 0xc0) >> 6); // 0xc0 = 1100 0000
    inpt_S[3] = 0x00 | ((E_R[2] & 0x3f));           // 0x3f = 0011 1111

    inpt_S[4] = 0x00 | ((E_R[3] & 0xfc) >> 2);      // 0xfc = 1111 1100
    inpt_S[5] = 0x00 | ((E_R[3] & 0x03) << 4);      // 0x03 = 0000 0011
    inpt_S[5] = inpt_S[1] | ((E_R[4] & 0xf0) >> 4); // 0xf0 = 1111 0000
    inpt_S[6] = 0x00 | ((E_R[4] & 0x0f) << 2);      // 0x0f = 0000 1111
    inpt_S[6] = inpt_S[6] | ((E_R[5] & 0xc0) >> 6); // 0xc0 = 1100 0000
    inpt_S[7] = 0x00 | ((E_R[5] & 0x3f));           // 0x3f = 0011 1111
}

void Apply_Sbox(unsigned char in, unsigned char *out, int (*S)[16]) // Working Fine
{
    *out = 0x00;
    unsigned char m, n;                   // m - row, n- column
    n = (in & 0x1e) >> 1;                 // 0x1e = 0001 1110
    m = (in & 0x01) | ((in & 0x20) >> 4); // 0x01 = 0000 0001, 0x20 =0010 0000
    *out = S[m][n];
    // printf("%x",*out);
}

void modify_Sout(unsigned char *in_4bit_8char, unsigned char *out_8bit_4char)
{
    int i;
    unsigned char b = 0x00; // working fine
    for (i = 0; i < 4; i++)
    {
        out_8bit_4char[i] = 0x00;               // initialization (not needed)
        b = (in_4bit_8char[2 * i] & 0x0f) << 4; // 0x0f = 0000 1111
        out_8bit_4char[i] = b | ((in_4bit_8char[2 * i + 1]) & 0x0f);
    }
}

void Apply_P_in_f(unsigned char *before_P, unsigned char *after_P)
{
    int i, j, k, m, n;
    unsigned char b;
    for (i = 0; i < 4; i++)
    {
        // printf("%d. ",i);
        after_P[i] = 0x00;
        for (j = 0; j < 8; j++)
        {
            k = P_in_f[i][j] - 1; /* Since in P there are elements between 1 and 32 not from 0 to 31 */
            m = k / 8;            // row or block
            n = 7 - k % 8;
            b = (before_P[m] >> n) & 1; // Extracting the bit value 0/1
            // printf(" %d. m = %d, n = %d, b = %x\n",j,m,n,b);
            after_P[i] = after_P[i] | (b << (7 - j));
        }
        // printf("\n");
    }
}

void trans56to28(unsigned char *input56, unsigned char *out28C, unsigned char *out28D)
{
    int i;
    for (i = 0; i < 3; i++)
    {
        out28C[2 * i] = (input56[i] >> 4) & 0x0f; // 0x0f = 0000 1111
        out28C[2 * i + 1] = input56[i] & 0x0f;    // 0x0f = 0000 1111
    }
    out28C[6] = (input56[3] >> 4) & 0x0f;

    out28D[0] = input56[3] & 0x0f;
    for (i = 1; i <= 3; i++)
    {
        out28D[2 * i - 1] = (input56[i + 3] >> 4) & 0x0f;
        out28D[2 * i] = input56[i + 3] & 0x0f;
    }
}

void leftshift(unsigned char *C, unsigned char *D)
{

    unsigned char bc, bd ;
    int i;
    bc = (C[0]>>3) & 1;
    bd = (D[0]>>3) & 1;
    for (i=6; i>=0; i-- )
    {
        C[i] = (C[i]<<1) | bc ;
        bc = (C[i]>>4) & 1;
        C[i] = C[i] & 15;

        D[i] = (D[i]<<1) | bd ;
        bd = (D[i]>>4) & 1;
        D[i] = D[i] & 15;
    }
}

void join_CD(unsigned char* S1, unsigned char* S2, unsigned char* S1S2)
{
    S1S2[0] = (S1[0] << 4) | S1[1];
    S1S2[1] = (S1[2] << 4) | S1[3];
    S1S2[2] = (S1[4] << 4) | S1[5];
    S1S2[3] = (S1[6] << 4) | S2[0];

    S1S2[4] = (S2[1] << 4) | S2[2];
    S1S2[5] = (S2[3] << 4) | S2[4];
    S1S2[6] = (S2[5] << 4) | S2[6];
}






void round_Key(int t)
{
    int i, j;
    /* Left shifting by one bit */
    if((t!=1)&&(t!=2)&&(t!=9)&&(t!=16))
        leftshift(C,D);    /* Since for round 1,2,9,16 there is only one left shift in DES */
    leftshift(C,D); /* At t-round key generation C[t-1] left-shifted to C[t] and D[t-1] left shifted to D[t] */
    unsigned char *CD;
    CD = (unsigned char *)calloc(7, sizeof(unsigned char));
    /* We make a 56-bit string of 7 unsigned char's (8-bits each) by concatenating Cand D */
    join_CD(C, D, CD);

    /* Apply permutation PC_2 on this 56-bit CD (7 unsigned char of 8-bits) */
    Apply_Per_PI(CD, K[t], PC2, 6, 8);
}








void Feistel(int t)
{

    int i, j;
    Expansion(R[t - 1], ER[t - 1], E); // Expansion of R[i] is ER[i]
    // for (i = 0; i<6; i++)
    //     printf("%x, ",ER[0][i]);
    // printf("\n");

    /* ------------------------------------------------------------- */

    /* For Printing S_boxes */
    // for (i = 0; i<4; i++){
    //     for(j=0;j<16;j++)
    //         printf("%d, ",S1[i][j]);
    // printf("\n");
    // }

    /* ---------------------------------------------------------------------------------------------------------  */
    /* do XOR of ER[t-1] and round_Key[t] and then feed it to S boxes */
//    printf(" Round key %d : ", t);
    round_Key(t); // round key K[t] generated and now it can be used
/*    for (i = 0; i < 6; i++)
        printf(" %x", K[t][i]);
    printf("\n");  */

    unsigned char *XORed_Key_ER;
    XORed_Key_ER = (unsigned char *)calloc(6, sizeof(unsigned char));
    XOR(ER[t - 1], K[t], XORed_Key_ER, 6);
    /* --------------------------------------------------------------------------------------------------------  */

    /* Break the ER character array of length 6 in another acharacter array of length of 8 */

    unsigned char *inputofS;
    inputofS = (unsigned char *)calloc(8, sizeof(unsigned char));
    makeInputS(XORed_Key_ER, inputofS); // this will make 6 bit input for S boxes as inputofS
    // for(i=0;i<8;i++)
    //     printf(" %x ",inputofS[i]);

    /* output of S box */
    unsigned char *output_S = calloc(8, sizeof(unsigned char));
    Apply_Sbox(inputofS[0], &output_S[0], S1);
    // printf("%x",output_S[0]);
    Apply_Sbox(inputofS[1], &output_S[1], S2);
    Apply_Sbox(inputofS[2], &output_S[2], S3);
    Apply_Sbox(inputofS[3], &output_S[3], S4);
    Apply_Sbox(inputofS[4], &output_S[4], S5);
    Apply_Sbox(inputofS[5], &output_S[5], S6);
    Apply_Sbox(inputofS[6], &output_S[6], S7);
    Apply_Sbox(inputofS[7], &output_S[7], S8);
    // for(i=0;i<8;i++)
    //     printf(" %x ",output_S[i]);

    /* modify output of S boxes i.e from these 8 output_S's (of 4-bit), make 4 8-bit blocks(of unsigned char ) */
    unsigned char *mod_out_S = calloc(4, sizeof(unsigned char));
    modify_Sout(output_S, mod_out_S);
    // for(i=0;i<4;i++)
    //     printf(" %x ",mod_out_S[i]);

    /* Apply permutation P on modified output of S box */
    unsigned char *pre_XOR_L = calloc(4, sizeof(unsigned char));
    Apply_P_in_f(mod_out_S, pre_XOR_L);
    // for(i=0;i<4;i++)
    //     printf(" %x ",pre_XOR_L[i]);

    unsigned char *post_XOR_L = calloc(4, sizeof(unsigned char));
    XOR(pre_XOR_L, L[t - 1], post_XOR_L, 4);
    // for(i=0;i<4;i++)
    //     printf(" %x ",post_XOR_L[i]);
    /* -------------------------------------------------------- */

    /* Final swaping */
//    printf("\n   L[%d]    R[%d]\n   ----    -----\n", t, t);
    for (int k = 0; k < 4; k++)
    {
        R[t][k] = post_XOR_L[k];
        L[t][k] = R[t - 1][k];
//        printf("    %x      %x\n", L[t][k], R[t][k]);
    }
}



void DES(unsigned char plaintext[8], unsigned char key[8], unsigned char *ciphertxt)
{
    /* Storing values of Initial Permutation IP  */
    int i, j;
    IP = (int **)malloc(8 * sizeof(int *));
    for (i = 0; i < 8; i++)
        IP[i] = (int *)malloc(8 * sizeof(int));
    int temp_IP = 0;
    for (j = 7; j >= 0; j--)
    {
        for (i = 0; i < 4; i++)
        {
            temp_IP += 2;
            IP[i][j] = temp_IP;
        }
    }
    temp_IP = -1;
    for (j = 7; j >= 0; j--)
    {
        for (i = 4; i < 8; i++)
        {
            temp_IP += 2;
            IP[i][j] = temp_IP;
        }
    }

    /* Printing IP : */
    // printf("\nIP :\n");
    // for (i = 0; i < 8; i++){
    //     for (j = 0; j < 8; j++)
    //         printf("%d ", IP[i][j]);
    //     printf("\n");
    // }

    /* --------------------------------------------- */

    /* Storing of values of Inverse of IP */

    InvIP = (int **)malloc(8 * sizeof(int *));
    for (i = 0; i < 8; i++)
        InvIP[i] = (int *)malloc(8 * sizeof(int));
    int temp_InvIP = 0;
    for (j = 1; j < 8; j = j + 2)
    {
        for (i = 7; i >= 0; i--)
        {
            temp_InvIP += 1;
            InvIP[i][j] = temp_InvIP;
        }
    }

    for (j = 0; j < 8; j = j + 2)
    {
        for (i = 7; i >= 0; i--)
        {
            temp_InvIP += 1;
            InvIP[i][j] = temp_InvIP;
        }
    }

    /* Printing InvIP : */
    // printf("\n InvIP :\n");
    // for (i = 0; i < 8; i++){
    //     for (j = 0; j < 8; j++)
    //         printf("%d ", InvIP[i][j]);
    //     printf("\n");
    // }

    /* ----------------------------------------------------- */





    unsigned char IP_Ptxt[8];  /* for storing string after applying permutation IP on plaintext */
    Apply_Per_PI(plaintext, IP_Ptxt, IP, 8, 8); /* Applying IP on Plaintext */
    // for (i = 0; i<8; i++)
    //     printf("%d, ",IP_Ptxt[i]);
    // printf("\n");

    // unsigned char check[8];
    // Apply_Per_PI(IP_Ptxt,check,InvIP,8,8);  /* For checking whether the Apply_Per_PI() function is working correctly or not. */
    // for (i = 0; i<8; i++){
    //     printf("%c, ",check[i]);
    // }

    /* ------------------------------------------------------ */

    /* Declaring L[0],L[1],....L[16] and R[0],R[1],...,R[16] , each are of 32 bit length, so we take them as character array of size 4, since 4*8=32 */
    L = (unsigned char **)malloc(17 * sizeof(unsigned char *));
    R = (unsigned char **)malloc(17 * sizeof(unsigned char *));
    for (i = 0; i <= 16; i++)
    {
        L[i] = (unsigned char *)calloc(4, sizeof(unsigned char));
        R[i] = (unsigned char *)calloc(4, sizeof(unsigned char));
    }
    /* Initialization of L[0] and R[0]  */
    // printf("\n Round : 0\n\n   L[0]    R[0]\n   ----    -----\n");
    for (i = 0; i < 4; i++)
    {
        L[0][i] = IP_Ptxt[i];
        R[0][i] = IP_Ptxt[4 + i];

        // printf("    %x      %x\n",L[0][i],R[0][i]);
    }

    /* Expansion matrix */
    E = (int **)malloc(8 * sizeof(int *));
    for (i = 0; i < 8; i++)
        E[i] = (int *)malloc(6 * sizeof(int));

    E[0][0] = 0;
    for (i = 1; i < 8; i++)
        E[i][0] = E[i - 1][0] + 4;
    for (i = 0; i < 8; i++)
        for (j = 1; j < 6; j++)
            E[i][j] = E[i][j - 1] + 1;
    E[0][0] = 32;
    E[7][5] = 1;
    // for (i = 0; i < 8; i++){
    //     for (j = 0; j < 6; j++)
    //         printf("%d ", E[i][j]);
    //     printf("\n");
    // }

    /* ER means expanded R from 32-bit to 48-bit   */
    ER = (unsigned char **)malloc(17 * sizeof(unsigned char *));
    for (i = 0; i <= 16; i++)
    {
        ER[i] = (unsigned char *)calloc(6, sizeof(unsigned char));
    }

    /* Key Scheduling */
    /* We have to generate 16 round keys of 48-bit from a 64-bit given key and then do XOR with respective ERs. We will make these round keys as global for further checking at the time of reverse key scheduling */

    K = (unsigned char **)malloc(17 * sizeof(unsigned char *));
    K[0] = (unsigned char *)malloc(8 * sizeof(unsigned char));
    for (i = 0; i < 8; i++)
    {
        K[0][i] = key[i]; /* Just for now to check correctness */
        //  printf("\n %x",K[0][i]);
    }
    for (i = 1; i <= 16; i++)
    {
        K[i] = (unsigned char *)calloc(6, sizeof(unsigned char)); /* same as ER[i] as we will XOR it with that */
    }

    unsigned char *afterPC_1;
    afterPC_1 = (unsigned char *)malloc(7 * sizeof(unsigned char)); /* 56-bit output of PC_1 will be stored in 7 unsigned char */
    /* Due to type casting we have to copy PC_1[][] 2-D array data to a double-pointer PC1 */
    PC1 = (int **)malloc(7 * sizeof(int *));
    for (i = 0; i < 7; i++)
    {
        PC1[i] = (int *)malloc(8 * sizeof(int));
        // printf("\n");
        for (j = 0; j < 8; j++)
        {
            PC1[i][j] = PC_1[i][j];
            // printf(" %d",PC1[i][j]);
        }
    }

    Apply_Per_PI(K[0], afterPC_1, PC1, 7, 8);

    /* Make C and D , two 28-bit of 56-bit afterPC_1 */
    C = (unsigned char *)calloc(7, sizeof(unsigned char));
    D = (unsigned char *)calloc(7, sizeof(unsigned char));
    
    trans56to28(afterPC_1, C, D);
    /* Due to type casting we have to copy PC_2[][] 2-D array data to a double-pointer PC2 */
    PC2 = (int **)malloc(7 * sizeof(int *));
    for (i = 0; i < 6; i++)
    {
        PC2[i] = (int *)malloc(8 * sizeof(int));
        // printf("\n");
        for (j = 0; j < 8; j++)
        {
            PC2[i][j] = PC_2[i][j];
            // printf(" %d",PC2[i][j]);
        }
    }

    // for (int t = 1; t <= 16; t++)
    // {
    //     // printf(" Round key %d : ", t);
    //     round_Key(t); /* round key K[t] generated and now it can be used */
    //     // for (i = 0; i <= 16; i++)
    //     //     printf(" %x", K[t][i]);
    //     // printf("\n");
    // }

    /* --------------------------------------------------------------------- */

    /* 16-round Feistel network  */
    for (i = 1; i <= 16; i++)
    {
//        printf("\n------------------------\n Round : %d\n\n", i);
        Feistel(i);
    }

    /* Finally invIP (inverse of initial permutation) on R[16]_L[16] */
    // strcat(R[16],L[16]);
    unsigned char *pre_cipher = calloc(8, sizeof(unsigned char));
    for (i = 0; i < 8; i++)
    {
        if (i < 4)
            pre_cipher[i] = R[16][i];
        else
            pre_cipher[i] = L[16][i - 4];
        // printf(" %x",pre_cipher[i]);
    }

    Apply_Per_PI(pre_cipher, ciphertxt, InvIP, 8, 8);

}