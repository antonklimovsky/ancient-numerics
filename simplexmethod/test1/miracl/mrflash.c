/*
 *  MIRACL floating-Slash arithmetic
 *  mrflash.c
 *
 *  Copyright (c) 1988-1997 Shamus Software Ltd.
 */

#include <stdio.h>
#include <miracl.h>

#ifdef MR_FLASH

void flop(_MIPD_ flash x,flash y,int *op,flash z)
{ /* Do basic flash operation - depending on   *
   * op[]. Performs operations of the form

              A.f1(x,y) + B.f2(x,y)
              -------------------
              C.f3(x,y) + D.f4(x,y)

   * Four functions f(x,y) are supported and    *
   * coded thus                                 *
   *              00 - Nx.Ny                    *
   *              01 - Nx.Dy                    *
   *              10 - Dx.Ny                    *
   *              11 - Dx.Dy                    *
   * where Nx is numerator of x, Dx denominator *
   * of x, etc.                                 *
   * op[0] contains the codes for f1-f4 in last *
   * eight bits =  00000000f1f2f3f4             *
   * op[1], op[2], op[3], op[4] contain the     *
   * single precision multipliers A, B, C, D    */
    int i,code;
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(69)
    numer(_MIPP_ x,mr_mip->w1);
    denom(_MIPP_ x,mr_mip->w2);
    numer(_MIPP_ y,mr_mip->w3);
    denom(_MIPP_ y,mr_mip->w4);

    mr_mip->check=OFF;
    for (i=1;i<=4;i++)
    {
        zero(mr_mip->w0);
        if (op[i]==0) continue;
        code=(op[0]>>(2*(4-i)))&3;
        switch (code)
        {
        case 0 : if (x==y) multiply(_MIPP_ mr_mip->w1,mr_mip->w1,mr_mip->w0);
                 else      multiply(_MIPP_ mr_mip->w1,mr_mip->w3,mr_mip->w0);
                 break;
        case 1 : multiply(_MIPP_ mr_mip->w1,mr_mip->w4,mr_mip->w0);
                 break;
        case 2 : multiply(_MIPP_ mr_mip->w2,mr_mip->w3,mr_mip->w0);
                 break;
        case 3 : if(x==y) multiply(_MIPP_ mr_mip->w2,mr_mip->w2,mr_mip->w0);
                 else     multiply(_MIPP_ mr_mip->w2,mr_mip->w4,mr_mip->w0);
                 break;
        }
        premult(_MIPP_ mr_mip->w0,op[i],mr_mip->w0);
        switch (i)
        {
        case 1 : copy(mr_mip->w0,mr_mip->w5);
                 break;
        case 2 : add(_MIPP_ mr_mip->w5,mr_mip->w0,mr_mip->w5);
                 break;
        case 3 : copy(mr_mip->w0,mr_mip->w6);
                 break;
        case 4 : add(_MIPP_ mr_mip->w6,mr_mip->w0,mr_mip->w6);
                 break;
        }
    }
    mr_mip->check=ON;
    round(_MIPP_ mr_mip->w5,mr_mip->w6,z);
    MR_OUT
}

void fmul(_MIPD_ flash x,flash y,flash z)
{ /* Flash multiplication - z=x*y */
    int op[5];
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(35)
    op[0]=0x0C;
    op[1]=op[3]=1;
    op[2]=op[4]=0;
    flop(_MIPP_ x,y,op,z);
    MR_OUT
}

void fdiv(_MIPD_ flash x,flash y,flash z)
{ /* Flash divide - z=x/y */
    int op[5];
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(36)
    op[0]=0x48;
    op[1]=op[3]=1;
    op[2]=op[4]=0;
    flop(_MIPP_ x,y,op,z);
    MR_OUT
}

void fadd(_MIPD_ flash x,flash y,flash z)
{ /* Flash add - z=x+y */
    int op[5];
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(37)
    op[0]=0x6C;
    op[1]=op[2]=op[3]=1;
    op[4]=0;   
    flop(_MIPP_ x,y,op,z);
    MR_OUT
}

void fsub(_MIPD_ flash x,flash y,flash z)
{ /* Flash subtract - z=x-y */
    int op[5];
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(38)
    op[0]=0x6C;
    op[1]=op[3]=1;
    op[2]=(-1);
    op[4]=0;
    flop(_MIPP_ x,y,op,z);
    MR_OUT
}

int fcomp(_MIPD_ flash x,flash y)
{ /* compares two Flash numbers             *
   * returns -1 if y>x; +1 if x>y; 0 if x=y */
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return 0;

    MR_IN(39)
    numer(_MIPP_ x,mr_mip->w1);
    denom(_MIPP_ y,mr_mip->w2);
    mr_mip->check=OFF;
    multiply(_MIPP_ mr_mip->w1,mr_mip->w2,mr_mip->w5);
    numer(_MIPP_ y,mr_mip->w1);
    denom(_MIPP_ x,mr_mip->w2);
    multiply(_MIPP_ mr_mip->w1,mr_mip->w2,mr_mip->w0);
    mr_mip->check=ON;
    MR_OUT
    return (compare(mr_mip->w5,mr_mip->w0));
}

void ftrunc(_MIPD_ flash x,big y,flash z)
{ /* sets y=int(x), z=rem(x) - returns *
   * y only for ftrunc(x,y,y)       */
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(45)
    numer(_MIPP_ x,mr_mip->w1);
    denom(_MIPP_ x,mr_mip->w2);
    divide(_MIPP_ mr_mip->w1,mr_mip->w2,mr_mip->w3);
    copy(mr_mip->w3,y);
    if (y!=z) fpack(_MIPP_ mr_mip->w1,mr_mip->w2,z);
    MR_OUT
}

void fmodulo(_MIPD_ flash x,flash y,flash z)
{ /* sets z=x mod y */
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(89)
    fdiv(_MIPP_ x,y,mr_mip->w8);
    ftrunc(_MIPP_ mr_mip->w8,mr_mip->w8,mr_mip->w8);
    fmul(_MIPP_ mr_mip->w8,y,mr_mip->w8);
    fsub(_MIPP_ x,mr_mip->w8,z);
    MR_OUT
}

void fconv(_MIPD_ int n,int d,flash x)
{ /* convert simple fraction n/d to Flash form */
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(40)
    if (d<0)
    {
        d=(-d);
        n=(-n);
    }
    convert(_MIPP_ n,mr_mip->w5);
    convert(_MIPP_ d,mr_mip->w6);
    fpack(_MIPP_ mr_mip->w5,mr_mip->w6,x);
    MR_OUT
}

void frecip(_MIPD_ flash x,flash y)
{ /* set y=1/x */
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(41)

    numer(_MIPP_ x,mr_mip->w1);
    denom(_MIPP_ x,mr_mip->w2);
    fpack(_MIPP_ mr_mip->w2,mr_mip->w1,y);

    MR_OUT
}

void fpmul(_MIPD_ flash x,int n,int d,flash y)
{ /* multiply x by small fraction n/d - y=x*n/d */
    int r,g;
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;
    if (n==0 || size(x)==0)
    {
        zero(y);
        return;
    }
    if (n==d)
    {
        copy(x,y);
        return;
    }

    MR_IN(42)

    if (d<0)
    {
        d=(-d);
        n=(-n);
    }
    numer(_MIPP_ x,mr_mip->w1);
    denom(_MIPP_ x,mr_mip->w2);
    r=subdiv(_MIPP_ mr_mip->w1,d,mr_mip->w3);
    g=igcd(d,r);
    r=subdiv(_MIPP_ mr_mip->w2,n,mr_mip->w3);
    g*=igcd(n,r);
    mr_mip->check=OFF;
    premult(_MIPP_ mr_mip->w1,n,mr_mip->w5);
    premult(_MIPP_ mr_mip->w2,d,mr_mip->w6);
    subdiv(_MIPP_ mr_mip->w5,g,mr_mip->w5);
    subdiv(_MIPP_ mr_mip->w6,g,mr_mip->w6);
    mr_mip->check=ON;
    if (fit(mr_mip->w5,mr_mip->w6,mr_mip->nib))
        fpack(_MIPP_ mr_mip->w5,mr_mip->w6,y);
    else
        round(_MIPP_ mr_mip->w5,mr_mip->w6,y);
    MR_OUT
}

void fincr(_MIPD_ flash x,int n,int d,flash y)
{ /* increment x by small fraction n/d - y=x+(n/d) */
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(43)

    if (d<0)
    {
        d=(-d);
        n=(-n);
    }
    numer(_MIPP_ x,mr_mip->w1);
    denom(_MIPP_ x,mr_mip->w2);
    mr_mip->check=OFF;
    premult(_MIPP_ mr_mip->w1,d,mr_mip->w5);
    premult(_MIPP_ mr_mip->w2,d,mr_mip->w6);
    premult(_MIPP_ mr_mip->w2,n,mr_mip->w0);
    add(_MIPP_ mr_mip->w5,mr_mip->w0,mr_mip->w5);
    mr_mip->check=ON;
    if (d==1 && fit(mr_mip->w5,mr_mip->w6,mr_mip->nib))
        fpack(_MIPP_ mr_mip->w5,mr_mip->w6,y);
    else
        round(_MIPP_ mr_mip->w5,mr_mip->w6,y);
    MR_OUT
}

#endif

