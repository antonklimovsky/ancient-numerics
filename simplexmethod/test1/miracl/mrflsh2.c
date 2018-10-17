/*
 *   MIRACL flash exponential and logs
 *   mrflsh2.c
 *
 *   Copyright (c) 1988-1995 Shamus Software Ltd.
 */

#include <stdio.h>
#include <math.h>
#include <miracl.h>

#ifdef MR_FLASH

#define mr_abs(x)  ((x)<0? (-(x)) : (x))

static int expon(_MIPD_ big w,int n)
{  /* generator for C.F. of e */ 
    if (n==0) return 2;
    if (n%3==2) return 2*(n/3)+2;
    else return 1;
}
 
void fexp(_MIPD_ flash x,flash y)
{ /* calculates y=exp(x) */
    int i,n,nsq,m,sqrn,op[5];
    BOOL minus,rem;
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;
    if (size(x)==0)
    {
        convert(_MIPP_ 1,y);
        return;
    }
    copy(x,y);

    MR_IN(54)

    minus=FALSE;
    if (size(y)<0)
    {
        minus=TRUE;
        negate(y,y);
    }
    ftrunc(_MIPP_ y,y,mr_mip->w9);
    n=size(y);
    if (n==MR_TOOBIG)
    {
        mr_berror(_MIPP_ MR_ERR_FLASH_OVERFLOW);
        MR_OUT
        return;
    }
    if (n==0) convert(_MIPP_ 1,y);
    else
    {
        build(_MIPP_ y,expon);
        if (minus)
        { /* underflow to zero - bit of a bodge */
            rem=mr_mip->ERCON;
            mr_mip->ERCON=TRUE;
            fpower(_MIPP_ y,n,y);
            mr_mip->ERCON=rem;
            if (mr_mip->ERNUM)
            {
                mr_mip->ERNUM=0;
                zero(y);
                MR_OUT
                return;
            }
        }
        else fpower(_MIPP_ y,n,y);
    }
    if (size(mr_mip->w9)==0)
    {
        if (minus) frecip(_MIPP_ y,y);
        MR_OUT
        return;
    }
    sqrn=isqrt(mr_mip->lg2b*mr_mip->workprec,mr_mip->lg2b);
    nsq=0;
    copy(mr_mip->w9,mr_mip->w8);
    frecip(_MIPP_ mr_mip->w9,mr_mip->w9);
    ftrunc(_MIPP_ mr_mip->w9,mr_mip->w9,mr_mip->w9);
    m=logb2(_MIPP_ mr_mip->w9);
    if (m<sqrn)
    { /* scale fraction down */
        nsq=sqrn-m;
        expint(_MIPP_ 2,nsq,mr_mip->w9);
        fdiv(_MIPP_ mr_mip->w8,mr_mip->w9,mr_mip->w8);
    }
    zero(mr_mip->w10);
    op[0]=0x4B;  /* set up for x/(C+y) */
    op[1]=1;
    op[2]=0;
    for (m=sqrn;m>0;m--)
    { /* Unwind C.F. expansion for exp(x)-1 */
        if (m%2==0) op[4]=2,op[3]=1;
        else        op[4]=m,op[3]=(-1);
        flop(_MIPP_ mr_mip->w8,mr_mip->w10,op,mr_mip->w10);
    }
    op[0]=0x2C;  /* set up for (x+2).y */
    op[1]=op[3]=1;
    op[2]=2;
    op[4]=0;
    for (i=0;i<nsq;i++)
    { /* now square it back up again */
        flop(_MIPP_ mr_mip->w10,mr_mip->w10,op,mr_mip->w10);
    }
    op[2]=1;
    flop(_MIPP_ mr_mip->w10,y,op,y);
    if (minus) frecip(_MIPP_ y,y);
    MR_OUT
}

void flog(_MIPD_ flash x,flash y)
{ /* calculate y=log(x) to base e */
    BOOL hack;
    int op[5];
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    copy(x,y);
    if (mr_mip->ERNUM) return;
    if (size(y)==1)
    {
        zero(y);
        return;
    }

    MR_IN(55)

    if (size(y)<=0)
    {
        mr_berror(_MIPP_ MR_ERR_NEG_LOG);
        MR_OUT
        return;
    }
    hack=FALSE;
    if (mr_lent(y)<=2)
    { /* for 'simple' values of y */
        hack=TRUE;
        build(_MIPP_ mr_mip->w11,expon);
        fdiv(_MIPP_ y,mr_mip->w11,y);
    }
    op[0]=0x68;
    op[1]=op[3]=1;
    op[2]=(-1);
    op[4]=0;
    mr_mip->workprec=mr_mip->stprec;
    dconv(_MIPP_ log(fdsize(_MIPP_ y)),mr_mip->w11);
    while (mr_mip->workprec!=mr_mip->nib)
    { /* Newtons iteration w11=w11+(y-exp(w11))/exp(w11) */
        if (mr_mip->workprec<mr_mip->nib) mr_mip->workprec*=2;
        if (mr_mip->workprec>=mr_mip->nib) mr_mip->workprec=mr_mip->nib;
        else if (mr_mip->workprec*2>mr_mip->nib) mr_mip->workprec=(mr_mip->nib+1)/2;
        fexp(_MIPP_ mr_mip->w11,mr_mip->w12);
        flop(_MIPP_ y,mr_mip->w12,op,mr_mip->w12);
        fadd(_MIPP_ mr_mip->w12,mr_mip->w11,mr_mip->w11);
    }
    copy(mr_mip->w11,y);
    if (hack) fincr(_MIPP_ y,1,1,y);
    MR_OUT
}
    
void fpowf(_MIPD_ flash x,flash y,flash z)
{ /* raise flash number to flash power *
   *     z=x^y  -> z=exp(y.log(x))     */
    int n;
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(56)

    n=size(y);
    if (mr_abs(n)<MR_TOOBIG)
    { /* y is small int */
        fpower(_MIPP_ x,n,z);
        MR_OUT
        return;
    }
    frecip(_MIPP_ y,mr_mip->w13);
    n=size(mr_mip->w13);
    if (mr_abs(n)<MR_TOOBIG)
    { /* 1/y is small int */
        froot(_MIPP_ x,n,z);
        MR_OUT
        return;
    }
    copy(x,z);
    flog(_MIPP_ z,z);
    fdiv(_MIPP_ z,mr_mip->w13,z);
    fexp(_MIPP_ z,z);
    MR_OUT
}

#endif

