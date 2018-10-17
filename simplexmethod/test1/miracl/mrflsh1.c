/*
 *   MIRACL flash roots and powers
 *   mrflsh1.c
 *
 *   Copyright (c) 1988-1997 Shamus Software Ltd.
 */

#include <stdio.h>
#include <math.h>
#include <miracl.h>

#ifdef MR_FLASH

static int quad(_MIPD_ big w,int n)
{ /* generator for C.F. of square root of small integer */
    int t;
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (n==0)
    {
        mr_mip->oldn=(-1);
        mr_mip->b=2*mr_mip->RD;
        mr_mip->c=mr_mip->b;
        mr_mip->a=1;
        mr_mip->d=mr_mip->RS;
        mr_mip->r=mr_mip->RD;
        if (mr_mip->r>=MR_TOOBIG)
        {
            convert(_MIPP_ mr_mip->r,w);
            return MR_TOOBIG;
        }
        return (mr_mip->r);
    }
    else if (n==mr_mip->oldn) return (mr_mip->r);
    t=mr_mip->a;
    mr_mip->a=mr_mip->r*(mr_mip->c-mr_mip->b)+mr_mip->d;
    mr_mip->d=t;
    mr_mip->r=mr_mip->b/mr_mip->a;
    mr_mip->c=mr_mip->b;
    mr_mip->b=2*mr_mip->RD-mr_mip->b%mr_mip->a;
    mr_mip->oldn=n;
    if (mr_mip->r>=MR_TOOBIG)
    {
        convert(_MIPP_ mr_mip->r,w);
        return (MR_TOOBIG);
    }
    return mr_mip->r;
}

void fpower(_MIPD_ flash x,int n,flash w)
{ /* raise floating-slash number to integer power w=x^n */
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    copy(x,mr_mip->w8);
    zero(w);
    if (mr_mip->ERNUM || size(mr_mip->w8)==0) return;
    convert(_MIPP_ 1,w);
    if (n==0) return;

    MR_IN(51)

    if (n<0)
    {
        n=(-n);
        frecip(_MIPP_ mr_mip->w8,mr_mip->w8);
    }
    if (n==1)
    {
        copy(mr_mip->w8,w);
        MR_OUT
        return;
    }
    forever
    {
        if (n%2!=0) fmul(_MIPP_ w,mr_mip->w8,w);
        n/=2;
        if (mr_mip->ERNUM || n==0) break;
        fmul(_MIPP_ mr_mip->w8,mr_mip->w8,mr_mip->w8);
    }
    MR_OUT
}
 
BOOL froot(_MIPD_ flash x,int n,flash w)
{ /* extract nth root of x  - w=x^(1/n) using Newtons method */
    BOOL minus,rn,rm,hack;
    int nm,dn,s,op[5];
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    copy(x,w);
    if (mr_mip->ERNUM || n==1) return TRUE;
    if (n==(-1))
    {
        frecip(_MIPP_ w,w);
        return TRUE;
    }

    MR_IN(52)

    minus=FALSE;
    if (n<0)
    {
        minus=TRUE;
        n=(-n);
    }
    s=exsign(w);
    if (n%2==0 && s==MINUS)
    {
        mr_berror(_MIPP_ MR_ERR_NEG_ROOT);
        MR_OUT
        return FALSE;
    }
    insign(PLUS,w);
    numer(_MIPP_ w,mr_mip->w8);
    denom(_MIPP_ w,mr_mip->w9);
    rn=nroot(_MIPP_ mr_mip->w8,n,mr_mip->w8);
    rm=nroot(_MIPP_ mr_mip->w9,n,mr_mip->w9);
    if (rn && rm)
    {
        fpack(_MIPP_ mr_mip->w8,mr_mip->w9,w);
        if (minus) frecip(_MIPP_ w,w);
        insign(s,w);
        MR_OUT
        return TRUE;
    }
    nm=size(mr_mip->w8);
    dn=size(mr_mip->w9);
    if (n==2 && ((nm<MR_TOOBIG) || rn) && ((dn<MR_TOOBIG) || rm))
    {
        if (!rn && nm<MR_TOOBIG)
        {
            multiply(_MIPP_ mr_mip->w8,mr_mip->w8,mr_mip->w8);
            numer(_MIPP_ w,mr_mip->w7);
            subtract(_MIPP_ mr_mip->w7,mr_mip->w8,mr_mip->w8);
            mr_mip->RS=(int)(mr_mip->w8[1]+mr_mip->base*mr_mip->w8[2]);
            mr_mip->RD=nm;
            build(_MIPP_ mr_mip->w8,quad);
        }
        if (!rm && dn<MR_TOOBIG)
        {
            multiply(_MIPP_ mr_mip->w9,mr_mip->w9,mr_mip->w9);
            denom(_MIPP_ w,mr_mip->w7);
            subtract(_MIPP_ mr_mip->w7,mr_mip->w9,mr_mip->w9);
            mr_mip->RS=(int)(mr_mip->w9[1]+mr_mip->base*mr_mip->w9[2]);
            mr_mip->RD=dn;
            build(_MIPP_ mr_mip->w9,quad);
        }
        if (size(mr_mip->w9)==1) copy(mr_mip->w8,w);
        else             fdiv(_MIPP_ mr_mip->w8,mr_mip->w9,w);
        if (minus) frecip(_MIPP_ w,w);
        insign(s,w);
        MR_OUT
        return FALSE;
    }
    hack=FALSE;
    if (mr_lent(w)<=2)
    { /* for 'simple' w only */
        hack=TRUE;
        fpi(_MIPP_ mr_mip->pi);
        fpmul(_MIPP_ mr_mip->pi,1,3,mr_mip->w10);
        fpower(_MIPP_ mr_mip->w10,n,mr_mip->w10);
        fmul(_MIPP_ w,mr_mip->w10,w);
    }
    op[0]=0x6C;  /* set up for [(n-1).x+y]/n */
    op[1]=n-1;
    op[2]=1;
    op[3]=n;
    op[4]=0;
    mr_mip->workprec=mr_mip->stprec;
    dconv(_MIPP_ pow(fdsize(_MIPP_ w),1.0/(double)n),mr_mip->w10);
    while (mr_mip->workprec!=mr_mip->nib)
    { /* Newtons iteration w10=(w/w10^(n-1)+(n-1)*w10)/n */
        if (mr_mip->workprec<mr_mip->nib) mr_mip->workprec*=2;
        if (mr_mip->workprec>=mr_mip->nib) mr_mip->workprec=mr_mip->nib;
        else if (mr_mip->workprec*2>mr_mip->nib) mr_mip->workprec=(mr_mip->nib+1)/2;
        fpower(_MIPP_ mr_mip->w10,n-1,mr_mip->w9);
        fdiv(_MIPP_ w,mr_mip->w9,mr_mip->w9);
        flop(_MIPP_ mr_mip->w10,mr_mip->w9,op,mr_mip->w10);
    }
    copy(mr_mip->w10,w);
    op[0]=0x48;
    op[1]=3;
    op[3]=1;
    op[2]=op[4]=0;
    if (hack) flop(_MIPP_ w,mr_mip->pi,op,w);
    if (minus) frecip(_MIPP_ w,w);
    insign(s,w);
    MR_OUT
    return FALSE;
}

#endif

