/*
 *   MIRACL flash number builder: uses generator of
 *   regular continued fraction expansion to create
 *   a flash number, rounded if necessary.
 *   mrbuild.c
 *
 *   Copyright (c) 1988-1997 Shamus Software Ltd.
 */

#include <stdio.h>
#include <miracl.h>

#ifdef MR_FLASH

#define abs(x)  ((x)<0? (-(x)) : (x))

void build(_MIPD_ flash x,int (*gen)(_MIPT_ big,int))
{ /* Build x from its regular c.f. *
   * generated by gen()            */
    mr_small ex1,ex2,ex,st,sr;
    int a,b,c,d,rm,q,n,prc,lw2,lw4,lz;
    BOOL finoff,last;
    big t;
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    if (mr_mip->ERNUM) return;

    MR_IN(48)

    zero(mr_mip->w1);
    convert(_MIPP_ 1,mr_mip->w2);
    convert(_MIPP_ 1,mr_mip->w3);
    zero(mr_mip->w4);
    finoff=FALSE;
    last=FALSE;
    n=0;
    q=(*gen)(_MIPP_ x,n);   /* Note - first quotient may be zero */
    ex=mr_mip->base-1;
    if (mr_mip->nib==mr_mip->workprec) prc=mr_mip->nib;
    else prc=mr_mip->workprec+1;
    while (!mr_mip->ERNUM && q>=0)
    {
        if (q==MR_TOOBIG || n==0 || finoff)
        {
            if (q!=MR_TOOBIG) convert(_MIPP_ q,x);
            else last=FALSE;
            mr_mip->check=OFF;
            multiply(_MIPP_ mr_mip->w2,x,mr_mip->w0);
            subtract(_MIPP_ mr_mip->w1,mr_mip->w0,mr_mip->w7);
            mr_mip->check=ON;
            if ((int)(mr_mip->w7[0]&MR_OBITS)>mr_mip->nib) break;
            copy(mr_mip->w7,mr_mip->w1);
            t=mr_mip->w1,mr_mip->w1=mr_mip->w2,mr_mip->w2=t;   /* swap(w1,w2) */
            mr_mip->check=OFF;
            multiply(_MIPP_ mr_mip->w4,x,mr_mip->w0);
            subtract(_MIPP_ mr_mip->w3,mr_mip->w0,mr_mip->w7);
            mr_mip->check=ON;
            if ((int)(mr_mip->w7[0]&MR_OBITS)>mr_mip->nib)
            { /* oops! */
                fpack(_MIPP_ mr_mip->w1,mr_mip->w4,x);
                negate(x,x);
                mr_mip->EXACT=FALSE;
                MR_OUT
                return;
            }
            copy(mr_mip->w7,mr_mip->w3);
            t=mr_mip->w3,mr_mip->w3=mr_mip->w4,mr_mip->w4=t;   /* swap(w3,w4) */
            n++;
        }
        lw2=(int)(mr_mip->w2[0]&MR_OBITS);
        lw4=(int)(mr_mip->w4[0]&MR_OBITS);
        lz=lw2+lw4;
        if (lz > prc) break;  /* too big - exit */
        if (last)
        {
            if (finoff) break;
            finoff=TRUE;
            q=(*gen)(_MIPP_ x,n);
            continue;
        }
        if (lz>=prc-1)
        { /* nearly finished - so be careful not to overshoot */
            if (mr_mip->base==0)
            {
                st=mr_mip->w2[lw2]+1;
                if (st==0) ex1=1;
                else ex1=muldvm((mr_small)1,(mr_small)0,st,&sr);
                st=mr_mip->w4[lw4]+1;
                if (st==0) ex2=1;
                else ex2=muldvm((mr_small)1,(mr_small)0,st,&sr);
            }
            else
            {
                ex1=mr_mip->base/(mr_mip->w2[lw2]+1);
                ex2=mr_mip->base/(mr_mip->w4[lw4]+1);
            }
            if (ex2>ex1) ex=ex1,ex1=ex2,ex2=ex;
            if (lz==prc) ex=ex2;
            else         ex=ex1;
            last=TRUE;
        }
        a=1;
        b=0;
        c=0;
        d=1;
        forever
        {
            q=(*gen)(_MIPP_ x,n);          
            if (q<0 || q>=MR_TOOBIG/abs(d))
            { /* there could be more.... *** V3.21 mod *** */
                last=FALSE;
                break;
            }
            rm=b-q*d;
            b=d;
            d=rm;
            rm=a-q*c;
            a=c;
            c=rm;
            n++;
            if ((mr_small)(abs(c-d))>ex) break;
        }
        premult(_MIPP_ mr_mip->w1,c,mr_mip->w7);
        premult(_MIPP_ mr_mip->w1,a,mr_mip->w1);
        premult(_MIPP_ mr_mip->w2,b,mr_mip->w0);
        premult(_MIPP_ mr_mip->w2,d,mr_mip->w2);
        add(_MIPP_ mr_mip->w1,mr_mip->w0,mr_mip->w1);
        add(_MIPP_ mr_mip->w2,mr_mip->w7,mr_mip->w2);
        premult(_MIPP_ mr_mip->w3,c,mr_mip->w7);
        premult(_MIPP_ mr_mip->w3,a,mr_mip->w3);
        premult(_MIPP_ mr_mip->w4,b,mr_mip->w0);
        premult(_MIPP_ mr_mip->w4,d,mr_mip->w4);
        add(_MIPP_ mr_mip->w3,mr_mip->w0,mr_mip->w3);
        add(_MIPP_ mr_mip->w4,mr_mip->w7,mr_mip->w4);
    }
    if (fit(mr_mip->w2,mr_mip->w4,mr_mip->nib)) fpack(_MIPP_ mr_mip->w2,mr_mip->w4,x);
    else                            fpack(_MIPP_ mr_mip->w1,mr_mip->w3,x);
    negate (x,x);
    if (q!=(-1)) mr_mip->EXACT=FALSE;
    MR_OUT
}

#endif

