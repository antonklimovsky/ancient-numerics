/*
 *   MIRACL Comba's method for ultimate speed modular multiplication
 *   mrcomba.tpl 
 *
 *   See "Exponentiation Cryptosystems on the IBM PC", IBM Systems
 *   Journal Vol. 29 No. 4 1990. Comba's method has been extended to 
 *   implement Montgomery reduction. 
 *
 *   Here the inner loops of the basic multiplication, squaring and 
 *   Montgomery's redc() functions are completely unravelled, and 
 *   reorganised for maximum possible speed. 
 *
 *   This approach is recommended for maximum speed where parameters
 *   are fixed and compute resources are constrained. The processor must 
 *   support an unsigned multiply instruction, and should have a carry flag.
 *
 *   This file is a template. To fill in the gaps and create mrcomba.c, 
 *   you must run the mex.c program to insert the C or assembly language 
 *   macros from the appropriate .mcs file. For use with C MR_NOASM must
 *   be defined in mirdef.h
 *
 *   This method would appear to be particularly useful for implementing 
 *   fast Elliptic Curve Cryptosystems over GF(p) and fast 1024-bit RSA
 *   decryption.
 *
 *   The #define MR_COMBA in mirdef.h determines the FIXED size of 
 *   modulus to be used. This *must* be determined at compile time. 
 *
 *   Note that this module can generate a *lot* of code for large values 
 *   of MR_COMBA. This should have a maximum value of 8-16. Any larger 
 *   that and you should define MR_KCM instead - see mrkcm.tpl
 *
 *   Note that on some processors it is *VITAL* that arrays be aligned on 
 *   4-byte boundaries
 *
 *   Copyright (c) 1988-2001 Shamus Software Ltd.
 */

#include <stdio.h>     
#include <miracl.h>    

#ifdef MR_COMBA
#if INLINE_ASM == 1    
#define N 2
#define POINTER WORD PTR  
#define PBP bp   
#define PBX bx   
#define PSI si   
#define PDI di   
#define DSI si   
#define DDI di   
#define DBP bp   
#define DAX ax   
#define DCX cx   
#define DDX dx   
#endif   
 
#if INLINE_ASM == 2    
#define N 4
#define POINTER DWORD PTR   
#define PBP bp   
#define PBX bx   
#define PSI si   
#define PDI di   
#define DSI esi  
#define DDI edi  
#define DBP ebp  
#define DAX eax  
#define DCX ecx  
#define DDX edx  
#endif           
  
#if INLINE_ASM == 3    
#define N 4
#define POINTER DWORD PTR   
#define PBP ebp   
#define PBX ebx   
#define PSI esi   
#define PDI edi   
#define DSI esi  
#define DDI edi  
#define DBP ebp  
#define DAX eax  
#define DCX ecx  
#define DDX edx  
#endif           
  
/* NOTE! z must be distinct from x and y */

void comba_mult(_MIPD_ big x,big y,big z) 
{ /* comba multiplier */
    int i;
    big a,b,c,w;
#ifdef MR_NOASM 
    mr_small extra;
    mr_large pp,sum;
#endif
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    w=mr_mip->w0; 
    for (i=2*MR_COMBA+1;i<=(w[0]&MR_OBITS);i++) w[i]=0;
  /*  zero(w); */
    w[0]=2*MR_COMBA;
    a=&x[1]; b=&y[1]; c=&w[1];
/*** MULTIPLY ***/      /* multiply a by b, result in c */

    if (w!=z) copy (w,z); 
}   
 
/* NOTE! z and x must be distinct */

void comba_square(_MIPD_ big x,big z)  
{ /* super comba squarer */
    int i;
    big a,c,w;
#ifdef MR_NOASM
    mr_small extra;
    mr_large pp,sum;
#endif
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    w=mr_mip->w0;       
    for (i=2*MR_COMBA+1;i<=(w[0]&MR_OBITS);i++) w[i]=0;  
  /*  zero(w); */
    w[0]=2*MR_COMBA;
    a=&x[1]; c=&w[1];
/*** SQUARE ***/    /* squares a, result in b */

    if (w!=z) copy (w,z);      
}                        
                         
/* NOTE! t and z must be distinct! */

void comba_redc(_MIPD_ big t,big z)     
{  /* super comba Montgomery redc() function */                      
    mr_small ndash,carry;
#ifdef MR_NOASM
    mr_small sp,extra;
    mr_large pp,sum,u;
#endif
    unsigned int i;
    big w,modulus;
    big a,b;
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
#ifdef MR_SPECIAL

/* !!! Implement here a "special" fast method for modular reduction, for
   for a particular modulus. Implemented here for 2#192-2#64-1       
   See for example "Software Implementation of the NIST Elliptic
   Curves Over Prime Fields", Brown et al., Report #36, 2000 available
   from www.cacr.math.uwaterloo.ca */
    
    int overshoot;
    mr_small k[6];
    big c;
    modulus=mr_mip->modulus;     
    for (i=MR_COMBA+1;i<=(z[0]&MR_OBITS);i++) z[i]=0;
  /*  zero(z);   */
    z[0]=MR_COMBA;
    a=&t[1]; b=k; c=&z[1];
    k[0]=a[6]; k[1]=a[7]; k[2]=a[6]; k[3]=a[7]; k[4]=k[5]=0; 
    
/*** ADDITION ***/
    overshoot=carry;  
    a=c;  c=&t[1];
    k[0]=k[1]=0; k[2]=c[8]; k[3]=c[9]; k[4]=c[8]; k[5]=c[9];

/*** INCREMENT ***/
    overshoot+=carry;
    k[0]=c[10]; k[1]=c[11]; k[2]=c[10]; k[3]=c[11]; k[4]=c[10]; k[5]=c[11];
    
/*** INCREMENT ***/
    overshoot+=carry;
    b=&modulus[1];
    while(overshoot!=0)
    {
/*** DECREMENT ***/
        overshoot-=carry;
    }
    if (compare(z,modulus)>=0)
    {
/*** DECREMENT ***/
    }
    mr_lzero(z);

#else
    modulus=mr_mip->modulus;  
    ndash=mr_mip->ndash;
    w=mr_mip->w0;
    if (t!=w) copy(t,w);       
    w[0]=2*MR_COMBA+1;
    a=&w[1]; b=&modulus[1];

/*** REDC ***/      /* reduces a mod b */
    
    for (i=MR_COMBA+2;i<=(z[0]&MR_OBITS);i++) z[i]=0;
 /*   zero(z); */
    z[0]=MR_COMBA+1;
    for (i=1;i<=MR_COMBA+1;i++) z[i]=w[i+MR_COMBA];
    mr_lzero(z);
    if (compare(z,modulus)>=0) 
    {
        a=&z[1]; b=&modulus[1];
/*** DECREMENT ***/    
        z[MR_COMBA+1]=0;
        mr_lzero(z);
    }
#endif
} 

void comba_add(_MIPD_ big x,big y,big w)
{ /* fast modular addition */
    int i;
    big modulus;
    big a,b,c;
    mr_small carry;  
#ifdef MR_NOASM
    mr_large u;
#endif
#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    modulus=mr_mip->modulus;
    if (w!=x && w!=y) 
    {
        for (i=MR_COMBA+1;i<=(w[0]&MR_OBITS);i++) w[i]=0;
        /* zero(w); */
    }
    
    a=&x[1]; b=&y[1]; c=&w[1];
/*** ADDITION ***/        /* add a and b, result in c */
    w[0]=MR_COMBA;
    if (carry || compare(w,modulus)>=0)
    {

        a=&w[1]; b=&modulus[1];
/*** DECREMENT ***/        /* decrement b from a */
    
    }
 
    mr_lzero(w);   
}

void comba_sub(_MIPD_ big x,big y,big w)
{ /* fast modular subtraction */
    int i;
    big modulus;
    big a,b,c;
    mr_small carry;  
#ifdef MR_NOASM
    mr_large u;
#endif

#ifndef MR_GENERIC_MT
    miracl *mr_mip=get_mip();
#endif
    modulus=mr_mip->modulus;
    if (x!=w && y!=w) 
    {
        for (i=MR_COMBA+1;i<=(w[0]&MR_OBITS);i++) w[i]=0;   
        /* zero(w); */
    }

    a=&x[1]; b=&y[1]; c=&w[1];
/*** SUBTRACTION ***/

    if (carry)
    {
        a=&w[1]; b=&modulus[1]; 
/*** INCREMENT ***/        /* add a and b, result in c */
    
    }
    w[0]=MR_COMBA;
    mr_lzero(w); 
}

#endif


