/*
 *   Proposed Digital Signature Standard (DSS)
 *
 *   See Communications ACM July 1992, Vol. 35 No. 7
 *   This new standard for digital signatures has been proposed by 
 *   the American National Institute of Standards and Technology (NIST)
 *   under advisement from the National Security Agency (NSA). 
 *
 *   This program asks for the name of a <file>, computes its message digest,
 *   signs it, and outputs the signature to a file <file>.dss. It is assumed 
 *   that the common values p, q and g, as well as the private key of the 
 *   signer have been previously generated by the dssgen program
 *
 *   Copyright (c) 1988-1997 Shamus Software Ltd.
 */

#include <stdio.h>
#include <miracl.h>
#include <stdlib.h>
#include <string.h>

void strip(char *name)
{ /* strip off filename extension */
    int i;
    for (i=0;name[i]!='\0';i++)
    {
        if (name[i]!='.') continue;
        name[i]='\0';
        break;
    }
}

static void hashing(FILE *fp,big hash)
{ /* compute hash function */
    char h[20];
    int i,ch;
    sha sh;
    shs_init(&sh);
    while ((ch=fgetc(fp))!=EOF) shs_process(&sh,ch);
    shs_hash(&sh,h);
    bytes_to_big(20,h,hash);
}

int main()
{
    FILE *fp;
    char ifname[13],ofname[13];
    big p,q,g,x,r,s,k,hash;
    long seed;
    int bits;
    miracl *mip;

/* get public data */
    fp=fopen("common.dss","r");
    if (fp==NULL)
    {
        printf("file common.dss does not exist\n");
        return 0;
    }

    fscanf(fp,"%d\n",&bits);
    mip=mirsys(3+bits/MIRACL,0);
    p=mirvar(0);
    q=mirvar(0);
    g=mirvar(0);
    x=mirvar(0);
    r=mirvar(0);
    s=mirvar(0);
    k=mirvar(0);
    hash=mirvar(0);

    mip->IOBASE=16;
    cinnum(p,fp);
    cinnum(q,fp);
    cinnum(g,fp);
    mip->IOBASE=10;
    fclose(fp);

/* randomise */
    printf("Enter 9 digit random number seed  = ");
    scanf("%ld",&seed);
    getchar();
    irand(seed);

/* calculate r - this can be done offline, 
   and hence amortized to almost nothing  */
    bigrand(q,k);
    powmod(g,k,p,r);   /* see brick.c for method to speed this up */
    divide(r,q,q);


/* get private key of signer */
    fp=fopen("private.dss","r");
    if (fp==NULL)
    {
        printf("file private.dss does not exist\n");
        return 0;
    }
    cinnum(x,fp);
    fclose(fp);

/* calculate message digest */
    printf("file to be signed = ");
    gets(ifname);
    strcpy(ofname,ifname);
    strip(ofname);
    strcat(ofname,".dss");
    if ((fp=fopen(ifname,"rb"))==NULL)
    {
        printf("Unable to open file %s\n",ifname);
        return 0;
    }
    hashing(fp,hash);
    fclose(fp);

/* calculate s */
    xgcd(k,q,k,k,k);
    mad(x,r,hash,q,q,s);
    mad(s,k,k,q,q,s);
    fp=fopen(ofname,"w");
    cotnum(r,fp);
    cotnum(s,fp);
    fclose(fp);
    mirexit();
    return 0;
}

