/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/21 10:56:58 $
 * $Revision: 1.1.2.9 $
 * $Id: zigzag.cpp,v 1.1.2.9 2006/07/21 10:56:58 syed Exp $
 ********************************************************************/

#include<iostream.h>
#include<fstream.h>
#include<stdio.h>
#include "coding/zigzag.h"
#include "tools/all.h"


//static int zz_id=0;

//#define ZIGZAG_DEBUG
//#define PROFILING_CHECK

WVector WZigZag::encode(WVector &data, int blocksize) {
    if (data.vect ==NULL || data.taille<1)
        throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if(blocksize<1)
        throw CZigZagException(INVALID_ZIGZAG_PARAMETER, __LINE__, __FILE__);
    if(initialized==0)
        throw CZigZagException(NOT_INITIALIZED, __LINE__, __FILE__);
    block_size=blocksize;
    I = blocksize/J+(blocksize%J==0?0:1);
    padsize = J*I-blocksize;
    int nbblock = data.taille/blocksize+(data.taille%blocksize==0?0:1);
    int outsize = I*H*nbblock, par=0, t;
//    superpadsize = blocksize-data.taille%blocksize;
    superpadsize = nbblock*blocksize-data.taille;
    WVector out(outsize,8), data1(nbblock*blocksize,padding_value);
    data1.insert(data,0);
    parity = WVector(outsize-data.taille);


    for(int i=0,k=0,t=0;i<nbblock;i++){
        par = initial_parity;
        for(int j=0;j<blocksize;j++) {
            if (j>0 && (j%J)==0) {
                parity.vect[t++]=par = ((par==parity_rule)?0:1);
                if(concat_mode==0) {
                    out.vect[k++]= par==0?1:-1;
                } else {
                    out.vect[k++]= par==0?0:1;
                }
            }
            par ^=data1.vect[i*blocksize+j];
            if(concat_mode==0) {
                out.vect[k++]=data1.vect[i*blocksize+j]==0?1:-1;
            } else {
                out.vect[k++]=data1.vect[i*blocksize+j]==0?0:1;
            }
        }

        for(int j=0; j<padsize;j++) {
            par ^= padding_value;
            if(concat_mode==0) {
                out.vect[k++] = padding_value==0?1:-1;
            } else {
                out.vect[k++] = padding_value==0?0:1;
            }
        }
        parity.vect[t++]=lastparity = par^((par==parity_rule)?0:1);
        if(concat_mode==0) {
            out.vect[k++]=(par==parity_rule)?1:-1;
        } else {
            out.vect[k++]=(par==parity_rule)?0:1;
        }
    }
    return out;
}

DVector WZigZag::decode(DVector &data){
    if (data.vect ==NULL || data.taille<=0) throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if(block_size<1) throw CZigZagException(INVALID_ZIGZAG_PARAMETER, __LINE__, __FILE__);
    if(initialized==0) throw CZigZagException(NOT_INITIALIZED, __LINE__, __FILE__);
    // end sanity check

    int blocksize=block_size;
    int codedsize = I*H;

    int nbblock = data.taille/codedsize+(data.taille%codedsize==0?0:1);
    int t;
    //  superpadsize = nbblock*blocksize-data.taille;

    double min, dv;

    DVector min1(I,M_HUGE_NUM), min2(I,M_HUGE_NUM), fv(I+1), bv(I+1);
    WVector min1_idx(I,0), min2_idx(I,0),
            sign(I,1), fsign(I,1),bsign(I,1);
    int dsign = 1;

    // data in = a very long sequemce
    // need to process block by block of length blocksize
    // in each block, need to process segment per segment of segment_length

    DVector segment, vLd(blocksize*nbblock);
    for(int k = 0; k<nbblock;k++) { // processing block by block of codedsize ..
        // copy block in matrix IxH representation
        DMatrix subblock = ReshapeLinewise(data.copy(k*codedsize,(k+1)*codedsize), I, H);
        //initialize MLA table
        DMatrix Ld(I,J); DVector Fp(I+1), Bp(I+1);
        Ld.reset();        Fp.reset();        Bp.reset();
        Fp.vect[0]=M_HUGE_NUM;
        Bp.vect[I]=subblock.mat[I-1][H-1];
        // processing segment by segment
        for(int ip=1,i=0;i<I;i++,ip++) { // get W(.) without F(p) and B(p)
            // calculate W(d(i,1),...,d(i,J)), identify two mins in the segment..
            sign.vect[i] = 1;
            min1.vect[i] = min2.vect[i] = M_HUGE_NUM; // arbitrary initialize a min
            min1_idx.vect[i] = min2_idx.vect[i] = -1; // identify the min position

            for(int j = 0; j < J;j++) { // look for another min and calculate W(.)
                sign.vect[i] *= (dv = subblock.mat[i][j]>=0.0?1:-1); // get the sign of current data
                dv *= subblock.mat[i][j]; // get the absolute value of segment.vect[j]
                if(dv<=min1.vect[i]) { // current value is min
                    min2.vect[i]=min1.vect[i]; min1.vect[i]=dv; // update min table
                    min2_idx.vect[i] = min1_idx.vect[i]; min1_idx.vect[i] = j;
                } else if (dv <=min2.vect[i]) { // second min value
                    min2.vect[i] = dv;            // update min table
                    min2_idx.vect[i] = j;
                }
            }
        }
        for(int ip=1,i=0;i<I;i++,ip++) { // get forward MLA
            // calculate F(p(i)) based on W(d(i,1),...,d(i,J)) calculated previously
            fsign.vect[i] = sign.vect[i]*(fv.vect[ip-1]=(Fp.vect[ip-1]>=0?1:-1)); // calculate the sign ..
            fv.vect[ip-1] *= Fp.vect[ip-1]; // fabs(F(p(i-1)))
            Fp.vect[ip] = subblock.mat[i][J]; // F(p(i)=p(i) {+W(F(p(i-1)),d(i,1),...,d(i,J))}
            if (fv.vect[ip-1]<=min1.vect[i]) { // if (F(p(i-1)) <= min1 => W(.) = sign*fabs(F(p(i-1))
                Fp.vect[ip] += fsign.vect[i]*fv.vect[ip-1];
            } else { // if (F(p(i-1)) > min1 => W(.) = sign*min1
                Fp.vect[ip] += fsign.vect[i]*min1.vect[i];
            }
        }
        for(int ip=I,i=I-1;i>=0;i--,ip--) { // get backward MLA
            // calculate B(p(i-1)) based on W(d(i,1),...,d(i,J)) calculated previously
            bsign.vect[i] = sign.vect[i]*(bv.vect[ip]=(Bp.vect[ip]>=0?1:-1)); // calculate the sign
            bv.vect[ip] *= Bp.vect[ip]; // fabs(B(p(i)))
            Bp.vect[ip-1] = i ? subblock.mat[i-1][J] :M_HUGE_NUM ; //  B(p(i-1)) = p(i-1) {+W(d(i,1),...,d(i,J),B(p(i)))}
            if(bv.vect[ip]<=min1.vect[i]) { // if B(p(i))<=min1 => W(.) = sign*fabs(B(p(i))
                Bp.vect[ip-1] += bsign.vect[i]*bv.vect[ip];
            } else { // if B(p(i)) > min1 => W(.) = sign*min1
                Bp.vect[ip-1] += bsign.vect[i]*min1.vect[i];
            }
        }
        for(int ip=1,i=0;i<I;i++,ip++) { // calculate L(d(i,j))
            int fpsign = (Fp.vect[ip-1]>=0?1:-1), bpsign = (Bp.vect[ip]>=0?1:-1);
            for (int j=0;j<J;j++) {
                // L(d(i,j)) =d(i,j) {+W(F(p(i-1)),d(i,1),...,d(i,j-1),d(i,j+1),...,d(i,J), B(p(i))}
                Ld.mat[i][j]=subblock.mat[i][j];

                // remove d(i,j) from the sign
                dsign = sign.vect[i]*(subblock.mat[i][j]>=0?1:-1);
                dsign *= fpsign*bpsign; //calculate new sign

                min = min1.vect[i];

                if(j==min1_idx.vect[i]) //check if d(i,j) is the current min ?
                    min = min2.vect[i];

                if(min>fv.vect[ip-1]) // check if F(p(i-1)) is the min
                    min = fv.vect[ip-1];

                if(min>bv.vect[ip])   // check if B(p(i)) is the min
                    min = bv.vect[ip];

                if(decoder_output_type==OUTPUT_TYPE_IS_BIT_LLR) {
                    Ld.mat[i][j] = min*dsign;
                } else {
                    Ld.mat[i][j] += min*dsign;
                }
            }
        }
        vLd.insert(MatrixAsLinePackedVector(Ld), k*blocksize);
    }
    return vLd.copy(0,vLd.taille-superpadsize);
}

// second variant of decoder: using llr input as well

DVector WZigZag::decode(DVector &data, DVector &llr){
    if (data.vect ==NULL || data.taille<=0)
        throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if(block_size<1)
        throw CZigZagException(INVALID_ZIGZAG_PARAMETER, __LINE__, __FILE__);
    if(initialized==0)
        throw CZigZagException(NOT_INITIALIZED, __LINE__, __FILE__);
    if(llr.vect==NULL || llr.taille<=0) {
        llr=DVector(data.taille);
    }else if(llr.taille!= data.taille) {
        llr=DVector(data.taille);
    }
    // end sanity check

    int blocksize=block_size;
    int codedsize = I*H;

    int nbblock = data.taille/codedsize+(data.taille%codedsize==0?0:1);
    int t;

    double min, dv;

    DVector min1(I,M_HUGE_NUM), min2(I,M_HUGE_NUM), fv(I+1), bv(I+1);
    WVector min1_idx(I,0), min2_idx(I,0),
            sign(I,1), fsign(I,1),bsign(I,1);
    int dsign = 1;

    // data in = a very long sequence
    // need to process block by block of length blocksize
    // in each block, need to process segment per segment of segment_length

    DVector segment, vLd(blocksize*nbblock), vllr(blocksize*nbblock);
    for(int k = 0; k<nbblock;k++) { // processing block by block of codedsize ..
        // copy block in matrix IxH representation
        DMatrix subblock = ReshapeLinewise(data.copy(k*codedsize,(k+1)*codedsize), I, H);
        DMatrix subblock_llr = ReshapeLinewise(llr.copy(k*codedsize,(k+1)*codedsize), I, H);
        //initialize MLA table
        DMatrix Ld(I,J); DVector Fp(I+1), Bp(I+1);
        Ld.reset();        Fp.reset();        Bp.reset();
        Fp.vect[0]=M_HUGE_NUM;
        Bp.vect[I]=subblock.mat[I-1][H-1];
        // processing segment by segment
        for(int ip=1,i=0;i<I;i++,ip++) { // get W(.) without F(p) and B(p)
            // calculate W(d(i,1),...,d(i,J)), identify two mins in the segment..
            sign.vect[i] = 1;
            min1.vect[i] = min2.vect[i] = M_HUGE_NUM; // arbitrary initialize a min
            min1_idx.vect[i] = min2_idx.vect[i] = -1; // identify the min position

            for(int j = 0; j < J;j++) { // look for another min and calculate W(.)
                sign.vect[i] *= (dv = subblock.mat[i][j]>=0.0?1:-1); // get the sign of current data
                dv *= subblock.mat[i][j]; // get the absolute value of segment.vect[j]
                if(dv<=min1.vect[i]) { // current value is min
                    min2.vect[i]=min1.vect[i]; min1.vect[i]=dv; // update min table
                    min2_idx.vect[i] = min1_idx.vect[i]; min1_idx.vect[i] = j;
                } else if (dv <=min2.vect[i]) { // second min value
                    min2.vect[i] = dv;            // update min table
                    min2_idx.vect[i] = j;
                }
            }
        }
        for(int ip=1,i=0;i<I;i++,ip++) { // get forward MLA
            // calculate F(p(i)) based on W(d(i,1),...,d(i,J)) calculated previously
            fsign.vect[i] = sign.vect[i]*(fv.vect[ip-1]=(Fp.vect[ip-1]>=0?1:-1)); // calculate the sign ..
            fv.vect[ip-1] *= Fp.vect[ip-1]; // fabs(F(p(i-1)))
            Fp.vect[ip] = subblock.mat[i][J]; // F(p(i)=p(i) {+W(F(p(i-1)),d(i,1),...,d(i,J))}
            if (fv.vect[ip-1]<=min1.vect[i]) { // if (F(p(i-1)) <= min1 => W(.) = sign*fabs(F(p(i-1))
                Fp.vect[ip] += fsign.vect[i]*fv.vect[ip-1];
            } else { // if (F(p(i-1)) > min1 => W(.) = sign*min1
                Fp.vect[ip] += fsign.vect[i]*min1.vect[i];
            }
        }
        for(int ip=I,i=I-1;i>=0;i--,ip--) { // get backward MLA
            // calculate B(p(i-1)) based on W(d(i,1),...,d(i,J)) calculated previously
            bsign.vect[i] = sign.vect[i]*(bv.vect[ip]=(Bp.vect[ip]>=0?1:-1)); // calculate the sign
            bv.vect[ip] *= Bp.vect[ip]; // fabs(B(p(i)))
            Bp.vect[ip-1] = i ? subblock.mat[i-1][J] :M_HUGE_NUM ; //  B(p(i-1)) = p(i-1) {+W(d(i,1),...,d(i,J),B(p(i)))}
            if(bv.vect[ip]<=min1.vect[i]) { // if B(p(i))<=min1 => W(.) = sign*fabs(B(p(i))
                Bp.vect[ip-1] += bsign.vect[i]*bv.vect[ip];
            } else { // if B(p(i)) > min1 => W(.) = sign*min1
                Bp.vect[ip-1] += bsign.vect[i]*min1.vect[i];
            }
        }
        for(int ip=1,i=0;i<I;i++,ip++) { // calculate L(d(i,j))
            int fpsign = (Fp.vect[ip-1]>=0?1:-1), bpsign = (Bp.vect[ip]>=0?1:-1);
            for (int j=0;j<J;j++) {
                // L(d(i,j)) =d(i,j) {+W(F(p(i-1)),d(i,1),...,d(i,j-1),d(i,j+1),...,d(i,J), B(p(i))}
                Ld.mat[i][j]=subblock.mat[i][j];

                dsign = sign.vect[i]*(subblock.mat[i][j]>=0?1:-1); // remove d(i,j) from the sign
                dsign *= fpsign*bpsign; //calculate new sign

                min = min1.vect[i];


                if (j==min1_idx.vect[i])//check if d(i,j) is the current min ?
                    min = min2.vect[i];
                if (min>fv.vect[ip-1])  // check if F(p(i-1)) is the min
                    min = fv.vect[ip-1];
                if (min>bv.vect[ip])    // check if B(p(i)) is the min
                     min = bv.vect[ip];
                Ld.mat[i][j] += subblock_llr.mat[i][j] + min*dsign;
            }
        }
        vLd.insert(MatrixAsLinePackedVector(Ld), k*blocksize);
        vllr.insert(MatrixAsLinePackedVector(subblock_llr), k*blocksize);
    }

    llr = vllr.copy(0, vllr.taille-superpadsize);
    return vLd.copy(0,vLd.taille-superpadsize);
}

#define TAKSIAP
#ifdef TAKSIAP

DVector WConcatZigZag::decode(DVector &data, int iter_max){
    if (data.vect ==NULL || data.taille<=0)
        throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if(block_size<1)
        throw CZigZagException(INVALID_ZIGZAG_PARAMETER, __LINE__, __FILE__);
    if(initialized==0)
        throw CZigZagException(NOT_INITIALIZED, __LINE__, __FILE__);
    // end sanity check

    int blocksize=block_size;
    int codedsize = I*(J+nb_encoder);

    int nbblock = data.taille/codedsize+(data.taille%codedsize==0?0:1);
    int t,e ;
//    const int iter_max=10;

    DVector output(nbblock*block_size);

    //allocate storage container for parity blocks
    UVector<DVector> parityblock(nb_encoder), Lmen(nb_encoder);
    DVector subdata(block_size), Lpmn(block_size), Lmn(block_size);
    for(e=0;e<nb_encoder;e++) {
        parityblock.vect[e]= DVector(I);
        Lmen.vect[e]= DVector(block_size,0.0);
    }

    // process block by block .. retrieve I*J+E*I data
    for(int n=0;n<nbblock;n++) {
        Lmn = subdata = data.copy(n*codedsize, n*codedsize+blocksize);
        for(e=0;e<nb_encoder;e++) {
            parityblock.vect[e]= data.copy(n*codedsize+blocksize+e*I, n*codedsize+blocksize+(e+1)*I);
            Lmen.vect[e].reset(0.0);
        }
        // have extracted relevant data .. now .. go process them, encoder by encoder .., iteration by iteration ..
        for(int it=0;it<iter_max;it++) {
            for(e=0;e<nb_encoder;e++) {
                Lpmn=Lmn-Lmen.vect[e];
                // getSoftValue(in, in, in, out, out)
                getSoftValues(Lpmn, parityblock, e, Lmn, Lmen);
            }
        }
        // Lmn is the last soft value ... send it home ...
        output.insert(Lmn, n*block_size);
    }
    return output.copy(0, output.taille-superpadsize);
}

//#define debugconcat

void WConcatZigZag::getSoftValues(DVector &Lpmn, UVector<DVector> par, int e, DVector &Lmn, UVector<DVector> Lmen) {
    DVector in = interleaver.vect[e].Apply(Lpmn);

    double min, dv;

    DVector min1(I,M_HUGE_NUM), min2(I,M_HUGE_NUM), fv(I+1), bv(I+1);
    WVector min1_idx(I,0), min2_idx(I,0),
            sign(I,1), fsign(I,1),bsign(I,1);

    /*DMatrix Ld(I,J);*/ DVector Fp(I+1), Bp(I+1);
    /*Ld.reset(); */       Fp.reset();        Bp.reset();
    Fp.vect[0]=M_HUGE_NUM;
    Bp.vect[I]=par.vect[e].vect[I-1];
        int dsign = 1;
    // processing segment by segment
    for(int ip=1,i=0;i<I;i++,ip++) { // get W(.) without F(p) and B(p)
        // calculate W(d(i,1),...,d(i,J)), identify two mins in the segment..
        sign.vect[i] = 1;
        min1.vect[i] = min2.vect[i] = M_HUGE_NUM; // arbitrary initialize a min
        min1_idx.vect[i] = min2_idx.vect[i] = -1; // identify the min position

        for(int j = 0; j < J;j++) { // look for another min and calculate W(.)
            sign.vect[i] *= (dv = in.vect[i*J+j]>=0.0?1:-1); // get the sign of current data
            dv *= in.vect[i*J+j]; // get the absolute value of segment.vect[j]
            if(dv<=min1.vect[i]) { // current value is min
                min2.vect[i]=min1.vect[i]; min1.vect[i]=dv; // update min table
                min2_idx.vect[i] = min1_idx.vect[i]; min1_idx.vect[i] = j;
            } else if (dv <=min2.vect[i]) { // second min value
                min2.vect[i] = dv;            // update min table
                min2_idx.vect[i] = j;
            }
        }
    }
    for(int ip=1,i=0;i<I;i++,ip++) { // get forward MLA
        // calculate F(p(i)) based on W(d(i,1),...,d(i,J)) calculated previously
        fsign.vect[i] = sign.vect[i]*(fv.vect[ip-1]=(Fp.vect[ip-1]>=0?1:-1)); // calculate the sign ..
        fv.vect[ip-1] *= Fp.vect[ip-1]; // fabs(F(p(i-1)))
        Fp.vect[ip] = par.vect[e].vect[i]; // F(p(i)=p(i) {+W(F(p(i-1)),d(i,1),...,d(i,J))}
        if (fv.vect[ip-1]<=min1.vect[i]) { // if (F(p(i-1)) <= min1 => W(.) = sign*fabs(F(p(i-1))
            Fp.vect[ip] += fsign.vect[i]*fv.vect[ip-1];
        } else { // if (F(p(i-1)) > min1 => W(.) = sign*min1
            Fp.vect[ip] += fsign.vect[i]*min1.vect[i];
        }
    }
    for(int ip=I,i=I-1;i>=0;i--,ip--) { // get backward MLA
        // calculate B(p(i-1)) based on W(d(i,1),...,d(i,J)) calculated previously
        bsign.vect[i] = sign.vect[i]*(bv.vect[ip]=(Bp.vect[ip]>=0?1:-1)); // calculate the sign
        bv.vect[ip] *= Bp.vect[ip]; // fabs(B(p(i)))
        Bp.vect[ip-1] = i ? par.vect[e].vect[i-1] : M_HUGE_NUM ; //  B(p(i-1)) = p(i-1) {+W(d(i,1),...,d(i,J),B(p(i)))}
        if(bv.vect[ip]<=min1.vect[i]) { // if B(p(i))<=min1 => W(.) = sign*fabs(B(p(i))
            Bp.vect[ip-1] += bsign.vect[i]*bv.vect[ip];
        } else { // if B(p(i)) > min1 => W(.) = sign*min1
            Bp.vect[ip-1] += bsign.vect[i]*min1.vect[i];
        }
    }
    for(int ip=1,i=0;i<I;i++,ip++) { // calculate L(d(i,j))
        int fpsign = (Fp.vect[ip-1]>=0?1:-1), bpsign = (Bp.vect[ip]>=0?1:-1);
        for (int j=0;j<J;j++) {
            // L(d(i,j)) =d(i,j) {+W(F(p(i-1)),d(i,1),...,d(i,j-1),d(i,j+1),...,d(i,J), B(p(i))}
            Lmn.vect[i*J+j]=in.vect[i*J+j];

            dsign = sign.vect[i]*(in.vect[i*J+j]>=0?1:-1); // remove d(i,j) from the sign
            dsign *= fpsign*bpsign; //calculate new sign

            min = min1.vect[i];

            if (j==min1_idx.vect[i])//check if d(i,j) is the current min ?
                min = min2.vect[i];
            if (min>fv.vect[ip-1])  // check if F(p(i-1)) is the min
                min = fv.vect[ip-1];
            if (min>bv.vect[ip])    // check if B(p(i)) is the min
                 min = bv.vect[ip];
            Lmen.vect[e].vect[i*J+j] = min*dsign;
            Lmn.vect[i*J+j] += Lmen.vect[e].vect[i*J+j];
        }
    }
    interleaver.vect[e].ipExtract(Lmn);
    interleaver.vect[e].ipExtract(Lmen.vect[e]);

}

#else
WVector WConcatZigZag::decode(WVector &data){
    return WVector();
}
#endif

WVector WConcatZigZag::getparity(WVector &data) {
    if (data.vect ==NULL || data.taille<=0) throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if(data.taille!=I*J) throw CZigZagException(INVALID_ZIGZAG_PARAMETER, __LINE__, __FILE__);

    WVector parity(I);

    int par=0,p=0;
    for(int i=0;i<data.taille;i++) {
        if(i && i%J==0) {
            par = parity.vect[p++]= par^data.vect[i];
        } else {
            par^=data.vect[i];
        }
    }
    parity.vect[p++] = par;

    return parity;
}

void WConcatZigZag::getparity(WVector &data, WVector &out, int parity_position) {
    if (data.vect ==NULL || data.taille<=0) throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if(data.taille!=I*J) throw CZigZagException(INVALID_ZIGZAG_PARAMETER, __LINE__, __FILE__);

    int par=0, p=parity_position;
    for(int i=0;i<data.taille;i++) {
        if(i && i%J==0) {
            out.vect[p++]= par;
        }
        par^=data.vect[i];
    }
    out.vect[p++]= par;
}

WVector WConcatZigZag::encode(WVector &data) {
    if (data.vect ==NULL || data.taille<=0) throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    block_size=I*J; padsize=0;
    int nbblock = data.taille/block_size+(data.taille%block_size==0?0:1), paritysize=I*nb_encoder;
    int outsize = (I*J+paritysize)*nbblock, t, pbsize=block_size+paritysize;
    superpadsize = nbblock*block_size-data.taille;
    WVector out(outsize,8), data1(nbblock*block_size,padding_value), datablock;
    data1.insert(data,0);
    WVector tmp(block_size) ;
    for(int i=0;i<nbblock;i++){
        datablock = data1.copy(i*block_size, ( i+1)*block_size);
        out.insert(datablock, i*pbsize);
        for(int e=0;e<nb_encoder;e++) {
            tmp = interleaver.vect[e].Apply(datablock);
            getparity(tmp , out, i*pbsize+block_size+e*I);
        }
    }
    return out;
}


WVector WAConcatZigZag::encode(WVector &data) {
    if (data.vect ==NULL || data.taille<=0) throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    block_size=I*J*K; padsize=0;
    int subcodedsize = I*(J+nb_encoder);
    int codedblocksize = subcodedsize*K;
    int subblock_size = I*J;
    int nbblock = data.taille/block_size+(data.taille%block_size==0?0:1), paritysize=I*nb_encoder;
    int outsize = K*(I*J+paritysize)*nbblock, t, pbsize=subblock_size+paritysize;
    superpadsize = nbblock*block_size-data.taille;
    WVector out(outsize,0), data1(nbblock*block_size,padding_value), datablock;
    data1.insert(data,0);
    WVector subdata(subblock_size), edatablock(block_size) ;
    for(int i=0;i<nbblock;i++){
        datablock = data1.copy(i*block_size, ( i+1)*block_size);
//        WVector testblock(block_size,10+i);
        for(int k=0;k<K;k++) {
            subdata = datablock.copy(k*subblock_size, (k+1)*subblock_size);
            out.insert(subdata, i*codedblocksize+k*subcodedsize); //+k*pbsize);
        }
//        cout << "+++ " << out << "---" << endl;
        for(int e=0;e<nb_encoder;e++) {
            edatablock = interleaver.vect[e].Apply(datablock);
            // break into smaller chunk
            for(int k=0;k<K;k++) {
                subdata = edatablock.copy(k*subblock_size,(k+1)*subblock_size);
                getparity(subdata, out, i*codedblocksize+k*subcodedsize+subblock_size+e*I);
            }
        }
    }
    return out;
}

void WAConcatZigZag::getparity(WVector &data, WVector &out, int parity_position) {
    if (data.vect ==NULL || data.taille<=0) throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if(data.taille!=I*J) throw CZigZagException(INVALID_ZIGZAG_PARAMETER, __LINE__, __FILE__);

    int par=0, p=parity_position;
    for(int i=0;i<data.taille;i++) {
        if(i && i%J==0) {
            out.vect[p++]= par;
        }
        par^=data.vect[i];
    }
    out.vect[p++]= par;
}

void WAConcatZigZag::getSoftValues(DVector &Lpmn, UVector<UDVector> & par, int e, int k, DVector &Lmn, UVector<DVector> & Lmen) {
    int subblock_size = I*J;
    DVector in = Lpmn.copy(k*subblock_size,(k+1)*subblock_size);

    double min, dv;

    DVector min1(I,M_HUGE_NUM), min2(I,M_HUGE_NUM), fv(I+1), bv(I+1);
    WVector min1_idx(I,0), min2_idx(I,0),
            sign(I,1), fsign(I,1),bsign(I,1);

    /*DMatrix Ld(I,J);*/ DVector Fp(I+1), Bp(I+1);
    /*Ld.reset(); */       Fp.reset();        Bp.reset();
    Fp.vect[0]=M_HUGE_NUM;
    Bp.vect[I]=par.vect[e].vect[k].vect[I-1];
        int dsign = 1;
    // processing segment by segment
    for(int ip=1,i=0;i<I;i++,ip++) { // get W(.) without F(p) and B(p)
        // calculate W(d(i,1),...,d(i,J)), identify two mins in the segment..
        sign.vect[i] = 1;
        min1.vect[i] = min2.vect[i] = M_HUGE_NUM; // arbitrary initialize a min
        min1_idx.vect[i] = min2_idx.vect[i] = -1; // identify the min position

        for(int j = 0; j < J;j++) { // look for another min and calculate W(.)
            sign.vect[i] *= (dv = in.vect[i*J+j]>=0.0?1:-1); // get the sign of current data
            dv *= in.vect[i*J+j]; // get the absolute value of segment.vect[j]
            if(dv<=min1.vect[i]) { // current value is min
                min2.vect[i]=min1.vect[i]; min1.vect[i]=dv; // update min table
                min2_idx.vect[i] = min1_idx.vect[i]; min1_idx.vect[i] = j;
            } else if (dv <=min2.vect[i]) { // second min value
                min2.vect[i] = dv;            // update min table
                min2_idx.vect[i] = j;
            }
        }
    }
    for(int ip=1,i=0;i<I;i++,ip++) { // get forward MLA
        // calculate F(p(i)) based on W(d(i,1),...,d(i,J)) calculated previously
        fsign.vect[i] = sign.vect[i]*(fv.vect[ip-1]=(Fp.vect[ip-1]>=0?1:-1)); // calculate the sign ..
        fv.vect[ip-1] *= Fp.vect[ip-1]; // fabs(F(p(i-1)))
        Fp.vect[ip] = par.vect[e].vect[k].vect[i]; // F(p(i)=p(i) {+W(F(p(i-1)),d(i,1),...,d(i,J))}
        if (fv.vect[ip-1]<=min1.vect[i]) { // if (F(p(i-1)) <= min1 => W(.) = sign*fabs(F(p(i-1))
            Fp.vect[ip] += fsign.vect[i]*fv.vect[ip-1];
        } else { // if (F(p(i-1)) > min1 => W(.) = sign*min1
            Fp.vect[ip] += fsign.vect[i]*min1.vect[i];
        }
    }
    for(int ip=I,i=I-1;i>=0;i--,ip--) { // get backward MLA
        // calculate B(p(i-1)) based on W(d(i,1),...,d(i,J)) calculated previously
        bsign.vect[i] = sign.vect[i]*(bv.vect[ip]=(Bp.vect[ip]>=0?1:-1)); // calculate the sign
        bv.vect[ip] *= Bp.vect[ip]; // fabs(B(p(i)))
        Bp.vect[ip-1] = i ? par.vect[e].vect[k].vect[i-1] : M_HUGE_NUM ; //  B(p(i-1)) = p(i-1) {+W(d(i,1),...,d(i,J),B(p(i)))}
        if(bv.vect[ip]<=min1.vect[i]) { // if B(p(i))<=min1 => W(.) = sign*fabs(B(p(i))
            Bp.vect[ip-1] += bsign.vect[i]*bv.vect[ip];
        } else { // if B(p(i)) > min1 => W(.) = sign*min1
            Bp.vect[ip-1] += bsign.vect[i]*min1.vect[i];
        }
    }
    for(int ip=1,i=0;i<I;i++,ip++) { // calculate L(d(i,j))
        int fpsign = (Fp.vect[ip-1]>=0?1:-1), bpsign = (Bp.vect[ip]>=0?1:-1);
        for (int j=0;j<J;j++) {
            // L(d(i,j)) =d(i,j) {+W(F(p(i-1)),d(i,1),...,d(i,j-1),d(i,j+1),...,d(i,J), B(p(i))}
            Lmn.vect[k*I*J+i*J+j]=in.vect[i*J+j];

            dsign = sign.vect[i]*(in.vect[i*J+j]>=0?1:-1); // remove d(i,j) from the sign
            dsign *= fpsign*bpsign; //calculate new sign

            min = min1.vect[i];

            if (j==min1_idx.vect[i])//check if d(i,j) is the current min ?
                min = min2.vect[i];
            if (min>fv.vect[ip-1])  // check if F(p(i-1)) is the min
                min = fv.vect[ip-1];
            if (min>bv.vect[ip])    // check if B(p(i)) is the min
                 min = bv.vect[ip];
            Lmen.vect[e].vect[k*I*J+i*J+j] = min*dsign;
            Lmn.vect[k*I*J+i*J+j] += Lmen.vect[e].vect[k*I*J+i*J+j];
        }
    }
}


DVector WAConcatZigZag::decode(DVector &data, int iter_max){
    if (data.vect ==NULL || data.taille<=0)
        throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if(block_size<1)
        throw CZigZagException(INVALID_ZIGZAG_PARAMETER, __LINE__, __FILE__);
    if(initialized==0)
        throw CZigZagException(NOT_INITIALIZED, __LINE__, __FILE__);
    // end sanity check

    int subblock_size =I*J;
    //int pbsize = I*(J+nb_encoder);
    int codedsize = K*I*(J+nb_encoder);

    int nbblock = data.taille/codedsize+(data.taille%codedsize==0?0:1);
    int  t,e, pbsize=I*(J+nb_encoder) ;
//    const int iter_max=10;

    DVector output(nbblock*block_size);

    //allocate storage container for parity blocks
    UVector<UDVector> parityblock(nb_encoder);
    UVector<DVector> bLmen(nb_encoder);
    DVector bigdata(codedsize), rawdata(block_size), bLmn(block_size);

    DVector subdata(block_size), bLpmn(block_size), eLpmn(block_size);
#ifdef PROFILING_CHECK
       cout << __FILE__ <<"@"<<__LINE__<< ":" << time(NULL) <<endl;
#endif       
    
    for(e=0;e<nb_encoder;e++) {
        parityblock.vect[e]= UDVector(K);
        bLmen.vect[e] = DVector(block_size);
        for(int k=0;k<K;k++) {
#ifdef PROFILING_CHECK
       cout << __FILE__ <<"@"<<__LINE__<< ":" << time(NULL) <<endl;
#endif       
            parityblock.vect[e].vect[k] = DVector(I);
        }
    }
#ifdef PROFILING_CHECK
       cout << __FILE__ <<"@"<<__LINE__<< ":" << time(NULL) <<endl;
#endif       

    // process block by block .. retrieve I*J+E*I data
    for(int n=0;n<nbblock;n++) {
#ifdef PROFILING_CHECK
       cout << __FILE__ <<"@"<<__LINE__<< ":" << time(NULL) <<endl;
#endif       
        bigdata = data.copy(n*codedsize, (n+1)*codedsize); // big block extracted
#ifdef PROFILING_CHECK
       cout << __FILE__ <<"@"<<__LINE__<< ":" << time(NULL) <<endl;
#endif       
        // extract information bit from big block ... initialize working chunks
        for(int k=0;k<K;k++) {
            rawdata.insert(bigdata.copy(k*pbsize, k*pbsize+subblock_size), k*subblock_size);
            for(int e=0;e<nb_encoder;e++) {
                parityblock.vect[e].vect[k] = bigdata.copy(k*pbsize+subblock_size+e*I, k*pbsize+subblock_size+(e+1)*I);
                bLmen.vect[e].reset(0.0);
            }
        }
#ifdef PROFILING_CHECK
       cout << __FILE__ <<"@"<<__LINE__<< ":" << time(NULL) <<endl;
#endif       

        bLmn = rawdata;
        for(int it=0;it<iter_max;it++) {
            for(int e=0;e<nb_encoder;e++) {
                bLpmn = bLmn - bLmen.vect[e];
                eLpmn = interleaver.vect[e].Apply(bLpmn);
                for(int k=0;k<K;k++) { // extract kecil kecilan
                    getSoftValues(eLpmn, parityblock, e,k, bLmn, bLmen);
                }
                interleaver.vect[e].ipExtract(bLmn);
                interleaver.vect[e].ipExtract(bLmen.vect[e]);
            }
        }
#ifdef PROFILING_CHECK
       cout << __FILE__ <<"@"<<__LINE__<< ":" << time(NULL) <<endl;
#endif       

        output.insert(bLmn, n*block_size);
    }
    return output.copy(0, output.taille-superpadsize);
}



#undef  ZIGZAG_DEBUG
WVector zigzag_encode(WVector &data, int seglen, int &initial_parity, int parity, int padvalue){
    if (data.vect ==NULL || data.taille<=0) throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if (seglen<3) throw CZigZagException(INVALID_SEGMENT_LENGTH, __LINE__, __FILE__);
    int data_per_segment = seglen-2;
    int nb_segment = data.taille/data_per_segment+((data.taille%data_per_segment)!=0?1:0);
    int padding_size = data_per_segment*nb_segment-data.taille;
    int outsize = nb_segment*(seglen-1);
#ifdef ZIGZAG_DEBUG
    fprintf(stderr, "data.taille:seglen:data_per_segment:padding_size:outsize:nb_segment\r\n      %d    :  %d  :          %d    :     %d     :   %d    :     %d\r\n", data.taille, seglen, data_per_segment, padding_size, outsize, nb_segment);
#endif
    WVector out(outsize);
    int i=0,j=0, par=initial_parity;
    for(i=0,j=0, par=initial_parity;i<data.taille;i++) {
        if (i>0 && (i%data_per_segment)==0) {
#ifdef ZIGZAG_DEBUG
            cout << " eiiii " << par << endl;
#endif
            par= (out.vect[j++] = ((par==parity) ? 0: 1));
        }
        par ^= (out.vect[j++] = data.vect[i]);

#ifdef ZIGZAG_DEBUG
        cout << out;
#endif
    }

    for(i=0;i<padding_size;i++) {
        par ^= (out.vect[j++] = padvalue);
    }
    initial_parity = par^ (out.vect[j++] = (par==parity) ? 0: 1);
    return out;
}

void show_zigzag_code(WVector &data, int segment_length) {
    if (data.vect ==NULL || data.taille<=0) throw CZigZagException(EMPTY_VECTOR,  __LINE__, __FILE__);
    if (segment_length<3) throw CZigZagException(INVALID_SEGMENT_LENGTH, __LINE__, __FILE__);
    cout << endl << endl;
    for(int i = 0; i<data.taille;i++) {
        if(i>0 && (i%(segment_length-1))==0)cout << endl;
        cout << data.vect[i] << " ";

    }
}


