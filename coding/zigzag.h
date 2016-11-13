/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/21 10:56:58 $
 * $Revision: 1.1.2.10 $
 * $Id: zigzag.h,v 1.1.2.10 2006/07/21 10:56:58 syed Exp $
 ********************************************************************/

#ifndef _ZIGZAG_H_
#define _ZIGZAG_H_
#include "tools/all.h"
#include <time.h>

#define INVALID_SEGMENT_LENGTH                  40
#define INVALID_ZIGZAG_PARAMETER                41
#define INVALID_NUMBER_OF_ENCODER               42
#define NOT_INITIALIZED                         43
#define INVALID_PACKET_SIZE                     44
#define INVALID_PARITY_SIZE                     45
#define CREATE_OPTIMIZED_PCZZ_FAILED            46

#define OUTPUT_TYPE_IS_BIT_LLR      10021
#define OUTPUT_TYPE_IS_BIT_APP      10022

static int zz_id=0;

typedef UVector<DVector> UDVector;

class CZigZagException {
    protected:
        void print() {
            cerr << "CZigZagException : " << err_mesg;
            if(line!=-1) cerr << " at "<< filename << "@" << line;
            cerr << endl;
            switch (err_mesg) {
                case EMPTY_VECTOR :
                    cerr << "CZigZagException : The data vector is empty." << endl;
                    break;
                case INVALID_SEGMENT_LENGTH	:
                    cerr << "CZigZagException : The zigzag code segment size is invalid" << endl;
                    break;
                case INVALID_ZIGZAG_PARAMETER :
                    cerr << "CZigZagException : The zigzag code parameter is invalid" << endl;
                    break;
                case INVALID_NUMBER_OF_ENCODER :
                    cerr << "CZigZagException : Invalid number of zigzag encoder  " << endl;
                    break;
                case NOT_INITIALIZED :
                    cerr << "CZigZagException : The zigzag code is uninitialized" << endl;
                    break;
                case INVALID_PACKET_SIZE :
                    cerr << "CZigZagException : The data packet size is incompatible with predefineed I and J" << endl;
                    break;
                case INVALID_PARITY_SIZE :
                    cerr << "CZigZagException : The parity packet size is incompatible with predefineed I" << endl;
                    break;
                default:
                    cerr << "CZigZagException : Unknown error. Refer to error definition" << endl;
                    break;
            }
        }

    public :
        int err_mesg, line;
        char filename[255];
		CZigZagException(int erm, int line, char *filename) { err_mesg = erm; this->line=line, strcpy(this->filename,filename); print();}
		CZigZagException(int erm) { err_mesg = erm; line=-1; strcpy(filename,"Unspecified"); print(); }
};

class WZigZag{
    protected:
        int I,J,H, padsize, block_size, superpadsize, initialized, myid;

    // I = number of segment , J = data per segment, H = segment_length-1
        WVector parity;
    public:
        int segment_length, initial_parity, parity_rule, padding_value, lastparity;
        int concat_mode, decoder_output_type;
        WZigZag(){
            initialized=0;
        }
        WZigZag(int seglen, int initpar=0, int parrule=0, int padvalue=0){
            if(seglen<3) throw CZigZagException(INVALID_SEGMENT_LENGTH);
            if((initpar&0xfffffffe) || (parrule&0xfffffffe) || (padvalue&0xfffffffe)) throw CZigZagException(INVALID_ZIGZAG_PARAMETER);
            segment_length = seglen;
            lastparity=initial_parity=initpar;
            parity_rule = parrule;
            padding_value=padvalue;
            J = segment_length-2;
            I=padsize=0;
            H=J+1;
            initialized=1;
            concat_mode = 0;
            decoder_output_type = OUTPUT_TYPE_IS_BIT_APP;
            //cout << "New zigzag decoder installed" <<endl;
            myid=zz_id++;
        }
        WVector encode(WVector &data, int blocksize);
        DVector decode(DVector &data);
        DVector decode(DVector &data, DVector &llr);
};

class WConcatZigZag: public WZigZag {
    private:
        int initialized;
    protected:
        UVector<CInterleaver> interleaver;
        WVector getparity(WVector &data);
        void getparity(WVector &data, WVector &out, int parity_position);
        void getSoftValues(DVector &Lpmn, UVector<DVector> par, int e, DVector &Lmn, UVector<DVector> Lmen);

    public:
        int nb_encoder;

        WConcatZigZag(int encoder_count, int I, int J, int initpar=0, int parrule=0, int padvalue=0):
                    WZigZag(J+2,initpar,parrule,padvalue){

            if(encoder_count<1) throw CZigZagException(INVALID_NUMBER_OF_ENCODER, __LINE__, __FILE__);

            int xtime = time(NULL);

            interleaver = UVector<CInterleaver>(nb_encoder=encoder_count);
            block_size=I*J;
            this->I=I;

            for(int i = 0;i<nb_encoder;i++) {
                interleaver.vect[i]=CInterleaver(block_size, xtime+block_size*i*2)  ;
            }
            initialized=1;
            myid=zz_id++;
        }

        WConcatZigZag(char type, int encoder_count, int I, int J, int initpar=0, int parrule=0, int padvalue=0):
                    WZigZag(J+2,initpar,parrule,padvalue){
            if(encoder_count<1) throw CZigZagException(INVALID_NUMBER_OF_ENCODER, __LINE__, __FILE__);

//             int xtime = time(NULL);

            interleaver = UVector<CInterleaver>(nb_encoder=encoder_count);
            block_size=I*J;
            this->I=I;

            WVector basic(block_size, 0, 'i'), tmp2(block_size,0,'i');
            CInterleaver tmp;
            tmp = interleaver.vect[0]=CInterleaver('o', I,J);

            int i=1,n=0;
            for(i = 1,n=0;i<nb_encoder && n<2*nb_encoder;n++) {
                tmp.ipApply(basic);
//                cout << basic;
//                cout << tmp.GetMap()-tmp2;
                tmp2=tmp.Apply(basic);
  //              cout << tmp.GetMap()-tmp2;
                tmp=CInterleaver(tmp2);
    //            cout << tmp.GetMap()-tmp2;
                bool valid = true;
                for(int t=0;t<i;t++) {
      //              cout << interleaver.vect[t].GetMap()-tmp2;
                    if(tmp==interleaver.vect[t]) { valid=false; break; }
                }
                if(valid) {
                    interleaver.vect[i++] = tmp;
                    n=-1;
                }
            }
            if(i!=nb_encoder) {
  //               cerr << "Oouppss" << endl << basic;
   //              cerr << " Alerte a malibu : i = " << i << " , n = " << n << endl;
                throw CZigZagException(CREATE_OPTIMIZED_PCZZ_FAILED, __LINE__, __FILE__);

            }
            initialized=1;
            myid=zz_id++;
        }

        void  init(int blocksize) {
            block_size = blocksize;
            initialized = 1;
        }

        WVector encode(WVector &data);
        DVector decode(DVector &data, int iter_max=20);
        void printInterleaver(){ cout << interleaver; }

};

/*
class WOConcatZigZag : public WConcatZigZag {
    protected:
        void OptimizedInterleaver();
};
*/



class WAConcatZigZag:public WZigZag {
    private:
        int initialized;

    protected:
        UVector<CInterleaver> interleaver;
        WVector getparity(WVector &data);
        void getparity(WVector &data, WVector &out, int parity_position);
        void getSoftValues(DVector &Lpmn, UVector<UDVector>& par, int e, int k, DVector &Lmn, UVector<DVector>& Lmen) ;
        int K; // interleaver multiplier

    public:
        int nb_encoder;

        WAConcatZigZag(int encoder_count, int I, int J, int K, int initpar=0, int parrule=0, int padvalue=0):
                    WZigZag(J+2,initpar,parrule,padvalue){

            if(encoder_count<1) throw CZigZagException(INVALID_NUMBER_OF_ENCODER, __LINE__, __FILE__);

            int xtime = time(NULL);

            interleaver = UVector<CInterleaver>(nb_encoder=encoder_count);
            block_size=I*J;
            this->I=I;
            this->K=K;

            for(int i = 0;i<nb_encoder;i++) {
                interleaver.vect[i]=CInterleaver(block_size*K, xtime+block_size*i*2)  ;
            }
            initialized=1;
            myid=zz_id++;
        }

        WAConcatZigZag(char type, int encoder_count, int I, int J, int K, int initpar=0, int parrule=0, int padvalue=0):
                    WZigZag(J+2,initpar,parrule,padvalue){
            if(encoder_count<1) throw CZigZagException(INVALID_NUMBER_OF_ENCODER, __LINE__, __FILE__);

//            int xtime = time(NULL);

            interleaver = UVector<CInterleaver>(nb_encoder=encoder_count);
            block_size=I*J;
            this->I=I;
            this->K=K;

            WVector basic(block_size*K, 0, 'i');
            CInterleaver tmp;
            tmp = interleaver.vect[0]=CInterleaver('o', I*K,J);

            int i=1,n=0;
            for(i = 1,n=0;i<nb_encoder && n<nb_encoder;n++) {
                tmp.ipApply(basic);
                tmp=CInterleaver(tmp.Apply(basic));
                bool valid = true;
                for(int t=0;t<i;t++) {
                    //if(tmp.GetMap().diff_count(interleaver.vect[t].GetMap())!=0) { valid=false; break; }
                    if(tmp==interleaver.vect[t]) { valid=false; break; }
                }
                if(valid) {
                    interleaver.vect[i++] = tmp;
                    n=-1;
                }
            }
            if(i!=nb_encoder) {
                throw CZigZagException(CREATE_OPTIMIZED_PCZZ_FAILED, __LINE__, __FILE__);

            }
            initialized=1;
            myid=zz_id++;
        }

        void  init(int blocksize) {
            block_size = blocksize;
            initialized = 1;
        }

        WVector encode(WVector &data);
        DVector decode(DVector &data, int iter_max=20);



};
WVector zigzag_encode(WVector &data, int seglen, int &initial_parity, int parity=0, int padvalue=0);
void show_zigzag_code(WVector &data, int segment_length);

#endif
