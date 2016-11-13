#include <iostream.h>
#include <time.h>
#include <tools/all.h>
#include <coding/zigzag.h>


//#define CZIGZAG_DEBUG
//#define CZIGZAG_VERBOSE

//#define VERBOSE

#ifdef VERBOSE2
    #ifndef VERBOSE
        #define VERBOSE
    #endif
#endif

//#define TEST_DATA

#define ERROR_BASED_NB_SIM 10000
#define SIM_UNTIL_MIN_ERROR_ACHIEVED
#define MIN_ERROR_REQUIRED 100
#define MIN_ERROR 10

int main1(int argc, char ** argv) {
    int J = 4, I  = 127, initpar=0, multiplier=30;
    int blocksize = I*J, data_stream_size = multiplier*blocksize;
    int nb_data_stream = 100, nberr=0,nbbit=0, err=0;
    double noise_variance=0.0025, snr_min=0,snr_max=12, snr_step=.5, snr;
    WVector data(data_stream_size), coded, decoded;
    DVector noise, noisy, ld;
   // WZigZag zzcoder(J+2,initpar);

    TRandBitStatePtr bitptr;
    bitptr.iseed = time(NULL);
    dRandGausStatePtr nptr;
    char filename1[100], filename2[100];
    sprintf(filename1, "bsk_perf-%lu.log", bitptr.iseed);
    sprintf(filename2, "bsk_perf_final-%lu.log", bitptr.iseed);
    ofstream os(filename1), os2(filename2);

    cout << " HELLO DARLING ... arrrrrrrrrrrrrrrrr .... " << endl;
    os << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os << "# J                = " << J << endl;
    os << "# I                = " << I << endl;
    os << "# initpar          = " << initpar << endl;
    os << "# multiplier       = " << multiplier << endl;
    os << "# blocksize        = " << blocksize << endl;
    os << "# data stream size = " << data_stream_size << endl;
    os << "# mn data stream   = " << nb_data_stream << endl;
    os << "# snr min          = " << snr_min << endl;
    os << "# snr max          = " << snr_max << endl;
    os << "# snr step         = " << snr_step << endl;
    os << "# **** End parameter *** " << endl;
    os << "# **** DATA **** " << endl;

    os << "snr:noise variance:coded.taille:decoded.taille:err:nberr:nbbit"<<endl;

    os2 << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os2 << "# J                = " << J << endl;
    os2 << "# I                = " << I << endl;
    os2 << "# initpar          = " << initpar << endl;
    os2 << "# multiplier       = " << multiplier << endl;
    os2 << "# blocksize        = " << blocksize << endl;
    os2 << "# data stream size = " << data_stream_size << endl;
    os2 << "# mn data stream   = " << nb_data_stream << endl;
    os2 << "# snr min          = " << snr_min << endl;
    os2 << "# snr max          = " << snr_max << endl;
    os2 << "# snr step         = " << snr_step << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;

    os2 << "J:I:initpar:blocksize:snr:noise variance:nberr:nbbit:BER"<<endl;


    for(snr = snr_min; snr<=snr_max;snr+=snr_step){
        double rate = J/(J+1.0);
        cout << endl << " Simulating for snr = " << snr << endl;
        noise_variance = convertsnrtosigma2(rate, snr);
        dRandGausInit(&nptr,  0, noise_variance);
        nberr = nbbit = 0;
        for(int t=0;t<nb_data_stream;t++){
            bRandBit1NRZi(&bitptr, data);
#ifdef VERBOSE
            cout << "+";
#endif
#ifdef VERBOSE
            cout <<"|";
#endif
            noise=DVector(data.taille);
            dbRandGaus(&nptr, noise);
            noisy = noise+data;
#ifdef VERBOSE
            cout <<"/";
#endif
#ifdef VERBOSE
            cout <<"-";
#endif
            decoded=WVector(noisy.taille);
            for (int p=0;p<decoded.taille;p++) decoded.vect[p] = noisy.vect[p]<0?-1:1;
#ifdef VERBOSE
            cout <<"*";
#endif
            nbbit+=decoded.taille;
            nberr+=err = data.diff_count(decoded);
            os << snr << ":" << noise_variance<< ":" <<coded.taille<< ":" << decoded.taille<< ":" << err<< ":" << nberr<< ":" << nbbit<<endl;
        }
        os2 << J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
    }
    os.close();
    os2.close();
    cout << endl << " * * * * DONE * * * *" << endl;
    cin >> blocksize;
    return -1;
}

int main4_2(int argc, char ** argv) {
    int J, Jmin=2, Jmax=10, Jstep=2,  I  = 127, initpar=0, multiplier=30;
    double noise_variance=0.0025, snr_min=0,snr_max=12, snr_step=1, snr;
    DVector noise, noisy, ld;
    int nb_data_stream = 5000, nberr=0,nbbit=0, err=0;
    WVector data, coded, decoded;

    TRandBitStatePtr bitptr;
    bitptr.iseed = time(NULL);
    dRandGausStatePtr nptr;
    char filename1[100], filename2[100];
    sprintf(filename1, "zigzag_perf-%lu.log", bitptr.iseed);
    sprintf(filename2, "zigzag_perf_final-%lu.log", bitptr.iseed);
    ofstream os(filename1), os2(filename2);

    cout << " HELLO DARLING ... arrrrrrrrrrrrrrrrr .... " << endl;
    os << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os << "# I                              = " << I << endl;
    os << "# initpar                        = " << initpar << endl;
    os << "# multiplier                     = " << multiplier << endl;
    os << "# blocksize(min,max,step)        = (" << I*Jmin << "," << I*Jmax << "," << I*Jstep << ")" <<endl;
    os << "# data stream size(min,max,step) = (" << multiplier*I*Jmin << "," << multiplier*I*Jmax << "," << multiplier*I*Jstep << ")" <<endl;
    os << "# mn data stream                 = " << nb_data_stream << endl;
    os << "# snr min                        = " << snr_min << endl;
    os << "# snr max                        = " << snr_max << endl;
    os << "# snr step                       = " << snr_step << endl;
    os << "# **** End parameter *** " << endl;
    os << "# **** DATA **** " << endl;

    os << "snr:noise variance:coded.taille:decoded.taille:err:nberr:nbbit"<<endl;

    os2 << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I                              = " << I << endl;
    os2 << "# initpar                        = " << initpar << endl;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max,step)        = (" << I*Jmin << "," << I*Jmax << "," << I*Jstep << ")" <<endl;
    os2 << "# data stream size(min,max,step) = (" << multiplier*I*Jmin << "," << multiplier*I*Jmax << "," << multiplier*I*Jstep << ")" <<endl;
    os2 << "# mn data stream                 = " << nb_data_stream << endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;

    os2 << "J:I:initpar:blocksize:snr:noise variance:nberr:nbbit:BER"<<endl;

    for(J=Jmin;J<=Jmax;J+=Jstep) {
        WZigZag zzcoder(J+2,initpar);
        int blocksize = I*J, data_stream_size = multiplier*blocksize;
        data = WVector(data_stream_size);
        for(snr = snr_min; snr<=snr_max;snr+=snr_step){
            double rate = J/(J+1.0);
            noise_variance = convertsnrtosigma2(rate, snr);
            cout << endl << " Simulating for snr = " << snr << "@ J = " << J << endl;
            dRandGausInit(&nptr,  0, noise_variance);
            nberr = nbbit = 0;
            for(int t=0;t<nb_data_stream;t++){
                bRandBit1(&bitptr, data);
#ifdef VERBOSE
                cout << "+";
#endif
                coded = zzcoder.encode(data, blocksize);
#ifdef VERBOSE
                cout <<"|";
#endif
                noise=DVector(coded.taille);
                dbRandGaus(&nptr, noise);
                noisy = noise+coded;
#ifdef VERBOSE
                cout <<"/";
#endif
                ld = zzcoder.decode(noisy);
#ifdef VERBOSE
                cout <<"-";
#endif
                decoded=WVector(ld.taille);
                for (int p=0;p<decoded.taille;p++) decoded.vect[p] = ld.vect[p]<0?1:0;
#ifdef VERBOSE
                cout <<"*";
#endif
                nbbit+=decoded.taille;
                nberr+=err = data.diff_count(decoded);
                os << snr << ":" << noise_variance<< ":" <<coded.taille<< ":" << decoded.taille<< ":" << err<< ":" << nberr<< ":" << nbbit<<endl;
            }
            os2 << J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
            if(nberr==0) break;
        }
    }
    os.close();
    os.close();
    cout << endl << " * * * * DONE * * * *" << endl;
    cin >> initpar;
    return -1;
}
//**** I, J sweep start
int main4_3(int argc, char ** argv) {
    int J, Jmin=2, Jmax=10, Jstep=2,
        I=1, Imin=2, Imax=1024, Istep=2,
        initpar=0, multiplier=30;
    double noise_variance=0.0025, snr_min=0,snr_max=12, snr_step=1, snr;
    DVector noise, noisy, ld;
    int nb_data_stream = 10000, nberr=0,nbbit=0, err=0;
    unsigned long starttime, stoptime;
    WVector data, coded, decoded;

    TRandBitStatePtr bitptr;
    starttime=bitptr.iseed = time(NULL);
    dRandGausStatePtr nptr;
    char filename1[100], filename2[100];
    sprintf(filename1, "zigzag_perf-%lu.log", bitptr.iseed);
    sprintf(filename2, "zigzag_perf_final-%lu.log", bitptr.iseed);
    ofstream os(filename1), os2(filename2);

    cout << " HELLO DARLING ... arrrrrrrrrrrrrrrrr .... " << endl;
    os << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os << "# I                              = " << I << endl;
    os << "# initpar                        = " << initpar << endl;
    os << "# multiplier                     = " << multiplier << endl;
    os << "# blocksize(min,max,step)        = (" << I*Jmin << "," << I*Jmax << "," << I*Jstep << ")" <<endl;
    os << "# data stream size(min,max,step) = (" << multiplier*I*Jmin << "," << multiplier*I*Jmax << "," << multiplier*I*Jstep << ")" <<endl;
    os << "# mn data stream                 = " << nb_data_stream << endl;
    os << "# snr min                        = " << snr_min << endl;
    os << "# snr max                        = " << snr_max << endl;
    os << "# snr step                       = " << snr_step << endl;
    os << "# **** End parameter *** " << endl;
    os << "# **** DATA **** " << endl;

    os << "snr:noise variance:coded.taille:decoded.taille:err:nberr:nbbit"<<endl;

    os2 << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# initpar                        = " << initpar << endl;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# mn data stream                 = " << nb_data_stream << endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;

    os2 << "J:I:initpar:blocksize:snr:noise variance:nberr:nbbit:BER"<<endl;

    for(I=Imin; I<=Imax;I*=Istep){
        for(J=Jmin;J<=Jmax;J+=Jstep) {
            WZigZag zzcoder(J+2,initpar);
            int blocksize = I*J, data_stream_size = multiplier*blocksize;
            data = WVector(data_stream_size);
            for(snr = snr_min; snr<=snr_max;snr+=snr_step){
                double rate = J/(J+1.0);
                noise_variance = convertsnrtosigma2(rate, snr);
                cout << endl << " Simulating for snr = " << snr << "@ J = " << J << ", I=" << I <<endl;
                dRandGausInit(&nptr,  0, noise_variance);
                nberr = nbbit = 0;
                for(int t=0;t<nb_data_stream;t++){
                    bRandBit1(&bitptr, data);
#ifdef VERBOSE
                    cout << "+";
#endif
                    coded = zzcoder.encode(data, blocksize);
#ifdef VERBOSE
                    cout <<"|";
#endif
                    noise=DVector(coded.taille);
                    dbRandGaus(&nptr, noise);
                    noisy = noise+coded;
#ifdef VERBOSE
                    cout <<"/";
#endif
                    ld = zzcoder.decode(noisy);
#ifdef VERBOSE
                    cout <<"-";
#endif
                    decoded=WVector(ld.taille);
                    for (int p=0;p<decoded.taille;p++) decoded.vect[p] = ld.vect[p]<0?1:0;
#ifdef VERBOSE
                    cout <<"*";
#endif
                    nbbit+=decoded.taille;
                    nberr+=err = data.diff_count(decoded);
                    os << snr << ":" << noise_variance<< ":" <<coded.taille<< ":" << decoded.taille<< ":" << err<< ":" << nberr<< ":" << nbbit<<endl;
                }
                os2 << J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
                if(nberr==0) break;
            }
        }
    }
    stoptime = time(NULL);
    os  << "#stop time              = " << stoptime<<endl;
    os2 << "#stop time              = " << stoptime<<endl;
    os.close();
    os.close();
    cout << endl << " * * * * DONE * * * *" << endl;

    cout << " Stopped at : " << stoptime << endl;
    cout << " Elapsed time (seconds) : " << stoptime - starttime << endl;
    cin >> initpar;
    return -1;
}

//**** I,J sweep end

//  direct concat
int main4_4(int argc, char ** argv) {
    int J, Jmin=4, Jmax=4, Jstep=4,
        I=1, Imin=4, Imax=64, Istep=4,
        initpar=0, multiplier=30;
    int nb_encoder = 4;
    UVector<CInterleaver> interleaver(nb_encoder);
    UVector<WZigZag> encoder(nb_encoder);

    double noise_variance=0.0025, snr_min=0,snr_max=12, snr_step=2, snr;
    DVector noise, noisy, ld;
    int nb_data_stream = 5000, nberr=0,nbbit=0, err=0;
    unsigned long starttime, stoptime;
    WVector data, coded, decoded, test;

    TRandBitStatePtr bitptr;
    starttime=bitptr.iseed = time(NULL);
    dRandGausStatePtr nptr;
    char filename1[100], filename2[100];
    sprintf(filename1, "zigzag_perf-%lu.log", bitptr.iseed);
    sprintf(filename2, "zigzag_perf_final-%lu.log", bitptr.iseed);
    ofstream os(filename1), os2(filename2);

    cout << " HELLO DARLING ... arrrrrrrrrrrrrrrrr .... " << endl;
    os << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os << "# I                              = " << I << endl;
    os << "# initpar                        = " << initpar << endl;
    os << "# multiplier                     = " << multiplier << endl;
    os << "# blocksize(min,max,step)        = (" << I*Jmin << "," << I*Jmax << "," << I*Jstep << ")" <<endl;
    os << "# data stream size(min,max,step) = (" << multiplier*I*Jmin << "," << multiplier*I*Jmax << "," << multiplier*I*Jstep << ")" <<endl;
    os << "# mn data stream                 = " << nb_data_stream << endl;
    os << "# snr min                        = " << snr_min << endl;
    os << "# snr max                        = " << snr_max << endl;
    os << "# snr step                       = " << snr_step << endl;
    os << "# nb encoder                     = " << nb_encoder << endl;
    os << "# **** End parameter *** " << endl;
    os << "# **** DATA **** " << endl;

    os << "snr:noise variance:coded.taille:decoded.taille:err:nberr:nbbit"<<endl;

    os2 << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# initpar                        = " << initpar << endl;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# mn data stream                 = " << nb_data_stream << endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
    os2 << "# nb encoder                     = " << nb_encoder << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;

    os2 << "J:I:initpar:blocksize:snr:noise variance:nberr:nbbit:BER"<<endl;

    for(I=Imin; I<=Imax;I*=Istep){
        for(J=Jmin;J<=Jmax;J+=Jstep) {
            for(int x=0;x<nb_encoder;x++) {
                encoder.vect[x] = WZigZag(J+2,initpar);
                encoder.vect[x].concat_mode = 1;
            }
            encoder.vect[nb_encoder-1].concat_mode=0;
            int blocksize = I*J, data_stream_size = multiplier*blocksize;
#ifndef TEST_DATA
            data = WVector(data_stream_size);
#else
            data = WVector(data_stream_size,1) ;
#endif
            for(snr = snr_min; snr<=snr_max;snr+=snr_step){
                double rate = J/(J+(double)nb_encoder);
                noise_variance = convertsnrtosigma2(rate, snr);
                cout << endl << " Simulating for snr = " << snr << "@ J = " << J << ", I=" << I <<endl;
                dRandGausInit(&nptr,  0, noise_variance);
                nberr = nbbit = 0;
                for(int t=0;t<nb_data_stream;t++){
#ifndef TEST_DATA
                    bRandBit1(&bitptr, data);
#else
                    cerr << " Using test data stream ..."<< endl;
#endif
#ifdef VERBOSE2
                    cout << "+";
#endif
#ifndef CZIGZAG_DEBUG
                    coded = data;
                    for(int x = 0;x<nb_encoder;x++) {
                        interleaver.vect[x] = CInterleaver(coded.taille, starttime+x);
                        test = interleaver.vect[x].Apply(coded);
                        // problem with g++ v3.4.4 .. does not compile if line above is embedded in the next line
                        coded = encoder.vect[x].encode( test, blocksize);
//                        cout << test << coded;
                    }
#else
                    coded = WVector(data.taille);
                    for(int y=0;y<coded.taille;y++) coded.vect[y]=y;
                    for(int x = 0;x<nb_encoder;x++) {
                        interleaver.vect[x] = CInterleaver(coded.taille, starttime+x);
                        coded = interleaver.vect[x].Apply(coded);
                    }
#endif
#ifdef VERBOSE2
                    cout <<"|";
#endif
#ifndef CZIGZAG_DEBUG
                    noise=DVector(coded.taille,0.0);
#ifndef TEST_DATA
                    dbRandGaus(&nptr, noise);
                    noisy = noise+coded;
#else
                    noisy = coded;
#endif
#ifdef CZIGZAG_VERBOSE
                    cout << noise << noisy;
#endif
#else
                    noise=DVector(coded.taille);
                    noisy = coded;
#endif

#ifdef VERBOSE2
                    cout <<"/";
#endif
                    ld = noisy;
#ifndef CZIGZAG_DEBUG
                    for(int x=nb_encoder-1;x>=0;x--){
#ifdef CZIGZAG_VERBOSE
                        cout << ld.copy(0,blocksize);
#endif
                        ld = interleaver.vect[x].Extract(encoder.vect[x].decode(ld));
                    }
#ifdef CZIGZAG_VERBOSE
                    cout << ld.copy(0,blocksize);
                    cin >> blocksize;
                    return 0;
#endif

#else
                    for(int x=nb_encoder-1;x>=0;x--){
                        ld = interleaver.vect[x].Extract(ld);
                    }
                    cout << ld;
                    cin >> blocksize;
                    return 0;
#endif

#ifdef VERBOSE2
                    cout <<"-";
#endif
                    decoded=WVector(ld.taille);
                    for (int p=0;p<decoded.taille;p++) decoded.vect[p] = ld.vect[p]<0?1:0;
#ifdef VERBOSE2
                    cout <<"*";
#endif
                    nbbit+=decoded.taille;
                    nberr+=err = data.diff_count(decoded);
                    os << snr << ":" << noise_variance<< ":" <<coded.taille<< ":" << decoded.taille<< ":" << err<< ":" << nberr<< ":" << nbbit<<endl;
                }
                os2 << J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
#ifdef VERBOSE
                cout << J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
#endif
                if(nberr==0) break;
            }
        }
    }
    stoptime = time(NULL);
    os  << "#stop time              = " << stoptime<<endl;
    os2 << "#stop time              = " << stoptime<<endl;
    os.close();
    os.close();
    cout << endl << " * * * * DONE * * * *" << endl;

    cout << " Stopped at : " << stoptime << endl;
    cout << " Elapsed time (seconds) : " << stoptime - starttime << endl;
    cin >> initpar;
    return -1;
}

// end direct concat


//  direct concat
int main4_5(int argc, char ** argv) {
    int J, Jmin=64, Jmax=64, Jstep=4,
        I, Imin=64, Imax=64, Istep=4,
        initpar=0, multiplier=30;
    int nb_encoder = 10;
    UVector<CInterleaver> interleaver(nb_encoder);
    UVector<WZigZag> encoder(nb_encoder);

    double noise_variance=0.0025, snr_min=0,snr_max=10, snr_step=0.5, snr;
    DVector noise, noisy, ld;
    int nb_data_stream = 50*nb_encoder, nberr=0,nbbit=0, err=0;
    unsigned long starttime, stoptime;
    WVector data, coded, decoded, test;

    TRandBitStatePtr bitptr;
    starttime=bitptr.iseed = time(NULL);
    dRandGausStatePtr nptr;
    char filename1[100], filename2[100];
    sprintf(filename1, "zigzag_perf-%lu.log", bitptr.iseed);
    sprintf(filename2, "zigzag_perf_final-%lu.log", bitptr.iseed);
#ifdef ZIGZAG_FILE_VERBOSE    
    ofstream os(filename1);
#endif
    ofstream os2(filename2);

    cout << " HELLO DARLING ... arrrrrrrrrrrrrrrrr .... " << endl;
#ifdef ZIGZAG_FILE_VERBOSE    
    os << "# **** main4_5 SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os << "# I                              = " << I << endl;
    os << "# initpar                        = " << initpar << endl;
    os << "# multiplier                     = " << multiplier << endl;
    os << "# blocksize(min,max,step)        = (" << I*Jmin << "," << I*Jmax << "," << I*Jstep << ")" <<endl;
    os << "# data stream size(min,max,step) = (" << multiplier*I*Jmin << "," << multiplier*I*Jmax << "," << multiplier*I*Jstep << ")" <<endl;
    os << "# mn data stream                 = " << nb_data_stream << endl;
    os << "# snr min                        = " << snr_min << endl;
    os << "# snr max                        = " << snr_max << endl;
    os << "# snr step                       = " << snr_step << endl;
    os << "# nb encoder                     = " << nb_encoder << endl;
    os << "# **** End parameter *** " << endl;
    os << "# **** DATA **** " << endl;
    os << "snr:noise variance:coded.taille:decoded.taille:err:nberr:nbbit"<<endl;
#endif

    os2 << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# initpar                        = " << initpar << endl;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# mn data stream                 = " << nb_data_stream << endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
    os2 << "# nb encoder                     = " << nb_encoder << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;

    os2 << "nb_encoder:J:I:initpar:blocksize:snr:noise variance:nberr:nbbit:BER"<<endl;

    for(I=Imin; I<=Imax;I*=Istep){
        for(J=Jmin;J<=Jmax;J*=Jstep) {
            for(int x=0;x<nb_encoder;x++) {
                encoder.vect[x] = WZigZag(J+2,initpar);
                encoder.vect[x].concat_mode = 1;
            }
            encoder.vect[nb_encoder-1].concat_mode=0;
            int blocksize = I*J, data_stream_size = multiplier*blocksize;
#ifndef TEST_DATA
            data = WVector(data_stream_size);
#else
            data = WVector(data_stream_size,1) ;
#endif
            for(snr = snr_min; snr<=snr_max;snr+=snr_step){
                double rate = J/(J+(double)nb_encoder);
                noise_variance = convertsnrtosigma2(rate, snr);
                cout << endl << " Simulating for snr = " << snr << "@ J = " << J << ", I=" << I <<endl;
                dRandGausInit(&nptr,  0, noise_variance);
                nberr = nbbit = 0;
                for(int t=0;t<nb_data_stream;t++){
#ifndef TEST_DATA
                    bRandBit1(&bitptr, data);
#else
                    cerr << " Using test data stream ..."<< endl;
#endif
#ifdef VERBOSE2
                    cout << "+";
#endif
#ifndef CZIGZAG_DEBUG
                    coded = data;
                    for(int x = 0;x<nb_encoder;x++) {
                        interleaver.vect[x] = CInterleaver(coded.taille, starttime+x);
                        test = interleaver.vect[x].Apply(coded);
                        // problem with g++ v3.4.4 .. does not compile if line above is embedded in the next line
                        coded = encoder.vect[x].encode( test, blocksize);
//                        cout << test << coded;
                    }
#else
                    coded = WVector(data.taille);
                    for(int y=0;y<coded.taille;y++) coded.vect[y]=y;
                    for(int x = 0;x<nb_encoder;x++) {
                        interleaver.vect[x] = CInterleaver(coded.taille, starttime+x);
                        coded = interleaver.vect[x].Apply(coded);
                    }
#endif
#ifdef VERBOSE2
                    cout <<"|";
#endif
#ifndef CZIGZAG_DEBUG
                    noise=DVector(coded.taille,0.0);
#ifndef TEST_DATA
                    dbRandGaus(&nptr, noise);
                    noisy = noise+coded;
#else
                    noisy = coded;
#endif
#ifdef CZIGZAG_VERBOSE
                    cout << noise << noisy;
#endif
#else
                    noise=DVector(coded.taille);
                    noisy = coded;
#endif

#ifdef VERBOSE2
                    cout <<"/";
#endif
                    ld = noisy;
#ifndef CZIGZAG_DEBUG
                    for(int x=nb_encoder-1;x>=0;x--){
#ifdef CZIGZAG_VERBOSE
                        cout << ld.copy(0,blocksize);
#endif
                        ld = interleaver.vect[x].Extract(encoder.vect[x].decode(ld));
                    }
#ifdef CZIGZAG_VERBOSE
                    cout << ld.copy(0,blocksize);
                    cin >> blocksize;
                    return 0;
#endif

#else
                    for(int x=nb_encoder-1;x>=0;x--){
                        ld = interleaver.vect[x].Extract(ld);
                    }
                    cout << ld;
                    cin >> blocksize;
                    return 0;
#endif

#ifdef VERBOSE2
                    cout <<"-";
#endif
                    decoded=WVector(ld.taille);
                    for (int p=0;p<decoded.taille;p++) decoded.vect[p] = ld.vect[p]<0?1:0;
#ifdef VERBOSE2
                    cout <<"*";
#endif
                    nbbit+=decoded.taille;
                    nberr+=err = data.diff_count(decoded);
#ifdef ZIGZAG_FILE_VERBOSE    
                    os << snr << ":" << noise_variance<< ":" <<coded.taille<< ":" << decoded.taille<< ":" << err<< ":" << nberr<< ":" << nbbit<<endl;
#endif
                }
                os2 <<nb_encoder<<":"<< J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
#ifdef VERBOSE
                cout << J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
#endif
                if(nberr==0) break;
            }
        }
    }
    stoptime = time(NULL);
#ifdef ZIGZAG_FILE_VERBOSE    
    os  << "#stop time              = " << stoptime<<endl;
    os.close();
#endif
    os2 << "#stop time              = " << stoptime<<endl;
    os2.close();
    cout << endl << " * * * * DONE * * * *" << endl;

    cout << " Stopped at : " << stoptime << endl;
    cout << " Elapsed time (seconds) : " << stoptime - starttime << endl;
    //cin >> initpar;
    return -1;
}

// end direct concat



int main4(int argc, char ** argv) {
    int J = 4, I  = 127, initpar=0, multiplier=30;
    int blocksize = I*J, data_stream_size = multiplier*blocksize;
    int nb_data_stream = 10000, nberr=0,nbbit=0, err=0;
    unsigned starttime;//, stoptime;
    double noise_variance=0.0025, snr_min=0,snr_max=12, snr_step=.2, snr;
    WVector data(data_stream_size), coded, decoded;
    DVector noise, noisy, ld;
    WZigZag zzcoder(J+2,initpar);

    TRandBitStatePtr bitptr;
    starttime = bitptr.iseed = time(NULL);
    dRandGausStatePtr nptr;
    char filename1[100], filename2[100];
    sprintf(filename1, "zigzag_perf-%lu.log", bitptr.iseed);
    sprintf(filename2, "zigzag_perf_final-%lu.log", bitptr.iseed);
    ofstream os(filename1), os2(filename2);

    cout << " HELLO DARLING ... arrrrrrrrrrrrrrrrr .... " << endl;
    os << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os << "# J                = " << J << endl;
    os << "# I                = " << I << endl;
    os << "# initpar          = " << initpar << endl;
    os << "# multiplier       = " << multiplier << endl;
    os << "# blocksize        = " << blocksize << endl;
    os << "# data stream size = " << data_stream_size << endl;
    os << "# mn data stream   = " << nb_data_stream << endl;
    os << "# snr min          = " << snr_min << endl;
    os << "# snr max          = " << snr_max << endl;
    os << "# snr step         = " << snr_step << endl;
    os << "# **** End parameter *** " << endl;
    os << "# **** DATA **** " << endl;

    os << "snr:noise variance:coded.taille:decoded.taille:err:nberr:nbbit"<<endl;

    os2 << "# **** SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os2 << "# J                = " << J << endl;
    os2 << "# I                = " << I << endl;
    os2 << "# initpar          = " << initpar << endl;
    os2 << "# multiplier       = " << multiplier << endl;
    os2 << "# blocksize        = " << blocksize << endl;
    os2 << "# data stream size = " << data_stream_size << endl;
    os2 << "# mn data stream   = " << nb_data_stream << endl;
    os2 << "# snr min          = " << snr_min << endl;
    os2 << "# snr max          = " << snr_max << endl;
    os2 << "# snr step         = " << snr_step << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;

    os2 << "J:I:initpar:blocksize:snr:noise variance:nberr:nbbit:BER"<<endl;


    for(snr = snr_min; snr<=snr_max;snr+=snr_step){
        cout << endl << " Simulating for snr = " << snr << endl;
        noise_variance = convertsnrtosigma2(snr);
        dRandGausInit(&nptr,  0, noise_variance);
        nberr = nbbit = 0;
        for(int t=0;t<nb_data_stream;t++){
            bRandBit1(&bitptr, data);
#ifdef VERBOSE
            cout << "+";
#endif
            coded = zzcoder.encode(data, blocksize);
#ifdef VERBOSE
            cout <<"|";
#endif
            noise=DVector(coded.taille);
            dbRandGaus(&nptr, noise);
            noisy = noise+coded;
#ifdef VERBOSE
            cout <<"/";
#endif
            ld = zzcoder.decode(noisy);
#ifdef VERBOSE
            cout <<"-";
#endif
            decoded=WVector(ld.taille);
            for (int p=0;p<decoded.taille;p++) decoded.vect[p] = ld.vect[p]<0?1:0;
#ifdef VERBOSE
            cout <<"*";
#endif
            nbbit+=decoded.taille;
            nberr+=err = data.diff_count(decoded);
            os << snr << ":" << noise_variance<< ":" <<coded.taille<< ":" << decoded.taille<< ":" << err<< ":" << nberr<< ":" << nbbit<<endl;
        }
        os2 << J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
    }
    os.close();
    os.close();
    cout << endl << " * * * * DONE * * * *" << endl;
    cin >> blocksize;
    return -1;
}


int main3(int argc, char **argv) {
    int nbseg = 18, segsize = 9, initpar = 0;
    int blocksize=177;
    WVector data( 0+1*(177+(nbseg*(segsize-2))),1);
    srand(time(NULL));
    dRandGausStatePtr nptr;
    dRandGausInit(&nptr,  0, 0.025);

    for(int i = 0; i < data.taille;i++) data.vect[i]=  rand()%2;
//   cout << ReshapeLinewise(data, nbseg, segsize-2); cin >> blocksize;
    WVector out = zigzag_encode(data, segsize, initpar);
    CInterleaver interleaver1(out.taille, time(NULL));
    WVector i_out = interleaver1.Apply(out);
    show_zigzag_code(out, segsize);
    cout << endl << " Taille : " << out.taille << endl;
    //cout << out << i_out<<interleaver1;
    WZigZag zcode(segsize, initpar);
    WVector out2 = zcode.encode(data, blocksize);
    DVector noise(out2.taille), noisy;
    dbRandGaus(&nptr, noise);
    noisy = out2;//+noise;
#ifdef ZIGZAG_DEBUG
    ofstream os("decode_debug.out");
    os.precision(5);
    os <<data << out << out2 << noisy;
    os.close();
#endif
    noisy.vect[13] *= -1;
    DVector ld =  zcode.decode(noisy);
    WVector decoded(ld.taille);
    for (int p=0;p<decoded.taille;p++) decoded.vect[p] = ld.vect[p]<0?1:0;

    cout << data << decoded << endl <<"Diff with padding stripped is : " << data.diff_count(decoded.copy(0,data.taille));

    int test;
    cin >> test;

/*    WConcatZigZag concatenated_zigzag(2, segsize, initpar);
    WVector data3 = concatenated_zigzag.encode(data, blocksize);
    cout <<"*****"<<endl;
    cout << data << data3;
  */  


    return 1;
}

// all sweep
int main4_6(int argc, char ** argv) {
    int J, Jmin=32, Jmax=128, Jstep=4,
        I, Imin=64, Imax=64, Istep=4,
        nb_encoder, nb_encoder_min=7, nb_encoder_max=7, nb_encoder_step=3,
        initpar=0, multiplier=150;
    

    double noise_variance=0.0025, snr_min=3.0,snr_max=6.0, snr_step=0.5, snr;
    DVector noise, noisy, ld, llr;
    unsigned long starttime, stoptime;
    WVector data, coded, decoded, test;

    TRandBitStatePtr bitptr;
    starttime=bitptr.iseed = time(NULL);
    dRandGausStatePtr nptr;

    cout << " HELLO DARLING ... arrrrrrrrrrrrrrrrr .... " << starttime <<endl;

    char filename2[100];
    sprintf(filename2, "zigzag_perf_final-%lu.log", bitptr.iseed);
    ofstream os2(filename2);


    os2 << "# **** main4_6_"<<__TIME__<< "SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os2 << "# nb_encoder(min,max,step)       = (" << nb_encoder_min << "," << nb_encoder_max << "," << nb_encoder_step << ")" <<endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# initpar                        = " << initpar << endl;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# nb data stream                 = 100*nb_encoder"<< endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;

    os2 << "nb_encoder:J:I:initpar:blocksize:snr:noise variance:nberr:nbbit:BER"<<endl;
    
    for(nb_encoder=nb_encoder_min;nb_encoder<=nb_encoder_max;nb_encoder+=nb_encoder_step) {
        UVector<CInterleaver> interleaver(nb_encoder);
        UVector<WZigZag> encoder(nb_encoder);
        int nberr=0,nbbit=0, err=0;
    
        for(I=Imin; I<=Imax;I*=Istep){
            for(J=Jmin;J<=Jmax;J*=Jstep) {
                int nb_data_stream = 10*(Imax*Jmax/(I*J)) ;
                cout << "Nb of data stream : " << nb_data_stream<<endl;
                for(int x=0;x<nb_encoder;x++) {
                    encoder.vect[x] = WZigZag(J+2,initpar);
                    encoder.vect[x].concat_mode = 1;
                    encoder.vect[x].decoder_output_type = OUTPUT_TYPE_IS_BIT_APP;
                }
                encoder.vect[nb_encoder-1].concat_mode=0;
                encoder.vect[0].decoder_output_type = OUTPUT_TYPE_IS_BIT_APP;

                int blocksize = I*J, data_stream_size = multiplier*blocksize;
                data = WVector(data_stream_size);
                for(snr = snr_min; snr<=snr_max;snr+=snr_step){
                    double rate = J/(J+(double)nb_encoder);
                    noise_variance = convertsnrtosigma2(rate, snr);
                    cout << endl << " Simulating for nb encoder = " << nb_encoder << ", snr = " << snr << "@ J = " << J << ", I=" << I <<endl;
                    dRandGausInit(&nptr,  0, noise_variance);
                    nberr = nbbit = 0;
                    for(int t=0;t<nb_data_stream;t++){
                        bRandBit1(&bitptr, data);
                        coded = data;
                        for(int x = 0;x<nb_encoder;x++) {
                            interleaver.vect[x] = CInterleaver(coded.taille, starttime+x);
                            test = interleaver.vect[x].Apply(coded);
                            // problem with g++ v3.4.4 .. does not compile if line above is embedded in the next line
                            coded = encoder.vect[x].encode( test, blocksize);
                        }
                        noise=DVector(coded.taille,0.0);
                        dbRandGaus(&nptr, noise);
                        noisy = noise+coded;
                        ld = noisy;
                        for(int x=nb_encoder-1;x>=0;x--){
                            ld = interleaver.vect[x].Extract(encoder.vect[x].decode(ld, llr));
                            llr = interleaver.vect[x].Extract(llr);
                        }
    
                        decoded=WVector(ld.taille);
                        for (int p=0;p<decoded.taille;p++) decoded.vect[p] = ld.vect[p]<0?1:0;
                        nbbit+=decoded.taille;
                        nberr+=err = data.diff_count(decoded);
#ifdef ERROR_BASED_NB_SIM
                        if(nberr>ERROR_BASED_NB_SIM) break;
#endif                            
                    }
                    os2 <<nb_encoder<<":"<< J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
                    if(nberr==0) break;
                }
            }
        }
    }
    stoptime = time(NULL);
    os2 << "#stop time              = " << stoptime<<endl;
    os2.close();
    cout << endl << " * * * * DONE * * * *" << endl;

    cout << " Stopped at : " << stoptime << endl;
    cout << " Elapsed time (seconds) : " << stoptime - starttime << endl;
    //cin >> initpar;
    return -1;
}

// end direct concat

int main4_7(int argc, char ** argv) {
    int J, Jmin=8, Jmax=128, Jstep=2,
        I, Imin=16, Imax=16, Istep=4,
        nb_encoder, nb_encoder_min=1, nb_encoder_max=11, nb_encoder_step=3,
        initpar=0, multiplier=Imax*2;


    double noise_variance=0.0025, snr_min=0,snr_max=10.0, snr_step=0.5, snr;
    DVector noise, noisy, ld, llr;
    unsigned long starttime, stoptime;
    WVector data, coded, decoded, test;

    TRandBitStatePtr bitptr;
    starttime=bitptr.iseed = time(NULL);
    dRandGausStatePtr nptr;

    cout << " HELLO DARLING ... arrrrrrrrrrrrrrrrr .... " << starttime <<endl;

    char filename2[100];
    sprintf(filename2, "zigzag_perf_final-%lu.log", bitptr.iseed);
    ofstream os2(filename2);


    os2 << "# **** main4_6_"<<__TIME__<< "SIMULATION PARAMETER (runtime : "<< bitptr.iseed<<") ****"<< endl;
    os2 << "# Interleaver                    = NO" <<endl;
    os2 << "# nb_encoder(min,max,step)       = (" << nb_encoder_min << "," << nb_encoder_max << "," << nb_encoder_step << ")" <<endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# initpar                        = " << initpar << endl;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# nb data stream                 = 100*nb_encoder"<< endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;

    os2 << "nb_encoder:J:I:initpar:blocksize:snr:noise variance:nberr:nbbit:BER"<<endl;

    for(nb_encoder=nb_encoder_min;nb_encoder<=nb_encoder_max;nb_encoder+=nb_encoder_step) {
//        UVector<CInterleaver> interleaver(nb_encoder);
        UVector<WZigZag> encoder(nb_encoder);
        int nberr=0,nbbit=0, err=0;
    
        for(I=Imin; I<=Imax;I*=Istep){
            for(J=Jmin;J<=Jmax;J*=Jstep) {
                int nb_data_stream = 500*(Imax*Jmax/(I*J)) ;
                cout << "Nb of data stream : " << nb_data_stream<<endl;
                for(int x=0;x<nb_encoder;x++) {
                    encoder.vect[x] = WZigZag(J+2,initpar);
                    encoder.vect[x].concat_mode = 1;
                    encoder.vect[x].decoder_output_type = OUTPUT_TYPE_IS_BIT_APP;
                }
                encoder.vect[nb_encoder-1].concat_mode=0;
                encoder.vect[0].decoder_output_type = OUTPUT_TYPE_IS_BIT_APP;

                int blocksize = I*J, data_stream_size = multiplier*blocksize;
                data = WVector(data_stream_size);
                for(snr = snr_min; snr<=snr_max;snr+=snr_step){
                    double rate = J/(J+(double)nb_encoder);
                    noise_variance = convertsnrtosigma2(rate, snr);
                    cout << endl << " Simulating for nb encoder = " << nb_encoder << ", snr = " << snr << "@ J = " << J << ", I=" << I <<endl;
                    dRandGausInit(&nptr,  0, noise_variance);
                    nberr = nbbit = 0;
                    for(int t=0;t<nb_data_stream;t++){
                        bRandBit1(&bitptr, data);
                        coded = data;
                        for(int x = 0;x<nb_encoder;x++) {
  //                          interleaver.vect[x] = CInterleaver(coded.taille, starttime+x);
    //                        test = interleaver.vect[x].Apply(coded);
                            // problem with g++ v3.4.4 .. does not compile if line above is embedded in the next line
                            coded = encoder.vect[x].encode(coded, blocksize);
                        }
                        noise=DVector(coded.taille,0.0);
                        dbRandGaus(&nptr, noise);
                        noisy = noise+coded;
                        ld = noisy;
                        for(int x=nb_encoder-1;x>=0;x--){
//                            ld = interleaver.vect[x].Extract(encoder.vect[x].decode(ld, llr));
                            ld = encoder.vect[x].decode(ld, llr);
//                            llr = interleaver.vect[x].Extract(llr);
                        
                        }
    
                        decoded=WVector(ld.taille);
                        for (int p=0;p<decoded.taille;p++) decoded.vect[p] = ld.vect[p]<0?1:0;
                        nbbit+=decoded.taille;
                        nberr+=err = data.diff_count(decoded);
#ifdef ERROR_BASED_NB_SIM
                        if(nberr>ERROR_BASED_NB_SIM) break;
#endif
                    }
                    os2 <<nb_encoder<<":"<< J << ":" << I << ":" << initpar << ":" << blocksize << ":" << snr << ":" << noise_variance<< ":" << nberr<< ":" << nbbit<<":"<< nberr/(double)nbbit<<endl;
                    if(nberr==0) break;
                }
            }
        }
    }
    stoptime = time(NULL);
    os2 << "#stop time              = " << stoptime<<endl;
    os2.close();
    cout << endl << " * * * * DONE * * * *" << endl;

    cout << " Stopped at : " << stoptime << endl;
    cout << " Elapsed time (seconds) : " << stoptime - starttime << endl;
    //cin >> initpar;
    return -1;
}

// end direct concat

int pczz1(int argc, char **argv) {
    int K=4, I=7, J=16;
    WConcatZigZag pczz(K, I, J);
    WVector data(K*I*J+1,1);
    TRandBitStatePtr bitptr;
    bitptr.iseed=time(NULL);
//    bRandBit1(&bitptr, data);
    WMatrix output = ReshapeColumnwise(pczz.encode(data),(K+1)*I,J+K);
    cout << output;
  //  for(int i = 0;i<data.taille;i++) data.vect[i]=i;
   // cout << ReshapeColumnwise(data, K*I,J);
    cin >> K;
    return -1;


}

int vtest(){
    DVector noise(51124887);
    dRandGausStatePtr nptr;
    srand(time(NULL));
    for(int x=1;x<11;x++) {
        cout << " Sample " << x << endl;
        double variance = 5.0*rand()/(double)RAND_MAX;
        dRandGausInit(&nptr, (time(NULL)%255411) * (time(NULL)%411+x), 0.0, variance);
        dbRandGaus(&nptr, noise);
        cout << "Given variance : " << variance << endl << " Calculated : " << noise.variance() << endl;
        cout << " * * * "<< endl << endl;
    }
    int var;
    cout << " DONE" << endl;
    cin >>var;
    return 1;
}

WVector genprime(int nb) {
    /*WVector v(nb,0,'p');
    int i,x,y, l;
    bool isprime;

    if(nb<=0) return WVector();
    if(nb<4) return WVector(nb,1,'I');

    v.vect[0]=1;v.vect[1]=2;v.vect[2]=3;
    cout << endl << "Initializing prime vector with 1,2 and 3" << endl;
    for(i=5,x=3;x<nb;i+=2) {
        l = (int) sqrt(i)+1;
        isprime=true;

        for(y=1;y<x && v.vect[y]<=l;y++) {
            if(i%v.vect[y]==0){ isprime=false; break; }
        }
        if(isprime) {
            v.vect[x++]=i;
//            cout << "."; // Found " << i << " as a prime number ... (found " << x << " out of " << nb << ")" << endl;
        }
    }
    return v;
    */
    return WVector(nb,0,'p');
}

int itest() {
    const int VSIZE=9900000;
    WVector v(VSIZE, 0, 'd');
    unsigned long int starttime=time(NULL), stoptime;
    CInterleaver i(VSIZE, 2501);
    stoptime=time(NULL);
    cout <<endl << "Elapsed time ... " << stoptime-starttime << " seconds" << endl;


    cout << "done ";
    char makan;
    cin >> makan;
    return 0;
}


int pconcat_zigzag_test(){
    int I=310,J=5,E=4, itermax=5;
    WVector bitin(2*I*J, 77, 'b');
//    cout << bitin << ReshapeLinewise(bitin,I,J) << enc;
//    enc = 1-2*enc;
    dRandGausStatePtr nptr;
    for(double s=0.0;s<=3.04;s+=0.5) {
        dRandGausInit(&nptr, 777, 0, convertsnrtosigma2(s));
        for(int e=E;e<=E;e++) {
            WConcatZigZag pczz(e,I,J);
//            pczz.printInterleaver();
            WVector enc = 1-2*pczz.encode(bitin);
            DVector noise(enc.taille);
            dbRandGaus(&nptr, noise);
            for(int i = 1;i<=itermax;i+=1) {
                noise+=enc;
                DVector dec = pczz.decode(noise, i);
                WVector out(dec.taille);
                for(int v=0;v<out.taille;v++) out.vect[v]=dec.vect[v]<=0.0 ? 1 : 0;
                int diff = out.diff_count(bitin);
                printf("OI > (SNR=%.1fdB, encoder=%d, iteration=%2i) %d differences found out of %d bit (BER=%.4f)\r\n", s,e, i,diff,bitin.taille,diff/(double)bitin.taille);
                if(diff==0) break;
            }
        }
    }

    return 1;
}

int apconcat_zigzag_test(){
    int I=10,J=5,E=4,K=31, itermax=5;
    WVector bitin(2*I*J*K, 77, 'b');
//    cout << bitin << ReshapeLinewise(bitin,I,J) << enc;
//    enc = 1-2*enc;
    dRandGausStatePtr nptr;
    for(double s=0.0;s<=3.04;s+=0.5) {
        dRandGausInit(&nptr, 777, 0, convertsnrtosigma2(s));
        for(int e=E;e<=E;e++) {
            WAConcatZigZag pczz(e,I,J,K);
//            pczz.printInterleaver();
            WVector enc =1-2*pczz.encode(bitin);
//            cout << bitin <<enc;
            DVector noise(enc.taille, 0.0);
            dbRandGaus(&nptr, noise);
            for(int i = 1;i<=itermax;i+=1) {
                noise+=enc;
                DVector dec = pczz.decode(noise, i);
                WVector out(dec.taille);
                for(int v=0;v<out.taille;v++) out.vect[v]=dec.vect[v]<=0.0 ? 1 : 0;
                int diff = out.diff_count(bitin);
                printf("AI > (SNR=%.1fdB, encoder=%d, iteration=%2i) %d differences found out of %d bit (BER=%.4f)\r\n", s,e, i,diff,bitin.taille,diff/(double)bitin.taille);
                if(diff==0) break;
            }
        }

    }

    return 1;
}

int pconcat_zigzag_test2(){
    int j, Jmin=64, Jmax=64, Jstep=8,
        i, Imin=512, Imax=4092, Istep=2,
        itmin=2, itmax=3, itstep=1,
        nb_encoder_min=4, nb_encoder_max=5, nb_encoder_step=3,
        multiplier=10;
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    int nbbit_min=10000000, nbbit_max=500000000;
#else
     int nbblock=1000;
#endif
    double noise_variance=0.0025, snr_min=3.0,snr_max=7.01, snr_step=.3;
    i=0;
    WVector iteration(6);
    iteration.vect[i++] = 1;
    iteration.vect[i++] = 2;
    iteration.vect[i++] = 7;
    iteration.vect[i++] = 12;
    iteration.vect[i++] = 17;
    iteration.vect[i++] = 22;

    unsigned long starttime, stoptime;

    starttime=time(NULL);
    dRandGausStatePtr nptr;
    char filename2[255];

    sprintf(filename2, "pczz_perf_final-%lu.log", starttime);
    ofstream os2(filename2);


    os2 << "# **** pconcat_zigzag_test2_"<<__TIME__<< "SIMULATION PARAMETER (runtime : "<< starttime<<") ****"<< endl;
    os2 << "# nb_encoder(min,max,step)       = (" << nb_encoder_min << "," << nb_encoder_max << "," << nb_encoder_step << ")" <<endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# Iteration                      : " << iteration;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# nb data stream                 = 100*nb_encoder"<< endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    os2 << "# Full monte carlo enabled " << endl;
    os2 << "# Required minimum error                  = " << MIN_ERROR_REQUIRED<< endl;
    os2 << "# Required minimum bit simulated        = " << nbbit_min << endl;
    os2 << "# Safeguard maximum bit simulated       = " << nbbit_max << endl;
#endif
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;
    os2 << "nb_encoder;J;I;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;
    cout << "nb_encoder;J;I;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;


    WVector enc;
    DVector noisy, dec;
    for(i=Imin;i<=Imax;i*=Istep) { // 8 16 32 64 128    [5]
        for(j=Jmin;j<=Jmax;j+=Jstep) { // 8 16 32 64 128   [5]
            for(int e=nb_encoder_min; e<=nb_encoder_max;e+=nb_encoder_step) { //1 4 7 11 [4]
                WConcatZigZag pczz(e,i,j);
//                for(int it=itmin;it<=itmax;it+=itstep) { // 1 6 11 16 21 [5]
                for(int it=itmin;it<itmax;it+=itstep) { // 1 6 11 16 21 [5]
                    for(double s=snr_min; s<=snr_max; s+=snr_step) { //   0..10 [11]
                        noise_variance = convertsnrtosigma2((j/((double)j+e)), s);
                        dRandGausInit(&nptr, 0.0, noise_variance);
                        int nbbit=0,diff=0;
#ifndef SIM_UNTIL_MIN_ERROR_ACHIEVED
                        for(int n=0;n<nbblock;n++) { // [500]
#else
                        for(int n=0;(nbbit<nbbit_min || diff<MIN_ERROR_REQUIRED)&&nbbit<nbbit_max;n++) { // [500]
#endif
                            WVector bitin(multiplier*i*j, starttime+i*j+j*e +i*e + it*207 + n*1997, 'b');
                            enc = 1-2*pczz.encode(bitin);
                            DVector noise(enc.taille);
                            dbRandGaus(&nptr, noise);
                            noise += enc;
                            dec = pczz.decode(noise,iteration.vect[it]);
                            WVector out=WVector(dec.taille);
                            for(int x=0;x<out.taille;x++) out.vect[x]=dec.vect[x]<=0.0 ? 1 : 0;
                            diff += out.diff_count(bitin);
                            nbbit+= out.taille;
                            if(diff>=ERROR_BASED_NB_SIM) break;
                        }
                        os2 << e << ";" << j << ";" << i << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << diff << ";" << nbbit << ";" << diff/(double)nbbit << endl;
                        cout << e << ";" << j << ";" << i << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << diff << ";" << nbbit << ";" << diff/(double)nbbit << endl;
                        if(diff<MIN_ERROR_REQUIRED) break;
                    }
                }
            }
        }
    }

    os2.close();
    stoptime=time(NULL);
    cout << "Done \r\nElapsed time : " << stoptime - starttime << " seconds " << endl;

    return 1;
}

int pconcat_zigzag_test3(){
    int j=0, Jmin=32, Jmax=128, Jstep=2,
        i, Imin=128, Imax=4097, Istep=2,
        itmin=3, itmax=4, itstep=1,
        nb_encoder_min=2, nb_encoder_max=9, nb_encoder_step=3,
        multiplier=10;
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    int nbbit_min=10000000, nbbit_max=100000000;
#else
     int nbblock=1000    
#endif
    double noise_variance=0.0025, snr_min=1,snr_max=6.01, snr_step=.4;
    i=0;
    WVector iteration(6);
    iteration.vect[i++] = 1;
    iteration.vect[i++] = 2;
    iteration.vect[i++] = 7;
    iteration.vect[i++] = 12;
    iteration.vect[i++] = 17;
    iteration.vect[i++] = 22;

    unsigned long starttime, stoptime;

    starttime=time(NULL);
    dRandGausStatePtr nptr;
    char filename2[255];

    sprintf(filename2, "pczz_perf_final-%lu.log", starttime);
    ofstream os2(filename2);


    os2 << "# **** pconcat_zigzag_test3_"<<__TIME__<< "SIMULATION PARAMETER (runtime : "<< starttime<<") ****"<< endl;
    os2 << "# nb_encoder(min,max,step)       = (" << nb_encoder_min << "," << nb_encoder_max << "," << nb_encoder_step << ")" <<endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# Iteration                      : " << iteration;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# nb data stream                 = 100*nb_encoder"<< endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    os2 << "# Full monte carlo enabled " << endl;
    os2 << "# Required minimum error                  = " << MIN_ERROR_REQUIRED<< endl;
    os2 << "# Required minimum bit simulated        = " << nbbit_min << endl;
    os2 << "# Safeguard maximum bit simulated       = " << nbbit_max << endl;
#endif
    os2 << "# Using optimized interleaver           =  YES " << nbbit_max << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;
    os2 << "nb_encoder;J;I;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;
    cout << "nb_encoder;J;I;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;


    WVector enc;
    DVector noisy, dec;
    for(i=Imin;i<=Imax;i*=Istep) { // 8 16 32 64 128    [5]
        for(j=Jmin;j<=Jmax;j*=Jstep) { // 8 16 32 64 128   [5]
            for(int e=nb_encoder_min; e<=nb_encoder_max;e+=nb_encoder_step) { //1 4 7 11 [4]
                try {
                    WConcatZigZag pczz('o', e,i,j);
//                for(int it=itmin;it<=itmax;it+=itstep) { // 1 6 11 16 21 [5]
                    for(int it=itmin;it<itmax;it+=itstep) { // 1 6 11 16 21 [5]
                        for(double s=snr_min; s<=snr_max; s+=snr_step) { //   0..10 [11]
                            noise_variance = convertsnrtosigma2((j/((double)j+e)), s);
                            dRandGausInit(&nptr, 0.0, noise_variance);
                            int nbbit=0,diff=0;
#ifndef SIM_UNTIL_MIN_ERROR_ACHIEVED
                            for(int n=0;n<nbblock;n++) { // [500]
#else
                            for(int n=0;(nbbit<nbbit_min || diff<MIN_ERROR_REQUIRED)&&nbbit<nbbit_max;n++) { // [500]
#endif
                                WVector bitin(multiplier*i*j, starttime+i*j+j*e +i*e + it*207 + n*1997, 'b');
                                enc = 1-2*pczz.encode(bitin);
                                DVector noise(enc.taille);
                                dbRandGaus(&nptr, noise);
                                noise += enc;
                                dec = pczz.decode(noise,iteration.vect[it]);
                                WVector out=WVector(dec.taille);
                                for(   int x=0;x<out.taille;x++) out.vect[x]=dec.vect[x]<=0.0 ? 1 : 0;
                                diff += out.diff_count(bitin);
                                    nbbit+= out.taille;
                                if(diff>=ERROR_BASED_NB_SIM) break;
                            }
                            os2 << e << ";" << j << ";" << i << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << diff << ";" << nbbit << ";" << diff/(double)nbbit << endl;
                            cout << e << ";" << j << ";" << i << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << diff << ";" << nbbit << ";" << diff/(double)nbbit << endl;
                            if(diff<MIN_ERROR_REQUIRED) break;
                        }
                    }
                } catch (...) {
                    cerr << endl << " nbenc failed ... " << endl;
                    break;
                }
            }
        }
    }

    os2.close();
    stoptime=time(NULL);
    cout << "Done \r\nElapsed time : " << stoptime - starttime << " seconds " << endl;

    return 1;
}

int allpzzcc_test1(){
    int j=0, Jmin=32, Jmax=65, Jstep=2,
        i, Imin=5, Imax=11, Istep=5,
        k=0, Kmin=100, Kmax = 101, Kstep = 2,
        itmin=3, itmax=4, itstep=1,
        nb_encoder_min=8, nb_encoder_max=9, nb_encoder_step=3,
        multiplier=10;
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    int nbbit_min=10000000, nbbit_max=500000000;
#else
     int nbblock=1000    
#endif
    double noise_variance=0.0025, snr_min=1.4,snr_max=6.21, snr_step=.25;
    i=0;
    WVector iteration(6);
    iteration.vect[i++] = 1;
    iteration.vect[i++] = 2;
    iteration.vect[i++] = 4;
    iteration.vect[i++] = 7;
    iteration.vect[i++] = 17;
    iteration.vect[i++] = 22;

    unsigned long starttime, stoptime;

    starttime=time(NULL);
    dRandGausStatePtr nptr;
    char filename2[255];

    sprintf(filename2, "pczz_perf_final-%lu.log", starttime);
    ofstream os2(filename2);


    os2 << "# **** allpzzcc_test1_"<<__TIME__<< "SIMULATION PARAMETER (runtime : "<< starttime<<") ****"<< endl;
    os2 << "# nb_encoder(min,max,step)       = (" << nb_encoder_min << "," << nb_encoder_max << "," << nb_encoder_step << ")" <<endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# K(min,max,step)                = (" << Kmin << "," << Kmax << "," << Kstep << ")" <<endl;
    os2 << "# Iteration                      : " << iteration;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# nb data stream                 = 100*nb_encoder"<< endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    os2 << "# Full monte carlo enabled " << endl;
    os2 << "# Required minimum error                  = " << MIN_ERROR_REQUIRED<< endl;
    os2 << "# Required minimum bit simulated        = " << nbbit_min << endl;
    os2 << "# Safeguard maximum bit simulated       = " << nbbit_max << endl;
#endif
    os2 << "# Using optimized interleaver           =  o,a=YES, r=NO " << nbbit_max << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;
    os2 << "interleaver;nb_encoder;J;I;K;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;
    cout << "interleaver;nb_encoder;J;I;K;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;


    WVector /*oenc,*/ renc, aenc;
    DVector noisy, /*odec,*/rdec,adec;
    for(i=Imin;i<=Imax;i+=Istep) { // 8 16 32 64 128    [5]
        for(k=Kmin;k<Kmax;k*=Kstep) {
            for(j=Jmin;j<=Jmax;j*=Jstep) { // 8 16 32 64 128   [5]
                for(int e=nb_encoder_min; e<=nb_encoder_max;e+=nb_encoder_step) { //1 4 7 11 [4]
                    try {
//                        WConcatZigZag opczz('o', e,i*k,j);
                        WAConcatZigZag apczz(e,i,j,k);
                        WConcatZigZag rpczz(e,i*k,j);
    //                for(int it=itmin;it<=itmax;it+=itstep) { // 1 6 11 16 21 [5]
                        for(int it=itmin;it<itmax;it+=itstep) { // 1 6 11 16 21 [5]
                            for(double s=snr_min; s<=snr_max; s+=snr_step) { //   0..10 [11]
                                noise_variance = convertsnrtosigma2((j/((double)j+e)), s);
                                dRandGausInit(&nptr, 0.0, noise_variance);
                                int /*onbbit=0,odiff=0,*/rnbbit=0,rdiff=0, anbbit=0, adiff=0;
#ifndef SIM_UNTIL_MIN_ERROR_ACHIEVED
                                for(int n=0;n<nbblock;n++) { // [500]
#else
                                for(int n=0;(/*onbbit<nbbit_min || */rdiff<MIN_ERROR_REQUIRED)&& rnbbit<nbbit_max;n++) { // [500]
#endif
                                    WVector bitin(multiplier*i*j*k, starttime+i*j+j*e +i*e + it*207 + n*1997+k*18*i, 'b');
                            ///        oenc = 1-2*opczz.encode(bitin);
                                    renc = 1-2*rpczz.encode(bitin);
                                    aenc = 1-2*apczz.encode(bitin);
                                   // DVector onoise(oenc.taille);
// //                                   dbRandGaus(&nptr, onoise);
     //                               noisy =onoise+ oenc;
//                                    odec = opczz.decode(noisy,iteration.vect[it]);
  //                                  WVector oout=WVector(odec.taille);
    //                                for(int x=0;x<oout.taille;x++) oout.vect[x]=odec.vect[x]<=0.0 ? 1 : 0;

                                    DVector rnoise(renc.taille);
                                    dbRandGaus(&nptr, rnoise);
                                    noisy = rnoise+renc;
                                    rdec = rpczz.decode(noisy,iteration.vect[it]);
                                    WVector rout=WVector(rdec.taille);
                                    for(int x=0;x<rout.taille;x++) rout.vect[x]=rdec.vect[x]<=0.0 ? 1 : 0;

                                    DVector anoise(aenc.taille);
                                    dbRandGaus(&nptr, anoise);
                                    noisy = anoise+aenc;
                                    adec = apczz.decode(noisy,iteration.vect[it]);
                                    WVector aout=WVector(adec.taille);
                                    for(int x=0;x<aout.taille;x++) aout.vect[x]=adec.vect[x]<=0.0 ? 1 : 0;

      //                              odiff += oout.diff_count(bitin);
        //                            onbbit+= oout.taille;

                                    rdiff += rout.diff_count(bitin);
                                    rnbbit+= rout.taille;

                                    adiff += aout.diff_count(bitin);
                                    anbbit+= aout.taille;

                                    if(/*odiff>=ERROR_BASED_NB_SIM && */rdiff>=ERROR_BASED_NB_SIM && adiff>=ERROR_BASED_NB_SIM) break;
                                }
//                                os2 << "o;" << e << ";" << j << ";" << i << ";" << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << odiff << ";" << onbbit << ";" << odiff/(double)onbbit << endl;
  //                              cout << "o;" << e << ";" << j << ";" << i << ";" << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << odiff << ";" << onbbit << ";" << odiff/(double)onbbit << endl;
                                os2 << "r;" << e << ";" << j << ";" << i << ";" << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << rdiff << ";" << rnbbit << ";" << rdiff/(double)rnbbit << endl;
                                cout << "r;" << e << ";" << j << ";" << i << ";" << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << rdiff << ";" << rnbbit << ";" << rdiff/(double)rnbbit << endl;
                                os2 << "a;" << e << ";" << j << ";" << i << ";" << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << adiff << ";" << anbbit << ";" << adiff/(double)anbbit << endl;
                                cout << "a;" << e << ";" << j << ";" << i << ";" << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << adiff << ";" << anbbit << ";" << adiff/(double)anbbit << endl;
                                if(/*odiff<MIN_ERROR_REQUIRED && */rdiff<MIN_ERROR_REQUIRED) break;
                            }
                        }
                    } catch (...) {
                            cerr << endl << " nbenc failed ... (" << i << ";" <<j << ";" << k << ")" << endl;
                        break;
                    }
                }
            }
        }
    }

    os2.close();
    stoptime=time(NULL);
    cout << "Done \r\nElapsed time : " << stoptime - starttime << " seconds " << endl;

    return 1;
}

int allpzzcc_test2(){
    int j=0, Jmin=32, Jmax=129, Jstep=15,
        i, Imin=5, Imax=11, Istep=5,
        k=0, Kmin=25, Kmax = 101, Kstep = 2,
        itmin=2, itmax=4, itstep=1,
        nb_encoder_min=5, nb_encoder_max=9, nb_encoder_step=3,
        multiplier=10;
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    int nbbit_min=10000000, nbbit_max=100000000;
#else
     int nbblock=1000
#endif
    double noise_variance=0.0025, snr_min=1.4,snr_max=6.21, snr_step=.25;
    i=0;
    WVector iteration(6);
    iteration.vect[i++] = 1;
    iteration.vect[i++] = 2;
    iteration.vect[i++] = 4;
    iteration.vect[i++] = 7;
    iteration.vect[i++] = 17;
    iteration.vect[i++] = 22;

    unsigned long starttime, stoptime;

    starttime=time(NULL);
    dRandGausStatePtr nptr;
    char filename2[255];

    sprintf(filename2, "pczz_perf_final-%lu.log", starttime);
    ofstream os2(filename2);


    os2 << "# **** allpzzcc_test1_"<<__TIME__<< "SIMULATION PARAMETER (runtime : "<< starttime<<") ****"<< endl;
    os2 << "# nb_encoder(min,max,step)       = (" << nb_encoder_min << "," << nb_encoder_max << "," << nb_encoder_step << ")" <<endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# K(min,max,step)                = (" << Kmin << "," << Kmax << "," << Kstep << ")" <<endl;
    os2 << "# Iteration                      : " << iteration;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# nb data stream                 = 100*nb_encoder"<< endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    os2 << "# Full monte carlo enabled " << endl;
    os2 << "# Required minimum error                  = " << MIN_ERROR_REQUIRED<< endl;
    os2 << "# Required minimum bit simulated        = " << nbbit_min << endl;
    os2 << "# Safeguard maximum bit simulated       = " << nbbit_max << endl;
#endif
    os2 << "# Using optimized interleaver           =  o,a=YES, r=NO " << nbbit_max << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;
    os2 << "interleaver;nb_encoder;J;I;K;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;
    cout << "interleaver;nb_encoder;J;I;K;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;


    WVector /*oenc, renc, */aenc;
    DVector noisy, /*odec,rdec,*/adec;
    for(i=Imin;i<=Imax;i+=Istep) { // 8 16 32 64 128    [5]
        for(k=Kmin;k<Kmax;k*=Kstep) {
            for(j=Jmin;j<=Jmax;j+=Jstep) { // 8 16 32 64 128   [5]
                for(int e=nb_encoder_min; e<=nb_encoder_max;e+=nb_encoder_step) { //1 4 7 11 [4]
                    try {
                        WAConcatZigZag apczz(e,i,j,k);
                        for(int it=itmin;it<itmax;it+=itstep) { // 1 6 11 16 21 [5]
                            for(double s=snr_min; s<=snr_max; s+=snr_step) { //   0..10 [11]

                                noise_variance = convertsnrtosigma2((j/((double)j+e)), s);
                                dRandGausInit(&nptr, 0.0, noise_variance);
                                int anbbit=0, adiff=0;
#ifndef SIM_UNTIL_MIN_ERROR_ACHIEVED
                                for(int n=0;n<nbblock;n++) { // [500]
#else
                                for(int n=0;(adiff<MIN_ERROR_REQUIRED && anbbit<nbbit_max)  ;n++) { // [500]
#endif
#ifdef PROGRESS_WATCH
                                      printf("SNR:%f, i=%d, j=%d, k=%d, e=%d, it=%d, n=%d\r\n", s, i, j, k, e, it,n);
#endif                                      

                                    WVector bitin(multiplier*i*j*k, starttime+i*j+j*e +i*e + it*207 + n*1997+k*18*i, 'b');
#ifdef PROGRESS_WATCH
                                    cout << __FILE__ <<"@"<<__LINE__<<endl;
#endif                                      
                                    aenc = 1-2*apczz.encode(bitin);
#ifdef PROGRESS_WATCH
                                    cout << __FILE__ <<"@"<<__LINE__<<endl;
#endif                                      
                                    DVector anoise(aenc.taille);
                                    dbRandGaus(&nptr, anoise);
                                    noisy = anoise+aenc;
#ifdef PROGRESS_WATCH
                                    cout << __FILE__ <<"@"<<__LINE__<<endl;
#endif                                      
                                    adec = apczz.decode(noisy,iteration.vect[it]);
#ifdef PROGRESS_WATCH
                                    cout << __FILE__ <<"@"<<__LINE__<<endl;
#endif                                      
                                    WVector aout=WVector(adec.taille);
                                    for(int x=0;x<aout.taille;x++) aout.vect[x]=adec.vect[x]<=0.0 ? 1 : 0;
#ifdef PROGRESS_WATCH
                                    cout << __FILE__ <<"@"<<__LINE__<<endl;
#endif                                      

                                    adiff += aout.diff_count(bitin);
                                    anbbit+= aout.taille;

                                    if(adiff>=ERROR_BASED_NB_SIM) break;
                                }
                                os2 << "a;" << e << ";" << j << ";" << i << ";" << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << adiff << ";" << anbbit << ";" << adiff/(double)anbbit << endl;
                                cout << "a;" << e << ";" << j << ";" << i << ";" << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << adiff << ";" << anbbit << ";" << adiff/(double)anbbit << endl;
                                if(adiff<MIN_ERROR_REQUIRED) break;
                            }
                        }
                    } catch (...) {
                            cerr << endl << " nbenc failed ... (" << i << ";" <<j << ";" << k << ")" << endl;
                        break;
                    }
                }
            }
        }
    }

    os2.close();
    stoptime=time(NULL);
    cout << "Done \r\nElapsed time : " << stoptime - starttime << " seconds " << endl;

    return 1;
}



int aconcat_zigzag_test1(){
    int j=0, Jmin=32, Jmax=128, Jstep=2,
        i=0, Imin=5, Imax=30, Istep=5,
        k=0, Kmin=1, Kmax = 17, Kstep = 2,
        itmin=3, itmax=4, itstep=1,
        nb_encoder_min=2, nb_encoder_max=9, nb_encoder_step=3,
        multiplier=10; 
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    int nbbit_min=10000000, nbbit_max=100000000;
#else
     int nbblock=1000;    
#endif
    double noise_variance=0.0025, snr_min=1,snr_max=6.01, snr_step=.4;
    i=0;
    WVector iteration(6);
    iteration.vect[i++] = 1;
    iteration.vect[i++] = 2;
    iteration.vect[i++] = 7;
    iteration.vect[i++] = 12;
    iteration.vect[i++] = 17;
    iteration.vect[i++] = 22;

    unsigned long starttime, stoptime;

    starttime=time(NULL);
    dRandGausStatePtr nptr;
    char filename2[255];

    sprintf(filename2, "pczz_perf_final-%lu.log", starttime);
    ofstream os2(filename2);


    os2 << "# **** aconcat_zigzag_test1_"<<__TIME__<< "SIMULATION PARAMETER (runtime : "<< starttime<<") ****"<< endl;
    os2 << "# nb_encoder(min,max,step)       = (" << nb_encoder_min << "," << nb_encoder_max << "," << nb_encoder_step << ")" <<endl;
    os2 << "# J(min,max,step)                = (" << Jmin << "," << Jmax << "," << Jstep << ")" <<endl;
    os2 << "# I(min,max,step)                = (" << Imin << "," << Imax << "," << Istep << ")" <<endl;
    os2 << "# K(min,max,step)                = (" << Kmin << "," << Kmax << "," << Kstep << ")" <<endl;
    os2 << "# Iteration                      : " << iteration;
    os2 << "# multiplier                     = " << multiplier << endl;
    os2 << "# blocksize(min,max)             = (" << Imin*Jmin << "," << Imax*Jmax << ")" <<endl;
    os2 << "# data stream size(min,max) = (" << multiplier*Imin*Jmin << "," << multiplier*Imax*Jmax << ")" <<endl;
    os2 << "# nb data stream                 = 100*nb_encoder"<< endl;
    os2 << "# snr min                        = " << snr_min << endl;
    os2 << "# snr max                        = " << snr_max << endl;
    os2 << "# snr step                       = " << snr_step << endl;
#ifdef SIM_UNTIL_MIN_ERROR_ACHIEVED
    os2 << "# Full monte carlo enabled " << endl;
    os2 << "# Required minimum error                  = " << MIN_ERROR_REQUIRED<< endl;
    os2 << "# Required minimum bit simulated        = " << nbbit_min << endl;
    os2 << "# Safeguard maximum bit simulated       = " << nbbit_max << endl;
#endif
    os2 << "# Using optimized interleaver           =  NO " << endl;
    os2 << "# Using asynchronous concatenated zigzag code =  YES " << endl;
    os2 << "# **** End parameter *** " << endl;
    os2 << "# **** DATA **** " << endl;
    os2 << "nb_encoder;J;I;K;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;
    cout << "nb_encoder;J;I;K;iteration;snr;noise variance;nberr;nbbit;BER"<<endl;


    WVector enc;
    DVector noisy, dec;
    for(i=Imin;i<=Imax;i+=Istep) { // 8 16 32 64 128    [5]
        for(k=Kmin;k<Kmax;k+=Kstep){
            for(j=Jmin;j<=Jmax;j*=Jstep) { // 8 16 32 64 128   [5]
                for(int e=nb_encoder_min; e<=nb_encoder_max;e+=nb_encoder_step) { //1 4 7 11 [4]
                    try {
                        WAConcatZigZag pczz( e,i,j,k);
    //                for(int it=itmin;it<=itmax;it+=itstep) { // 1 6 11 16 21 [5]
                        for(int it=itmin;it<itmax;it+=itstep) { // 1 6 11 16 21 [5]
                            for(double s=snr_min; s<=snr_max; s+=snr_step) { //   0..10 [11]
                                noise_variance = convertsnrtosigma2((j/((double)j+e)), s);
                                dRandGausInit(&nptr, 0.0, noise_variance);
                                int nbbit=0,diff=0;
#ifndef SIM_UNTIL_MIN_ERROR_ACHIEVED
                                for(int n=0;n<nbblock;n++) { // [500]
#else
                                for(int n=0;(nbbit<nbbit_min || diff<MIN_ERROR_REQUIRED)&&nbbit<nbbit_max;n++) { // [500]
#endif
                                    WVector bitin(multiplier*i*j, starttime+i*j+j*e +i*e + it*207 + n*1997, 'b');
                                    enc = 1-2*pczz.encode(bitin);
                                    DVector noise(enc.taille);
                                    dbRandGaus(&nptr, noise);
                                    noise += enc;
                                    dec = pczz.decode(noise,iteration.vect[it]);
                                    WVector out=WVector(dec.taille);
                                    for(   int x=0;x<out.taille;x++) out.vect[x]=dec.vect[x]<=0.0 ? 1 : 0;
                                    diff += out.diff_count(bitin);
                                        nbbit+= out.taille;
                                    if(diff>=ERROR_BASED_NB_SIM) break;
                                }
                                os2 << e << ";" << j << ";" << i << ";"  << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << diff << ";" << nbbit << ";" << diff/(double)nbbit << endl;
                                cout << e << ";" << j << ";" << i << ";"  << k << ";" << iteration.vect[it] << ";" << s << ";" << noise_variance << ";" << diff << ";" << nbbit << ";" << diff/(double)nbbit << endl;
                                if(diff<MIN_ERROR_REQUIRED) break;
                            }
                        }
                    } catch (...) {
                        cerr << endl << " nbenc failed ... (" << i << ";" <<j << ";" << k << ")" << endl;
                        break;
                    }
                }
            }
        }
    }

    os2.close();
    stoptime=time(NULL);
    cout << "Done \r\nElapsed time : " << stoptime - starttime << " seconds " << endl;

    return 1;
}


void OptimizedInterleaver1(int I, int J){
    int a=-1,b=-1, c, /*as,bs,*/ ap,bp,/* is, js, */vmax, i2, imin;
    WMatrix m(I,J); WVector v(I,0,'i');

/*    for(int i=0,c=0;i<I;i++) {
        for(int j=0;j<J;j++){
            m.mat[i][j]=c++;
        }
    }             */

    //find b .. b must be relative prime to I and b = floor(I/J)
    vmax = I/J;
    for(c=vmax;c>1;c--) {
        if (isRelativePrime(c,I)) {b = c; break;}
    }
    if(b==-1) throw CZigZagException(CREATE_OPTIMIZED_INTERLEAVER_FAILED, __LINE__, __FILE__);
    //have b ready ..  get a now
    vmax = I-1;
    for(c=vmax;c>1;c--) {
        if (isRelativePrime(c,I) && isRelativePrime(c,b)) {a = c; break;}
    }
    if(a==-1) throw CZigZagException(CREATE_OPTIMIZED_INTERLEAVER_FAILED, __LINE__, __FILE__);
    // have a and b .. calculate ap and bp

    ap = a%I;
    bp = b%I;

    // get the first column
    for(int i=0;i<I;i++) {
        m.mat[i][0]=(v.vect[i]=((a*i)%I))*J;
    }

    imin=-1;

    for(i2=0;i2<I;i2++)
        if (((a*i2)%I)==b) {
            imin =i2;
            break;
        }

    if(imin==-1) throw CZigZagException(CREATE_OPTIMIZED_INTERLEAVER_FAILED, __LINE__, __FILE__);
            // have imin .. shift the first column imin*j up ..

    cout << v;

    int offset;

    for(int j=1;j<J;j++){
        offset=imin*j;
        for(i2=0;i2<I;i2++) {
            m.mat[i2][j] = v.vect[(J*I+i2-offset)%I]*J+j;
            cout << m;
        }
    }

    cout << m;
    WVector map = MatrixInVectorRepresentation(ReshapeLinewise(m,1,I*J));
    map.vect[17]=0;

    CInterleaver ti(map);
    cout << map << ti;
}

void opczz_test() {
    WConcatZigZag opczz('o', 16, 25,4);
}

int main(int argc, char ** argv) {
//    cout << genprime(2000000).vect[1999999];
//    pconcat_zigzag_test3();
//    aconcat_zigzag_test1();
      allpzzcc_test2();
//    OptimizedInterleaver1(23,4);
//    opczz_test();
//    apconcat_zigzag_test();
  //  pconcat_zigzag_test();
    int x;
    cin >> x;
 //   return main4_6(argc,argv);
     return 1;
}
