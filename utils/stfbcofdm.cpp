/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/12/08 08:36:32 $
 * $Revision: 1.1.2.3 $
 * $Id: stfbcofdm.cpp,v 1.1.2.3 2004/12/08 08:36:32 syed Exp $
 *******************************************************************/

#include <all.h>


int main() {
    TConstellation constel=psk8;
    int nb_ofdm=16000, nb_porteuse=128;
    int timestamp=time(NULL);
    int bsize=constel.bit_per_sym*nb_ofdm*nb_porteuse*4;
    double snr, snr_min=0, snr_max=25, snr_step=3;
    int snr_index; //, rxcnt;
    WVector bitin(bsize);
    TRandBitStatePtr bseed;
    bseed.iseed=time(NULL)+314159;
    srand(timestamp);
    bRandBit1(&bseed, bitin);
    double sqrt_05=sqrt(0.5);
    cout << "Modulating ... " <<endl;
    ZVector symtx = (*constel.modulate_function)(bitin)*sqrt_05;
    zRandGausStatePtr nptr;
    dRandUniStatePtr tseed, pseed;
    dRandUniInit(&tseed);
    dRandUniInit(&pseed, rand()%8872411);
    double max_doppler = 75;
    double symbol_time = 64;
    int filter_length=4;
//    ZMatrix xpath;
    char varname[450];
    ifstream pathfile[2];
    ofstream rx[40];
    WVector rxcnt(40,0);
    WVector rxflag(40,0);
    UVector<DVector> prefixed_tap_power(2);
    prefixed_tap_power.vect[0] = DVector(filter_length);
    int cp=0;
  //  0,661312823	0,233779161	0,030208484	0,074699532

    prefixed_tap_power.vect[0].vect[cp++]= 0.661312823; //155634 ;
    prefixed_tap_power.vect[0].vect[cp++]= 0.233779161; //466445 ;
    prefixed_tap_power.vect[0].vect[cp++]= 0.030208484; //754797 ;
    prefixed_tap_power.vect[0].vect[cp++]= 0.074699532; //424724 ;
//    prefixed_tap_power.vect[0].vect[cp++]= 0.0899921 ;
    
    prefixed_tap_power.vect[1] = DVector(filter_length);
    cp=0;
//0,616618893	0,001716341	0,006416902	0,375247864
    prefixed_tap_power.vect[1].vect[cp++]= 0.616618893;//126736 ;
    prefixed_tap_power.vect[1].vect[cp++]=0.001716341; //655173 ;
    prefixed_tap_power.vect[1].vect[cp++]=0.006416902; // 137991 ;
    prefixed_tap_power.vect[1].vect[cp++]=0.375247864; // 298397 ;
//    prefixed_tap_power.vect[1].vect[cp++]=0.668285 ;

    cout << "creating channel .. " << endl;
    for(int p =0;p<2;p++) {
        DVector tpower(filter_length,1);
//        dbRandUni(&pseed, tpower);
        tpower.normalizepower();
//        DVector tpower = prefixed_tap_power.vect[p];
        ZMatrix xpath = fadingtaps(max_doppler, symbol_time, filter_length, symtx.taille/nb_porteuse, 1024, &tseed, tpower);
//        ZMatrix xpath = fadingtaps(max_doppler, symbol_time, filter_length, symtx.taille/(nb_porteuse*4) /*it should be enough*/, 1024, &tseed, tpower);
        sprintf(varname, "%d-path_P%d-fr%d-l%d.dat", timestamp, p, nb_ofdm, filter_length);
        ofstream tapsofile(varname);
        cout << "Path " << p  << " power : " << tpower <<endl;
        cout << "saving channel in file .. " <<  varname << endl;
        
    	for (int t=0;t<xpath.line;t++) {
              tapsofile << filter_length << " : [ ";
              for (int l=0;l<filter_length;l++) {
                  tapsofile << xpath.mat[t][l] << " ";
              }
              tapsofile << "]" << endl;
         }
         tapsofile.close();    
         pathfile[p].open(varname);
    }    

    ZMatrix Q=ZMatrix::Q_DFT(nb_porteuse, filter_length);
    ZMatrix h(4,4), x(4,1);
    int symptr=0;
    ZMatrix noise(4,1);
    ZMatrix EYE4 = ZMatrix::eye(4);
    cout << "starting simulation .. "<<endl;
    for(int t=0;symptr<symtx.taille;t++) {
        if(t%1000==0) cout << "T"<<t; cout.flush();
        ZVector tmp1,tmp2;
        pathfile[0] >> tmp1;
/*        h.mat[3][0] = h.mat[0][3] = h.mat[0][0] = tmp.vect[0];
        h.mat[3][3] = -h.mat[3][0];
        h.mat[1][1] = h.mat[1][2] = h.mat[2][2] = -tmp.vect[0].conj();
        h.mat[2][1] = tmp.vect[0].conj();*/
        pathfile[1] >> tmp2;        
/*        h.mat[0][1] = h.mat[0][2] = h.mat[3][2] = tmp.vect[0];
        h.mat[3][1] = -h.mat[3][2] ;
        h.mat[1][0] = h.mat[1][3] = h.mat[2][0] = tmp.vect[0].conj();
        h.mat[2][3] = -tmp.vect[0].conj(); */
        ZVector hf1=Q*tmp1, hf2=Q*tmp2;
        int hptr=0;
        for(int ll=0;hptr<hf1.taille;ll++) {
           h.mat[0][3] = h.mat[0][0] = hf1.vect[hptr];
           h.mat[1][1] = h.mat[1][2] = -hf1.vect[hptr].conj();
           h.mat[0][1] = h.mat[0][2] = hf2.vect[hptr];
           h.mat[1][0] = h.mat[1][3] = hf2.vect[hptr++].conj();

           h.mat[3][3] = -hf1.vect[hptr];
           h.mat[3][0] = hf1.vect[hptr];
           h.mat[2][2] = -hf1.vect[hptr].conj();
           h.mat[2][1] = hf1.vect[hptr].conj();
           
           h.mat[3][1] = -hf2.vect[hptr];
           h.mat[3][2] = hf2.vect[hptr];
           h.mat[2][3] = -hf2.vect[hptr].conj();
           h.mat[2][0] = hf2.vect[hptr++].conj();
           
//        for(int k=0;k<(nb_ofdm>>1);k++) {
//            cout << "K"<<k; cout.flush();
            x.mat[0][0] = symtx.vect[symptr++];
            x.mat[1][0] = symtx.vect[symptr++].conj();
            x.mat[2][0] = symtx.vect[symptr++].conj();
            x.mat[3][0] = symtx.vect[symptr++];
            ZMatrix hx = h*x;
            
            for (snr=snr_min, snr_index=0;snr<=snr_max;snr+=snr_step, snr_index++) {
//                cout << "S"<<snr_index; cout.flush();
                double sig2=convertsnrtosigma2(snr);
                zRandGausInit(&nptr, 0, sig2, rand()%12345876);
                zmRandGaus(&nptr, noise);
                ZMatrix r = hx+noise;
                ZMatrix xest = (h.hermit()*h+sig2*EYE4).inv_nip()*h.hermit()*r;
                xest.mat[1][0] = xest.mat[1][0].conj();
                xest.mat[2][0] = xest.mat[2][0].conj();
//                cout<< "X " << x << endl<< "Xest " << xest <<endl;
                if(rxflag.vect[snr_index]==0) {
                    sprintf(varname, "%d_snr%d-fl%d_rx.dat", timestamp, snr_index, filter_length);
//                    cout << varname << endl;
                    rx[snr_index].open(varname);
                    rxflag.vect[snr_index]=1;
                }
                rx[snr_index] << xest << endl;
                rxcnt.vect[snr_index]++;
            }    
         }        
    }    
    ZMatrix rxtmp;
    cout << endl << "Interpretingresults .. "<<endl;
    sprintf(varname, "%d_result_fl%d_p%d_fr%d.log", timestamp, filter_length, nb_porteuse, nb_ofdm);
    ofstream result(varname);

    for (snr=snr_min, snr_index=0;snr<=snr_max;snr+=snr_step, snr_index++) {
        rx[snr_index].close();
        ZVector symrx(symtx.taille);
        int symrx_cnt=0;
        sprintf(varname, "%d_snr%d-fl%d_rx.dat", timestamp, snr_index, filter_length);
//        cout << varname << endl;
        ifstream rxo(varname);
        for(int c=0;c<rxcnt.vect[snr_index];c++) {
            rxo >> rxtmp;
            for (int l=0;l<rxtmp.line;l++) {
                symrx.vect[symrx_cnt++] = rxtmp.mat[l][0];
            }    
        }    
        rxo.close();
        WVector bitout = (*constel.demodulate_function)(symrx);
        int diff=bitin.diff_count(bitout);
        cout << "Diff for snr index " << snr_index << " is " <<diff << " out of " << bitin.taille << " bits" <<endl;
        result << snr_index << ";" << snr  << ";" << diff <<";" << bitin.taille<<endl;
    }    
    result.close();
    getchar();
    char toto;
    cin >> toto;
    return 0;
}    

