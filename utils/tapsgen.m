%cmplx_mat tapsgen(double  symbol_T, int filter_length, long nb_samples, long nb_rays, const DVector &fd, const DMatrix &g, dRandUniStatePtr *tapseed){
function [fadtaps] = tapsgen(symbol_T, filter_length, nb_samples, nb_rays, fd, g,  tapseed)
% int i,j,k;
%	cmplx sample=DCplx(0,0);
	symbol_T=symbol_T*1e-6;
    %cmplx_mat fadtaps(nb_samples, filter_length);

    fadtaps = zeros(nb_samples, filter_length);


%	 Generate the coefficients of tap-delay-line 

%    double coef=sqrt(2.0/nb_rays);
    coef=sqrt(2.0/nb_rays);
%	DVector theta(nb_rays);
    theta = zeros(nb_rays);

%	for(int t=0;t<theta.taille;t++) {
    for t=1:nb_rays
%		theta.vect[t]=2.0*M_PI*dRandUni(tapseed);
        theta(t)=2*pi*uni_rand(tapseed);
%	}
    end


%	DVector norm(filter_length);
    norm = zeros(filter_length);
%	for(j=0; j<nb_samples; j++){
    for jx=1:nb_samples
%		for(k=0; k<filter_length; k++){
        for k=1:filter_length
%			sample.re=0; sample.im=0;
			sample = 0;
			%for(i=0;i<nb_rays;i++){				
            for ix=1:nb_rays
%				double tau=theta.vect[i]+2*M_PI*fd.vect[i]*j*symbol_T; 
				tau=theta(ix)+2*pi*fd(ix)*jx*symbol_T; 
%				sample.re+=cos(tau)*g.mat[i][k];
%				sample.im+=sin(tau)*g.mat[i][k];
                sample = sample + cos(tau)*g(ix,k) + i*sin(tau)*g(ix,k);
%			}
            end;
%			sample*=coef;
			sample=sample*coef;
%			fadtaps.mat[j][k]=sample;
			fadtaps(jx,k)=sample;
%			norm.vect[k]+=sample.re*sample.re+sample.im*sample.im;
			norm(k)= norm(k)+abs(sample)^2;
%		}
        end
%	}
    end

%// normalization
%	for(register int kk=0;kk<filter_length;kk++) {
    for k=1:filter_length
%		norm.vect[kk]=1.0/sqrt(norm.vect[kk]/nb_samples);
%		for(register int jj=0;jj<nb_samples;jj++) {
%			fadtaps.mat[jj][kk] *= norm.vect[kk];
%		}
        fadtaps(:,k) = fadtaps(:,k)/sqrt(norm(k)/nb_samples);
%	}
    end
%	return fadtaps;
%}

