{
	double sigma1 = 9.007;
	double sigma2 = 2.808;
	double sigma = sigma1+sigma2;
	cout<<"--> Sigma = "<<sigma<<" [ub]"<<endl;
	const int nn = 4;
	double f1[] = {1.0e6,1.0e6,1.0e6,1.0e6};
	double f2[] = {4.617e6,4.622e6,4.615e6,4.618e6};
	double f[nn];
	double errf[nn];
	double favg = 0;
	double errfavg = 0;

	for(int i = 0; i < nn; i++)
	{
		f[i] 	= f1[i]/f2[i];
		errf[i] = TMath::Sqrt( TMath::Power(TMath::Sqrt(f1[i])/f2[i],2) + TMath::Power(TMath::Sqrt(f2[i])*f1[i]/(f2[i]*f2[i]),2) );
		favg 	+= f[i];
		errfavg += TMath::Power(errf[i],2);
	}
	favg	= favg/nn;
	errfavg = TMath::Sqrt(errfavg)/nn;
	cout<<"--> Fragmentation function = "<<favg<<" +/- "<<errfavg<<endl;
	cout<<"--> Total cross section = "<<sigma*favg<<" +/- "<<sigma*errfavg<<" [ub]"<<endl;
}
