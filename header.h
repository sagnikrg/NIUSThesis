/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Definition for Function Calls
   for all Analytical Functions
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*Amplitude and Angle*/



double A2 (double xf, double xi, double t)
    {
      double wnull;
      wnull=1.0;


    double D;
    D = xf*xf+xi*xi+2*xf*xi*cosh(wnull*t)/(sinh(wnull*t)*sinh(wnull*t));
    return D;
}


complex double phi (double xf, double xi, double t)
    {
      double wnull;
      wnull=1.0;

    complex double D;
    D = cacos(-I*xi*sinh(wnull*t)/csqrt(xf*xf+xi*xi-2*xf*xi*cosh(wnull*t)));
    return D;
}

/*Angular Frequency*/


complex double omega (double xf, double xi, double t)
    {
      double wnull;
      double K;
      double lambda;

      wnull=1.0;
      lambda=0.1;
      K=1;

    complex double D;

    D = csqrt(-wnull*wnull+3*K*lambda*A2(xf,xi,t)/2+3*K*lambda*A2(xf,xi,t)*(ccos(2*(I*wnull-3*I*lambda*A2(xf,xi,t)/(8*wnull))*t+2*phi(xf,xi,t)))/2);
    return D;
}


/*Classical Path and Action*/

complex double Xcl (double xf, double xi, double t)
    {
    complex double D;
    double wnull;
    double Kpot;
    double lambda;

     wnull=1.0;
     lambda=0.1;
     Kpot=1;

    D = csqrt(A2(xf,xi,t))*ccos(omega(xf,xi,t)*t+phi(xf,xi,t))+lambda*Kpot/(8*wnull*wnull)*csqrt(A2(xf,xi,t))*A2(xf,xi,t)*ccos(omega(xf,xi,t)*t+3*phi(xf,xi,t))*csin((omega(xf,xi,t)*t))*csin((omega(xf,xi,t)*t));
    return D;
}

complex double L (double xf, double xi, double t)
    {
    complex double D;
    double wnull;
    double Kpot;
    double lambda;
    double h;

     wnull=1.0;
     lambda=0.1;
     Kpot=1;
     h=0.05;

    D = 0.125*(Xcl(xf,xi,t+h)-Xcl(xf,xi,t-h))*(Xcl(xf,xi,t+h)-Xcl(xf,xi,t-h))-wnull*wnull*Xcl(xf,xi,t)*Xcl(xf,xi,t)/2-Kpot*lambda*Xcl(xf,xi,t)*Xcl(xf,xi,t)*Xcl(xf,xi,t)*Xcl(xf,xi,t)/24;
    return D;
}



complex double Scl (double xf, double xi, double t)
    {
    complex double D;
    complex double P;
    double tprime, h;

    h=0.1;
    P=0.1;
     for (tprime = 0; tprime < t; tprime=tprime+h) {
        P= P+h*L(xf,xi,tprime);
      }
    D=P;
    P=0.0; 

    return D;
}


/*Initial Wave Packet (Gaussian)*/


double Wpkt (double x)
    {
    double D;
    double sigma;
    double a;

    a=-1.5;
    sigma=1.0;

    D = 1/sqrt(2*M_PI)/sigma*exp(-(x-a)*(x-a)/sigma/sigma);
    return D;
}
