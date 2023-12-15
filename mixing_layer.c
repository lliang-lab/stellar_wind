/* to solve the four differential equations in Raga, Cabrit 1995*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define x0 1e-5
#define x1 1000
#define pi 3.1415926

#define para_x 0.04 //C-parameter 1/Lambda
#define cos_0  	0.8 //r_0 = r_d * cos_0 * cos_0 
#define epsilon 0
#define alpha	0
#define alpha1  0.03 // entrainment coeff
#define vwcs 53 // v_wind in unit of c_s
#define vdcs 1.78  // v_d in unit of c_s

#define cs1 1   //sound speed of the shocked wind in unit of cs
#define cs2 1  //sound speed of the shocked envelope in unit of cs
#define csl 2 //sound speed of the mixed material

int Nh=5e5; // number of integration step
double *rad, *omega, *sigv, *sigv2, *ram, *v1, *v2, *vram2,*vgam, *vr, *vth;
double dx,term1,term2,term3;
double *v1_ref, *v2_ref, *vl, *vl2, *h1, *h2, *rho1,*rho2,*rhol;
double omega_0;

double a,b,c,d,delta;
int flag,nmax;

//observing angle
double theta_obs[4];

int main(void)
{
	double rk4(double(*f1)(double), double(*f2)(double, double, double, double, double),double(*f3)(double, double, double, double), double(*f4)(double, double, double, double, double), double dx, double x, double y1, double y2, double y3, double y4, int num);
	double f1(double k1);
	double f2(double l, double l1, double l2, double l3, double l4);
	double f3(double m, double m1, double m2, double m3);
	double f4(double n, double n1, double n2, double n3, double n4);
    	double matrics_y1(double p1, double p2, double p3, double p4, double p5, double p6);
    	double matrics_y2(double q1, double q2, double q3, double q4, double q5, double q6);
    	double cal_theta(double radius1, double theta1);
	double cal_ini(double radius2, double theta2);
	int i,flag,index_min; 
    	double xi,vmin;
	double r_0, A, B, C, gamma,T,y1,y2,y3,beta;
    	double a1,a2,b1,b2,c1,c2;
	double *fterm1, *fterm2, *pl;
    	double *vtemp1,*vtemp2,*vtemp3;
    	double inten[100][4],intlow[100][4];
    	double mass,cos1,sin1,cos2,sin2;
    	FILE *velfile,*emfile, tempfile;
	//FILE *tmpfile;
    
    velfile=fopen("../data/ulrich_vel_c004w53d18c1c2c1y.out","w");
    emfile=fopen("../data/ulrich_em_c004w53d18c1c2c1y.out","w");
    //dougfile=fopen("cavity_c004d5","w")
    dx =(x1 - x0)/(1.*Nh);

	//tmpfile=fopen("entrain_bend.tmp","w");

	rad = (double *)malloc(sizeof(double)*Nh);
	omega = (double *)malloc(sizeof(double)*Nh);
	sigv = (double *)malloc(sizeof(double)*Nh);
	sigv2 = (double *)malloc(sizeof(double)*Nh);
	fterm1 = (double *)malloc(sizeof(double)*Nh);
    fterm2 = (double *)malloc(sizeof(double)*Nh);
	vr = (double *)malloc(sizeof(double)*Nh);
	vth = (double *)malloc(sizeof(double)*Nh);
	ram = (double *)malloc(sizeof(double)*Nh);
    v1 = (double *)malloc(sizeof(double)*Nh);
    v2 = (double *)malloc(sizeof(double)*Nh);
    vram2 = (double *)malloc(sizeof(double)*Nh);
    
	//a = para_x * (1.+epsilon*epsilon)*(1.+vdcs*vdcs*(1.-cos_0*cos_0)/cos_0/cos_0);
	//b = -2.*para_x*(1.+epsilon*epsilon)*pow(1.-cos_0*cos_0,0.5)/cos_0/cos_0*vdcs*vdcs;
	//c = para_x*(1.+epsilon*epsilon)*(1+vdcs*vdcs/cos_0/cos_0)-(1.-cos_0*cos_0)/pow(cos_0,3);
  
	omega_0 = 1./pow(1-cos_0*cos_0,0.5);

	fprintf(stderr,"cos_0=%g, omega_0=%g\n",cos_0,omega_0);

	flag=0;
	if(omega_0 > 1./pow(1-cos_0*cos_0,0.5)){flag=1; omega_0 = pow((1-cos_0*cos_0)/pow(cos_0,3)/para_x/(1+epsilon*epsilon)-1.,0.5);}
	fprintf(stderr,"omega_0=%g\n",omega_0);

	r_0 = cos_0*cos_0; //in unit of r_d

	rad[0] = r_0;
	omega[0] = omega_0;
	sigv[0]= 0;
	sigv2[0] = 0; //0.5*dMwinddt/2./pi/r_0/alpha;

	rad[1] = r_0 + omega_0 * dx;
	omega[1] =  omega_0 + dx * (pow(1.+omega_0*omega_0,2)/omega_0*(0.5*para_x*cos_0/(1.-cos_0*cos_0)*(1.+epsilon*epsilon)*(omega_0+1./pow(1.-cos_0*cos_0,1.5))-4.*omega_0/(1.+omega_0*omega_0)/cos_0/cos_0));
	sigv[1] = sigv[0] + (1./pow(cos_0,4) + para_x*vwcs/cos_0/(1.-cos_0*cos_0)*(1.+epsilon*epsilon)*alpha/vwcs*pow(1+omega_0*omega_0,0.5)) * dx;
	sigv2[1] = sigv2[0]+omega_0/pow(1+omega_0*omega_0,0.5)/pow(cos_0,4) * dx;	

//use RK method to find the analytic solution for the cavities
	for(i=2;i<Nh;i++){
        xi=x0+dx*i;
		rk4(f1,f2,f3,f4,dx, xi,rad[i-1],omega[i-1],sigv[i-1],sigv2[i-1],i);
		fterm1[i]=term1;
        //thermal pressure: rho_d c^2_s  in unit of [M_w v_w /4pi/r^2_d]
		fterm2[i]=term2;
        //wind ram pressure in unit of [M_w v_w /4pi/r^2_d]
	}
    
    fprintf(stderr,"RK method done.\n");
	
//+++++++++++++++++++++++++++++++++++++++++++++++++++++//
//	NOW WE NEED TO SOLVE THE STATE OF THE MIXING LAYER //	
//+++++++++++++++++++++++++++++++++++++++++++++++++++++//
	vl = (double *)malloc(sizeof(double)*Nh);
	vl2 = (double *)malloc(sizeof(double)*Nh);
	v1_ref = (double *)malloc(sizeof(double)*Nh);//in a different reference frame
	v2_ref = (double *)malloc(sizeof(double)*Nh);//in a different reference frame
	h1 = (double *)malloc(sizeof(double)*Nh);
 	h2 = (double *)malloc(sizeof(double)*Nh);
	rho1 = (double *)malloc(sizeof(double)*Nh);//gas density
	rho2 = (double *)malloc(sizeof(double)*Nh);
	rhol = (double *)malloc(sizeof(double)*Nh);
    vtemp1 = (double *)malloc(sizeof(double)*Nh);//
    vtemp2 = (double *)malloc(sizeof(double)*Nh);//
    vtemp3 = (double *)malloc(sizeof(double)*Nh);//
		
	for(i=0;i<Nh;i++){
		//xi = x0 + dx * i;
		xi=x0+dx*i;
        rho1[i] = pow(rad[i]-xi*omega[i],2)/(1.+omega[i]*omega[i])/pow(xi*xi+rad[i]*rad[i],2);//pressure
		 rho1[i] = rho1[i]/cs1/cs1;//\frac{\dot{M}_{w} v_w}{4\pi r^2_d}
		rho1[i] = rho1[i]*vwcs;//\frac{\dot{M}_{w}}{4\pi r^2_d * c_s}
		
		rho2[i] = rho1[i]*cs1*cs1/cs2/cs2;//density of the shocked inflow, in the same unit
		rhol[i] = rho1[i]*cs1*cs1/csl/csl;//the density of the gas in the mixing layer
	
        v1[i]=sigv2[i]/(1e-20+sigv[i]);//velocity of the shocked stellar wind
        //velocity of the shocked inflowing envelope
        //y=(a1(:,2)-a1(:,3).*a1(:,1)).^2./(1+a1(:,3).*a1(:,3))./(a1(:,1).^2+a1(:,2).^2).^2;
		//x=a1(:,1);x(1)=1e-5;
		//plot(x,y,'k-','linewidth',2);
        //velocity in unit of v_w
		h1[i]=0.;
		h2[i]=0.;
	}
    
    //mean velocity
    for(i=0;i<Nh;i++){
        v2[i]=0.;
        vtemp1[i]=0.;
        vtemp2[i]=0.;
        vtemp3[i]=0.;
    }
    for(i=0;i<Nh;i++){
        //xi=x0+dx*i;
        xi=x0+dx*i;
        
        vtemp1[i]=pow((omega[i]*rad[i]+xi)/pow(1.+omega[i]*omega[i],0.5)/(xi*xi+rad[i]*rad[i]),0.5);//v_m from the wind side
        //vtemp2[i]=fterm1[i]*ram[i]*(omega[i]*rad[i]+xi)/pow((xi*xi+rad[i]*rad[i]),0.5)*rad[i];//I MISTAKENLY USED RAM INSTEAD OF V2
        vtemp2[i]=fterm1[i]*ram[i]*(omega[i]*rad[i]+xi)/(rad[i]-omega[i]*xi)*pow(1+omega[i]+omega[i],0.5)*rad[i];
        vtemp2[i]*=vwcs;//pressure in unit of [mw*cs/4/pi/r^2_d]
                        //or density in unit of [mw/4/pi/rd^2/cs]
    }
    //y2(length(y2))=y2(length(y2))*(a1(length(y2),1)-a1(length(y2)-1,1));
    vtemp2[Nh-1]=vtemp2[Nh-1]*dx;
        //for n=length(y2)-1:-1:1
    for(i=Nh-1;i>0;i--){
        vtemp2[i-1]=vtemp2[i-1]*dx;
        vtemp2[i-1]=vtemp2[i-1]+vtemp2[i];
    }
    //y2=R*sigma_m*v^2_m [mw*cs/4/pi]
    for(i=0;i<Nh;i++){
        vtemp3[i]=fterm1[i]*pow(ram[i],0.5)*pow((1+omega[i]*omega[i]),0.5)*rad[i];
        vtemp3[i]*=vwcs;
    }//fterm1 pressure [mw*cs/4/pi/r^2_d]
     //vtemp3 in unit of [mw/4/pi/r_d]
    vtemp3[Nh-1]=vtemp3[Nh-1]*dx;
    //for n=length(y3)-1:-1:1
    for(i=Nh-1;i>0;i--){
        vtemp3[i-1]=vtemp3[i-1]*dx;
        vtemp3[i-1]=vtemp3[i-1]+vtemp3[i];
    }
    //y3=2*\pi*R*sigma_m*v_m [mw/4/pi]
    vmin=0.;
    for(i=0;i<Nh;i++){
        v2[i]=-vtemp2[i]/(vtemp3[i]+1e-20); //y2/y3=v_m (in unit of c_s)
        //if(i<100) fprintf(stderr,"v2=%g\n",v2[i]); //in unit of c_s
        if(vmin>v2[i] && i<1e3) {vmin=v2[i];index_min=i;}
    }//velocity of the shocked inflowing envelope
  
    //for(i=0;i<100;i++) fprintf(stderr,"v2=%g\n",v2[i]);
    //to define the average velocities
	for(i=1;i<Nh;i++){
		v1_ref[i]=v1[i]*vwcs-vmin;//in unit of cs
		v2_ref[i]=-vmin+v2[i];//make sure that in this ref frame, velocities at both sides are positive
        //v2_ref[i]=0;
		vl[i]=0.5*(v1_ref[i]+v2_ref[i]);
		vl2[i]=(pow(v1_ref[i],3)-pow(v2_ref[i],3))/3./(v1_ref[i]-v2_ref[i]);
	}
    
//To make the emission line intensity plot
//    for(i=0;i<100;i++){
//        for(j=0;j<4;j++) {inten[i][j]=0.;intlow[i][j]=0.;}
//    }
    
 //   mass = 0.;
    
//the viewing angles from the disc plane
//    theta_obs[0] = 0.;
//    theta_obs[1] = 0.524;
//    theta_obs[2] = 1.047;
//    theta_obs[3] = 1.571;

/////////////////////////////////////////////////////////////////////
/////      start calculating the central mixing layer           /////
/////////////////////////////////////////////////////////////////////
    
    for(i=2;i<Nh;i++){
        if(rad[i]<0) break;
        xi=x0+dx*i;
        a1=rho1[i-1]*v1_ref[i-1]-rhol[i-1]*vl[i-1];
        //rho density [\frac{Mw}{4\pi r^2_d * c_s}]
        b1=rho2[i-1]*v2_ref[i-1]-rhol[i-1]*vl[i-1];
        //b1=-rhol[i-1]*vl[i-1];
        a2=rho1[i-1]*v1_ref[i-1]*v1_ref[i-1]-rhol[i-1]*vl2[i-1];
        b2=rho2[i-1]*v2_ref[i-1]*v2_ref[i-1]-rhol[i-1]*vl2[i-1];
        c1=(h1[i-1]+h2[i-1])*vl[i-1]*(rhol[i]-rhol[i-1])/(dx)+rhol[i-1]*(h1[i-1]+h2[i-1])*(vl[i]*rad[i]-vl[i-1]*rad[i-1])/(dx)/rad[i-1]-alpha1*rhol[i-1]*csl*pow(1+omega[i-1]*omega[i-1],0.5);
        c2=(h1[i-1]+h2[i-1])*vl2[i-1]*(rhol[i]-rhol[i-1])/(dx)+rhol[i-1]*(h1[i-1]+h2[i-1])*(vl2[i]*rad[i]-vl2[i-1]*rad[i-1])/(dx)/rad[i-1]-alpha1*rhol[i-1]*csl*v2_ref[i-1]*pow(1+omega[i-1]*omega[i-1],0.5)+(h1[i-1]+h2[i-1])*(rhol[i]*csl*csl-rhol[i-1]*csl*csl)/(dx);
        h1[i]=h1[i-1]+dx*matrics_y1(a1,a2,b1,b2,c1,c2);
        h2[i]=h2[i-1]+dx*matrics_y2(a1,a2,b1,b2,c1,c2);
        if(i<1000||i%20==0) fprintf(velfile,"%g %g %g %g %g %g %g %g %g %g %g %g %g\n",xi, v1_ref[i], v2_ref[i],a1,a2,b1,b2,c1,c2,h1[i],h2[i],h1[i]+h2[i],rhol[i]);
    }

	for(i=2;i<Nh;i++){
		if(rad[i]<0) break;
		xi=x0+dx*i;
		term3=para_x/pow(xi*xi+r_0*r_0,0.75)/pow(1+xi/pow(xi*xi+r_0*r_0,0.5)/cal_theta(pow(xi*xi+r_0*r_0,0.5),atan(r_0/(xi+1e-20))),0.5)/(xi/pow(xi*xi+r_0*r_0,0.5)/2./cal_theta(pow(xi*xi+r_0*r_0,0.5),atan(r_0/(1e-20+xi)))+1./pow(xi*xi+r_0*r_0,0.5)*pow(cal_theta(pow(xi*xi+r_0*r_0,0.5),atan(r_0/(1e-20+xi))),2))*(1.+epsilon*epsilon);
        //if(i<100) {
          //  fprintf(stderr,"%g %g %g\n",v1_ref[i],v2_ref[i],vl[i]);
            //fprintf(stderr,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", x0 + dx * i, rad[i], omega[i], sigv[i], sigv2[i]/sigv[i], -sigv2[i]/sigv[i]*(omega[i]/pow(1+omega[i]*omega[i],0.5)),fterm1[i],fterm2[i],term3, ram[i], v2[i], vr[i],vth[i],h1[i],h2[i],rhol[i]);
        //}
		if(i<1000||i%20==0) fprintf(stdout,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", xi, rad[i], omega[i], sigv[i], sigv2[i]/sigv[i], -sigv2[i]/sigv[i]*(omega[i]/pow(1+omega[i]*omega[i],0.5)),fterm1[i],fterm2[i],term3, ram[i], vram2[i], vr[i],vth[i], v2[i], h1[i],h2[i],rhol[i]);//vr,vth in unit of cs
	}
    
//////////////////////////////////////////////////////////
/////////           Emission Properties       ////////////
//////////////////////////////////////////////////////////

//    for(i=1;i<Nh;i++){
  //      if(rad[i]<0) break;
    //    if(i%1000==0) fprintf(stderr,"i=%d\n",i);
        
      //  for(j=0;j<100;j++){
        //    phi = 2*pi*j/100.;
          //  cos1 = omega[i]/pow(1.+omega[i]*omega[i],0.5);
            //sin1 = 1./pow(1.+omega[i]*omega[i],0.5);
            
            //for(k=0;k<4;k++){
              //  for(m=0;m<20;m++){
                //    bin = (-0.05*m*sigv2[i]/(1.e-20+sigv[i])*cos(theta_obs[k])*cos1*cos(phi)-0.05*m*sigv2[i]/(1.e-20+sigv[i])*sin(theta_obs[k])*sin1+1.)/(2./100.);
                  //  if(bin<0) bin = 0;
                    //if(bin>99) bin = 99;
                    //inten[bin][k] += 2*pi*rad[i]*0.01*(dx*pow(1+omega[i]*omega[i],0.5))*sigv[i]*sigv[i]/(1.e-20+sigv2[i])*pow(rad[i]-omega[i]*(x0+dx*i),2)/(1.+omega[i]*omega[i])/pow(pow(x0+dx*i,2)+rad[i]*rad[i],2)/20.;
                    //the lower half of the cavity
                    
                    //bin = (-0.05*m*sigv2[i]/(1.e-20+sigv[i])*cos(theta_obs[k])*cos1*cos(phi)+0.05*m*sigv2[i]/(1.e-20+sigv[i])*sin(theta_obs[k])*sin1+1.)/(2./100.);
                    //if(bin<0) bin = 0;
                    //if(bin>99) bin = 99;
                    //intlow[bin][k] += 2*pi*rad[i]*0.01*(dx*pow(1+omega[i]*omega[i],0.5))*sigv[i]*sigv[i]/(1.e-20+sigv2[i])*pow(rad[i]-omega[i]*(x0+dx*i),2)/(1.+omega[i]*omega[i])/pow(pow(x0+dx*i,2)+rad[i]*rad[i],2)/20.;
                //}
            //}
        //}
        //mass += 2.*2.*pi*rad[i]*(dx*pow(1.+omega[i]*omega[i],0.5))*pow(sigv[i]*sigv[i]/(sigv2[i]+1e-20),1)*pow(rad[i]-omega[i]*(x0+dx*i),2)/(1.+omega[i]*omega[i])/pow(pow(x0+dx*i,2)+rad[i]*rad[i],2);
    //}
    
    //fprintf(stderr,"start normalization.\n");
    
    //for(i=0;i<100;i++){
      //  for(j=0;j<4;j++) {inten[i][j] = inten[i][j]/mass;intlow[i][j] = intlow[i][j]/mass;}
    //}
    
    //for(i=0;i<100;i++){
      //  fprintf(emfile,"%g %g %g %g %g %g %g %g\n",inten[i][0],inten[i][1],inten[i][2],inten[i][3],intlow[i][0],intlow[i][1],intlow[i][2],intlow[i][3]);
    //}

    fprintf(stderr,"%g\n", matrics_y1(1,1,1,-1,5,4));
    fprintf(stderr,"%g\n", matrics_y2(1,1,1,-1,5,4));
    fprintf(stderr,"vmin=%g (cs) index=%d\n", vmin,index_min);
	return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++//
//	THE SOLUTION FOR 2X2 MATRICS //	
//+++++++++++++++++++++++++++++++++++++++++++++++++++++//
double matrics_y1(double a1, double a2, double b1, double b2, double c1, double c2){
	double determin;//the determinant of the 2x2 matric
	
	determin = a1*b2-a2*b1;
	if(determin==0) fprintf(stderr,"determinant = 0!!!\n");
	return (c1*b2-c2*b1)/determin;
}
double matrics_y2(double a1, double a2, double b1, double b2, double c1, double c2){
	double determin;
	
	determin = a1*b2-a2*b1;
	if(determin==0) fprintf(stderr,"determinant = 0!!!\n");
	return (c2*a1-c1*a2)/determin;
}

//Runge-Kutta solver for the 4 OD equations
double rk4(double(*f1)(double), double(*f2)(double, double, double, double, double),double(*f3)(double, double, double, double), double(*f4)(double, double, double, double, double), double dx, double x, double y1, double y2, double y3, double y4, int num)
{
	double	k1 = dx * f1(y2),
		l1 = dx * f2(x, y1, y2, y3, y4),
		m1 = dx * f3(x, y1, y2, y3),
		n1 = dx * f4(x, y1, y2, y3, y4); 

	double	k2 = dx * f1(y2 + l1 / 2),
		l2 = dx * f2(x + dx / 2, y1 + k1 / 2, y2 + l1 / 2, y3 + m1 /2, y4 + n1 / 2),
		m2 = dx * f3(x + dx / 2, y1 + k1 / 2, y2 + l1 / 2, y3 + m1 /2),			
		n2 = dx * f4(x + dx / 2, y1 + k1 / 2, y2 + l1 / 2, y3 + m1 /2, y4 + n1 / 2);

	double	k3 = dx * f1(y2 + l2 / 2),
		l3 = dx * f2(x + dx / 2, y1 + k2 / 2, y2 + l2 / 2, y3 + m2 /2, y4 + n2 / 2),
		m3 = dx * f3(x + dx / 2, y1 + k2 / 2, y2 + l2 / 2, y3 + m2 /2),	
		n3 = dx * f4(x + dx / 2, y1 + k2 / 2, y2 + l2 / 2, y3 + m2 /2, y4 + n2 / 2);

	double	k4 = dx * f1(y2 + l3),
		l4 = dx * f2(x + dx, y1 + k3, y2 + l3, y3 + m3, y4 + n3),
	
	m4 = dx * f3(x + dx, y1 + k3, y2 + l3, y3 + m3),		
		n4 = dx * f4(x + dx, y1 + k3, y2 + l3, y3 + m3, y4 + n3);

	int n = (x - x0) / dx; 

	rad[num] = y1 + 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4); 
	omega[num] = y2 + 1. / 6. * (l1 + 2 * l2 + 2 * l3 + l4); 
	sigv[num] = y3 + 1. / 6. * (m1 + 2 * m2 + 2 * m3 + m4); 
	sigv2[num] = y4 + 1. / 6. * (n1 + 2 * n2 + 2 * n3 + n4); 
 	
	if(num<10){
        	fprintf(stderr,"k:%g %g %g %g\n",k1,k2,k3,k4);	
		fprintf(stderr,"l:%g %g %g %g\n",l1,l2,l3,l4);
	        fprintf(stderr,"m:%g %g %g %g\n",m1,m2,m3,m4);
	        fprintf(stderr,"n:%g %g %g %g\n",n1,n2,n3,n4);
		fprintf(stderr,"%g %g\n",1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4),dx*f1(y2));
		fprintf(stderr,"%g %g\n",1. / 6. * (l1 + 2 * l2 + 2 * l3 + l4),dx*f2(x,y1,y2,y3,y4));
		fprintf(stderr,"%g %g\n",1. / 6. * (m1 + 2 * m2 + 2 * m3 + m4),dx*f3(x,y1,y2,y3));
		fprintf(stderr,"%g %g\n",1. / 6. * (n1 + 2 * n2 + 2 * n3 + n4),dx*f4(x,y1,y2,y3,y4));

	}
}

double f1(double x)
{
	return x;
}

double f2(double x, double y1, double y2, double y3, double y4)
{
	//FILE *tmpfile;
	//tmpfile = fopen("entrain_ulrich.tmp","a");
	int num;
	double cal_theta(double radius1, double theta1);
	double vram,pram;

	num = (x - x0)/dx;

	if(x<100){
		//fprintf(stderr,"atan=%g cos=%g dw/dz=%g\n",atan(y1/(x+1e-20)),cal_theta(pow(x*x+y1*y1,0.5)/r_d,atan(y1/(x+1e-20))),pow(1+y2*y2,1.5)/y3/y4*(-dMinfdt/8./pi/pow(G*Mstar,0.5)/pow(x*x+y1*y1,0.75)/pow(1+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5)/r_d,atan(y1/(x+1e-20))),0.5)/(x/pow(x*x+y1*y1,0.5)/2./cal_theta(pow(x*x+y1*y1,0.5)/r_d,atan(y1/(1e-20+x)))+r_d/pow(x*x+y1*y1,0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5)/r_d,atan(y1/(1e-20+x))),2))*cs*cs+(dMwinddt/4./pi*v_wind)*pow(y1-y2*x,2)/pow(x*x+y1*y1,2)/(1+y2*y2)));
	}
    term1 = para_x/pow(x*x+y1*y1,0.75)/pow(1+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(x+1e-20))),0.5)/(x/pow(x*x+y1*y1,0.5)/2./cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))+1./pow(x*x+y1*y1,0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),2))*(1.+epsilon*epsilon);
	term2 = pow(y1-y2*x,2)/pow(x*x+y1*y1,2)/(1+y2*y2);
 
	//term3 = para_x/pow(x*x+y1*y1,0.75)/pow(1+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(x+1e-20))),0.5)/(x/pow(x*x+y1*y1,0.5)/2./cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))+1./pow(x*x+y1*y1,0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),2))*(1.+epsilon*epsilon);

        //fprintf(stderr,"%g %g\n",term1,term2);

	vr[num] = -pow(vdcs*vdcs/pow(x*x+y1*y1,0.5),0.5)*pow(1.+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),0.5);//in unit of cs

	vth[num] = -pow(vdcs*vdcs/pow(x*x+y1*y1,0.5),0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))) + x/pow(x*x+y1*y1,0.5),0.5)/pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),0.5)/(y1/pow(x*x+y1*y1,0.5))*(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))-x/pow(x*x+y1*y1,0.5));//in unit of cs

	vram = vr[num]*(x*y2-y1)/pow(x*x+y1*y1,0.5)/pow(1+y2*y2,0.5)+vth[num]*(y1*y2+x)/pow(x*x+y1*y1,0.5)/pow(1+y2*y2,0.5);//in unit of cs

	pram = vram *vram;//in unit of cs^2

	//fprintf(stderr,"z=%g,R=%g vr=%g vth=%g vram = %g pram=%g\n",x,y1,vr[num],vth[num],vram,pram);

	ram[num] = pram;
	vram2[num] = vr[num]*vr[num]+vth[num]*vth[num];

	if(vram > 0) return pow(1+y2*y2,1.5)/y4*(-para_x/pow(x*x+y1*y1,0.75)/pow(1+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(x+1e-20))),0.5)/(x/pow(x*x+y1*y1,0.5)/2./cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))+1./pow(x*x+y1*y1,0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),2))*(1.+epsilon*epsilon)*(1.+pram)+pow(y1-y2*x,2)/pow(x*x+y1*y1,2)/(1+y2*y2));
	else {ram[num]=0; return pow(1+y2*y2,1.5)/y4*(-para_x/pow(x*x+y1*y1,0.75)/pow(1+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(x+1e-20))),0.5)/(x/pow(x*x+y1*y1,0.5)/2./cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))+1./pow(x*x+y1*y1,0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),2))*(1.+epsilon*epsilon)+pow(y1-y2*x,2)/pow(x*x+y1*y1,2)/(1+y2*y2));}
//	return 
//	ratio = para_x/pow(x*x+y1*y1,0.75)/pow(1+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(x+1e-20))),0.5)/(x/pow(x*x+y1*y1,0.5)/2./cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))+1./pow(x*x+y1*y1,0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),2))*(1.+epsilon*epsilon)/pow(y1-y2*x,2)/pow(x*x+y1*y1,2)/(1+y2*y2);
//	fprintf(stderr,"%g\n",ratio);
}

double f3(double x, double y1, double y2, double y3)
{
	double cal_theta(double radius1, double theta2);
	return (y1-x*y2)/pow(x*x+y1*y1,1.5)-y3*y2/y1+para_x*vwcs/pow(x*x+y1*y1,0.75)/pow(1+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(x+1e-20))),0.5)/(x/pow(x*x+y1*y1,0.5)/2./cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))+1./pow(x*x+y1*y1,0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),2))*(1.+epsilon*epsilon)*alpha/vwcs*pow(1+y2*y2,0.5);

}

double f4(double x, double y1, double y2, double y3, double y4)
{
	return -y2/y1*y4+(y1-y2*x)*(y2*y1+x)/pow(1+y2*y2,0.5)/pow(x*x+y1*y1,2);

//	return 1./y3/pow(x*x+y1*y1,2)*((y1-x*y2)*(y2*y1+x)/pow(1+y2*y2,0.5)*dMwinddt*v_wind/4./pi-y4*(y1-y2*x)*pow(x*x+y1*y1,0.5)*dMwinddt/4./pi);
}
//r,omega, sigma*v, v as a function of z
//initial values
//	r = 1, omega = ctan(beta) = 0, sigma*v = 0, v=0
// k, l, m, n

double cal_theta(double radius, double theta){
	double y1, y2, y3, T, beta;

	c = radius -1 ;
	d = -radius * cos(theta);

	delta = 81.*radius*radius*cos(theta)*cos(theta)+12*pow(radius-1,3);
	
	if(delta > 0){

		y1 = 3./2.*(-9.*radius*cos(theta)+pow(delta,0.5));
		y2 = 3./2.*(-9.*radius*cos(theta)-pow(delta,0.5));
		if (theta<0.79 && theta >0.78) fprintf(stderr,"delta=%g y1=%g y2=%g radius=%g %g theta_0=%g\n",delta, y1,y2, radius, 81.*radius*radius*cos(theta)*cos(theta),-1./3.*(pow(y1,1/3.)-pow(-y2,1/3.)));

		if(y1<0 && y2<0) return  -1./3.*(-pow(-y1,1./3.)-pow(-y2,1./3.));
		else if(y1<0 && y2>=0) return  -1./3.*(-pow(-y1,1./3.)+pow(y2,1./3.));
		else if(y1>=0 && y2<0) return -1./3.*(pow(y1,1./3.)-pow(-y2,1./3.));
		else return  -1./3.*(pow(y1,1./3.)+pow(y2,1./3.));


	}else if(delta == 0){

		y1 = -3.*radius*cos(theta)/(radius-1);
		y2 = 3./2. *radius*cos(theta)/(radius-1);

		if(radius < 1){
			if(y1<=1 && y1>=0) return y1;
			else fprintf(stderr,"cosine > 1!, radius=%g theta=%g y1=%g y2=%g\n", radius, theta, y1, y2);
		}else if(radius == 1){
			fprintf(stderr,"Wrong value!");
			return -1;
		}else{
			if(y2<=1 && y2>=0) return y2;
			else fprintf(stderr,"cosine > 1!, radius=%g theta=%g y1=%g y2=%g\n", radius, theta, y1, y2);
		}
	}else{
		T = - 27./2.*radius*cos(theta)/pow(-27.*pow(radius-1,3),0.5);
		beta = acos(T);
		y1 = pow(-3*(radius-1),0.5)/3.*(cos(beta/3.)+1.732*sin(beta/3.));
		y2 = pow(-3*(radius-1),0.5)/3.*(cos(beta/3.)-1.732*sin(beta/3.));
		y3 = -2 * pow(-3*(radius-1),0.5)*cos(beta/3.)/3;

		if(y1>=0 && y1<=1 ) return y1;
		//else if(y2>=0 && y2<=1 && y2 >= cos(theta)) return y2;
		//else if(y3>=0 && y3<=1 && y3 >= cos(theta)) return y3;
		//else fprintf(stderr,"no solution!\n");
	}
}

double cal_ini(double radius, double theta){
	double y1, y2, y3, T, beta;

	c = radius -1 ;
	d = -radius * cos(theta);

	if(delta > 0){

		y1 = 3./2.*(-9.*a+pow(delta,0.5));
		y2 = 3./2.*(-9.*a-pow(delta,0.5));
		//if (theta<0.79 && theta >0.78) fprintf(stderr,"delta=%g y1=%g y2=%g radius=%g %g\n",delta, y1,y2, radius);

		if(y1<0 && y2<0) return  -1./3./a*(-pow(-y1,1./3.)-pow(-y2,1./3.));
		else if(y1<0 && y2>=0) return  -1./3./a*(-pow(-y1,1./3.)+pow(y2,1./3.));
		else if(y1>=0 && y2<0) return -1./3./a*(pow(y1,1./3.)-pow(-y2,1./3.));
		else return  -1./3./a*(pow(y1,1./3.)+pow(y2,1./3.));


	}else if(delta == 0){

		y1 = 9.*a;
		y2 = -9./2.*a;

		if(y1<=1 && y1>=0) return y1;
		else fprintf(stderr,"cosine > 1!, radius=%g theta=%g y1=%g y2=%g\n", radius, theta, y1, y2);
	}else{
		T = (2-27*a*a)/2.;
		beta = acos(T);

		y1 = 1./3./a*(-1+cos(beta/3.)+1.732*sin(beta/3.));
		y2 = 1./3./a*(-1+cos(beta/3.)-1.732*sin(beta/3.));
		y3 = (-1-2*cos(beta/3))/3./a;

		fprintf(stderr,"T=%g beta=%g y1=%g\n",T,beta, y1);

		if(y1>=0 && y1<=1 ) return y1;
		//else if(y2>=0 && y2<=1 && y2 >= cos(theta)) return y2;
		//else if(y3>=0 && y3<=1 && y3 >= cos(theta)) return y3;
		else {}//fprintf(stderr,"no solution!\n");
	}
}

