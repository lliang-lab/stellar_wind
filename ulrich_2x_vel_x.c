/* To solve equations 2-5 in Barral and Conto 81 by using the enviromment density characterized by eq 13 of ulrich 71 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define z0 1.5
#define Nh 4e5
#define x0 0
#define x1 5000
#define pi 3.1415926
#define dMinfdt 1e-6 
#define dMwinddt 1e-6
#define Mstar 1
#define G  3.9458e8
#define v_wind 6806.1 // 150km s-1
#define cs 453.74 // 10km s-1
#define r_d 2150.97 // 10AU
#define para_x 0.005
#define cos_0  	0.966
#define epsilon	0
#define alpha	0
#define vwcs 141.4
//#define omega_0 27.9
#define vdcs 7.07

double *rad, *omega, *sigv, *sigv2, *ram, *v2, *vgam, *vr, *vth;
double dx = (x1 - x0)/Nh,term1,term2,term3;
//double omega_0 = pow(((1-cos_0*cos_0)/pow(cos_0,3)/(1.+epsilon*epsilon)/para_x)-1.,0.5);
//double P0 = 0;
double omega_0;

double a,b,c,d,delta;
int flag;
//FILE *tmpfile;

//tmpfile=fopen("entrain_ulrich.tmp","w");

int main(void)
{
	double rk4(double(*f1)(double), double(*f2)(double, double, double, double, double),double(*f3)(double, double, double, double), double(*f4)(double, double, double, double, double), double dx, double x, double y1, double y2, double y3, double y4, int num);
	double f1(double k1);
	double f2(double l, double l1, double l2, double l3, double l4);
	double f3(double m, double m1, double m2, double m3);
	double f4(double n, double n1, double n2, double n3, double n4);
	double cal_theta(double radius1, double theta1);
	double cal_ini(double radius2, double theta2);
	int i,flag; 
	double xi;
	double r_0, A, B, C, gamma,T,y1,y2,y3,beta;
	double *fterm1, *fterm2;
	//FILE *tmpfile;

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
        v2 = (double *)malloc(sizeof(double)*Nh);

	//a = para_x*(1.+epsilon*epsilon)*(1.+vdcs*vdcs*(1-cos_0*cos_0)/cos_0/cos_0);
	//b = -2.*para_x*(1.+epsilon*epsilon)*vdcs*vdcs*pow(1-cos_0*cos_0,0.5)/cos_0/cos_0;
	//c = para_x * (1.+epsilon*epsilon)*(1.+vdcs*vdcs/cos_0/cos_0)-(1.-cos_0*cos_0)/pow(cos_0,3);
	//cos_0 = 0.5*(-para_x*vdcs*vdcs+pow(4+pow(para_x*vdcs*vdcs,2),0.5));
	a = para_x * (1.+epsilon*epsilon)*(1.+vdcs*vdcs*(1.-cos_0*cos_0)/cos_0/cos_0);
	b = -2.*para_x*(1.+epsilon*epsilon)*pow(1.-cos_0*cos_0,0.5)/cos_0/cos_0*vdcs*vdcs;
	c = para_x*(1.+epsilon*epsilon)*(1+vdcs*vdcs/cos_0/cos_0)-(1.-cos_0*cos_0)/pow(cos_0,3);
	//omega_0 = (-b+pow(b*b-4.*a*c,0.5))/2./a;

	//omega_0 = 0.;
	omega_0 = 1./pow(1-cos_0*cos_0,0.5);

	fprintf(stderr,"cos_0=%g, omega_0=%g\n",cos_0,omega_0);

	flag=0;
	if(omega_0 > 1./pow(1-cos_0*cos_0,0.5)){flag=1; omega_0 = pow((1-cos_0*cos_0)/pow(cos_0,3)/para_x/(1+epsilon*epsilon)-1.,0.5);}
	fprintf(stderr,"omega_0=%g\n",omega_0);

//	delta = 81.*a*a-12;
//	r_0 = r_d* (cal_ini(1,0)*cal_ini(1,0));
//	cos_0=cal_ini(1,0);

//	fprintf(stderr,"a=%g r_d=%g r_0=%g cos = %g\n",a, r_d, r_0, cal_ini(1,0));


	//find boundary conditions
//	alpha = dMinfdt*pow(r_0,1.5)*cs*cs*(1./cos_0/cos_0-1.)/2./pow(G*Mstar,0.5)/dMwinddt/v_wind;
//	a = 5.*alpha/2./r_0;
//	b = 1/2./r_0/pow(1.-cos_0*cos_0,0.5)+cos_0*cos_0/2./pow(r_0,3)/pow(1.-cos_0*cos_0,1.5);
//	b = b * alpha;
//	c = 5.* alpha / 2. / r_0 -4.;
//	d = b;

//	A = b*b-3.*a*c;
//	B=b*c-9.*a*d;
//	C=c*c-3.*b*d;
//	delta = B*B-4.*A*C;

//	if(delta<0){
//		T = (2.*A*b-3.*a*B)/2./pow(A,1.5);
//		beta = acos(T);

//		y1 = 1./3./a*(-b+pow(A,0.5)*(cos(beta/3.)+1.732*sin(beta/3.)));
//		y1 = 1./3./a*(-b+pow(A,0.5)*(cos(beta/3.)-1.732*sin(beta/3.)));
//		y3 = (-b-2*pow(A,0.5)*cos(beta/3))/3./a;

//		fprintf(stderr,"T=%g beta=%g y1=%g y2=%g y3=%g\n",T,beta, y1,y2,y3);
//	}

//	fprintf(stderr,"a=%g b=%g c=%g d=%g A=%g B=%g C=%g delta=%g\n",a,b,c,d,A,B,C,delta);

	 //r=r0+alpha*z
	r_0 = cos_0*cos_0; //in unit of r_d

//	fprintf(stderr,"correction=%g\n", alpha * sigv[0]*4*pi*r_0/dMwinddt);

	rad[0] = r_0;
	omega[0] = omega_0;
	sigv[0]= 0;
	sigv2[0] = 0; //0.5*dMwinddt/2./pi/r_0/alpha;
	//vel[0] = alpha/pow(1.+alpha*alpha,0.5)*v_wind;

	//fprintf(stderr,"correction=%g\n", alpha * sigv[0]*4*pi*r_0/dMwinddt);

	rad[1] = r_0 + omega_0 * dx;
	//if(flag==0 ) omega[1] = omega_0 + dx * (pow(1.+omega_0*omega_0,2)/omega_0*(0.5*para_x*cos_0/(1.-cos_0*cos_0)*(1.+epsilon*epsilon)*(omega_0+1./pow(1.-cos_0*cos_0,1.5))*(1.+vdcs*vdcs/cos_0/cos_0/(1.+omega_0*omega_0)*pow(1.-omega_0*pow(1-cos_0*cos_0,0.5),2))+para_x*cos_0/(1.-cos_0*cos_0)/(1.+omega_0*omega_0)*vdcs*vdcs/cos_0/cos_0*(pow(omega_0,3)*(1-cos_0*cos_0)-3.*omega_0*omega_0*pow(1-cos_0*cos_0,0.5)+2.*omega_0*cos_0*cos_0+2.*pow(1-cos_0*cos_0,0.5)-1./pow(1-cos_0*cos_0,0.5))));
	//omega[1]=omega[0];
	omega[1] =  omega_0 + dx * (pow(1.+omega_0*omega_0,2)/omega_0*(0.5*para_x*cos_0/(1.-cos_0*cos_0)*(1.+epsilon*epsilon)*(omega_0+1./pow(1.-cos_0*cos_0,1.5))-4.*omega_0/(1.+omega_0*omega_0)/cos_0/cos_0));

	//omega[1] = omega_0 + dx * (pow(1.+omega_0*omega_0,2)/omega_0*(0.5*para_x*cos_0/(1.-cos_0*cos_0)*(1.+epsilon*epsilon)*(omega_0+1./pow(1.-cos_0*cos_0,1.5))*(1.+vdcs*vdcs/cos_0/cos_0/(1.+omega_0*omega_0)*pow(1.-omega_0*pow(1-cos_0*cos_0,0.5),2))-4.*omega_0/(1.+omega_0*omega_0)/cos_0/cos_0));
	sigv[1] = sigv[0] + (1./pow(cos_0,4) + para_x*vwcs/cos_0/(1.-cos_0*cos_0)*(1.+epsilon*epsilon)*alpha/vwcs*pow(1+omega_0*omega_0,0.5)) * dx;

//	sigv[1] = sigv[0] + (dMwinddt/4./pi/r_0/r_0-alpha*sigv[0]/r_0) * dx;
	sigv2[1] = sigv2[0]+omega_0/pow(1+omega_0*omega_0,0.5)/pow(cos_0,4) * dx;	

	//a = 1;
	//b = 0;
	//delta = 81.*radius*radius*cos(theta)*cos(theta)+12*pow(radius-1,3);

	for(i=2;i<Nh;i++){
		rk4(f1,f2,f3,f4,dx, x0 + dx * (i-1),rad[i-1],omega[i-1],sigv[i-1],sigv2[i-1],i);
		fterm1[i]=term1;
		fterm2[i]=term2;
//		if(rad[i-1]<1) rad[i-1]=1;
	}

	for(i=0;i<Nh;i++){
		if(rad[i]<0) break;
		
		xi = x0 + dx * i;
		term3=para_x/pow(xi*xi+r_0*r_0,0.75)/pow(1+xi/pow(xi*xi+r_0*r_0,0.5)/cal_theta(pow(xi*xi+r_0*r_0,0.5),atan(r_0/(xi+1e-20))),0.5)/(xi/pow(xi*xi+r_0*r_0,0.5)/2./cal_theta(pow(xi*xi+r_0*r_0,0.5),atan(r_0/(1e-20+xi)))+1./pow(xi*xi+r_0*r_0,0.5)*pow(cal_theta(pow(xi*xi+r_0*r_0,0.5),atan(r_0/(1e-20+xi))),2))*(1.+epsilon*epsilon);
		fprintf(stdout, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n", x0 + dx * i, rad[i], omega[i], sigv[i], sigv2[i]/sigv[i], -sigv2[i]/sigv[i]*(omega[i]/pow(1+omega[i]*omega[i],0.5)),fterm1[i],fterm2[i],term3, ram[i], v2[i], vr[i],vth[i]);//vr,vth in unit of cs
	}
	return 0;
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
        term1 = para_x/pow(x*x+y1*y1,0.75)/pow(1+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(x+1e-20))),0.5)/(x/pow(x*x+y1*y1,0.5)/2./cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))+1./pow(x*x+y1*y1,0.5)*pow(cal_theta(pow(x*
x+y1*y1,0.5),atan(y1/(1e-20+x))),2))*(1.+epsilon*epsilon);
	term2 = pow(y1-y2*x,2)/pow(x*x+y1*y1,2)/(1+y2*y2);
 
	//term3 = para_x/pow(x*x+y1*y1,0.75)/pow(1+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(x+1e-20))),0.5)/(x/pow(x*x+y1*y1,0.5)/2./cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))+1./pow(x*x+y1*y1,0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),2))*(1.+epsilon*epsilon);

        //fprintf(stderr,"%g %g\n",term1,term2);

	vr[num] = -pow(vdcs*vdcs/pow(x*x+y1*y1,0.5),0.5)*pow(1.+x/pow(x*x+y1*y1,0.5)/cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),0.5);

	vth[num] = -pow(vdcs*vdcs/pow(x*x+y1*y1,0.5),0.5)*pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))) + x/pow(x*x+y1*y1,0.5),0.5)/pow(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x))),0.5)/(y1/pow(x*x+y1*y1,0.5))*(cal_theta(pow(x*x+y1*y1,0.5),atan(y1/(1e-20+x)))-x/pow(x*x+y1*y1,0.5));

	vram = vr[num]*(x*y2-y1)/pow(x*x+y1*y1,0.5)/pow(1+y2*y2,0.5)+vth[num]*(y1*y2+x)/pow(x*x+y1*y1,0.5)/pow(1+y2*y2,0.5);

	pram = vram *vram;

	//fprintf(stderr,"z=%g,R=%g vr=%g vth=%g vram = %g pram=%g\n",x,y1,vr[num],vth[num],vram,pram);

	ram[num] = pram;
	v2[num] = vr[num]*vr[num]+vth[num]*vth[num];

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

