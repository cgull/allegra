#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* complex number typedefs and assorted effluvia */

typedef struct{
  double re;
  double im;
} complex;

static inline complex add(complex m, complex n)
{
  complex out;
  out.re = m.re + n.re;
  out.im = m.im + n.im;
  return out;
}

static inline complex cdiff(complex m, complex n)
{
  complex out;
  out.re = m.re - n.re;
  out.im = m.im - n.im;
  return out;
}

static inline complex mult(complex m, complex n)
{
  complex out;
  out.re = m.re*n.re - m.im *n.im;
  out.im = m.im*n.re + m.re * n.im;
  return out;
}

static inline complex rmult(double u, complex m)
{
  complex out;
  out.re = u*m.re;
  out.im = u*m.im;
  return out;
}

static inline complex jcon(complex m)
{
  complex out;
  out.re = m.re;
  out.im = -m.im;
  return out;
}

static inline double norm2(complex m)
{
  double out;
  out = (m.re*m.re + m.im*m.im);
  return out;
}


static inline complex recip(complex m)
{
  complex out;
  out = rmult(1/norm2(m),jcon(m));
  return out;
}

static inline complex cdiv(complex m, complex v)
{
  complex out;
  out = mult(m,recip(v));
  return out;
}



/* we need to define a raw exponential -- that is something
   which takes a complex number and raises it to a complex power, and
   I need to do this at the moment without the benefit of libraries
   that might slow the computation down */

/* if the imaginary component is zero then all we have to do is
   to multiply exp(real)* (cos(imag) + i * sin(imag) -- if imag is
   nonzero we have to do some extra fiddling, and we should
   do it intelligently -- the less we actually have to deal with
   the power series for exp the better, and because of the n^2 term
   in the theta function definition, we're not going to need
   more than twenty or thirty terms before the difference between the
   actual value and the computer representation is small enough
   to make no difference to the phase of the complex number being
   plotted */

static inline complex expc(complex m)
{
  complex out;
  complex mesp, frim;
  double radius, theta;

  out.re = exp(m.re) * cos(m.im);
  out.im = exp(m.re) * sin(m.im);
  return out;
}


static inline complex powc(complex ag, complex bg)
{  
  complex out;
  complex mesp, frim;
  double radius, theta;
  /* get the proper polar form of the complex number */
  radius =  sqrt(ag.re*ag.re + ag.im*ag.im);
  theta = atan2(ag.im,ag.re);
  /* mesp gives R^(c+di) */
  mesp.re = pow(radius,bg.re)*cos(bg.im*log(radius));
  mesp.im = pow(radius,bg.re)*sin(bg.im*log(radius));
  /* frim gives e^(i theta (c+di)) */
  /* now since we already have the machinery
     for performing complex exponentiation (just exp), we
     can just call that here */
  frim.re = -1.0 * bg.im * theta;
  frim.im = bg.re * theta;
  frim = expc(frim);  
  out = mult(mesp,frim);
  return out;
}


typedef struct
{
  int red;
  int green;
  int blue;
} color;

color hsv2rgb(double hue, double saturation, double value)
{
  double qred;
  double qgreen;
  double qblue;
  int xred=0;
  int xgreen=0; 
  int xblue=0;
  int lexen;
  double nexel;
  double M;
  double N;
  double K;
  color out;
  if(saturation == 0)
    {
      // achromatic case
      xred = (int) (255.0 * value);
      xgreen = (int) (255.0 * value);
      xblue = (int) (255.0 * value);
    }
  else
    {
      if(hue >= 1.0){
	hue = 0.0;
      }
      else
	{
	  hue = hue * 6;
	}
  
      lexen = ((int) (hue));
      nexel = hue - ((double) (lexen));
 
      M = value * (1 - saturation);
      N = value * (1 - saturation * nexel);
      K = value * (1-saturation*(1-nexel));    
      if(lexen==0) { xred = (int) (255.0*value); xgreen = (int) (255.0*K); xblue = (int) (255.0*M);} 
      if(lexen==1) { xred = (int) (255.0*N); xgreen = (int) (255.0*value); xblue = (int) (255.0*M);} 
      if(lexen==2) { xred = (int) (255.0*M); xgreen = (int) (255.0*value); xblue = (int) (255.0*K);} 
      if(lexen==3) { xred = (int) (255.0*M); xgreen = (int) (255.0*N); xblue = (int) (255.0*value);} 
      if(lexen==4) { xred = (int) (255.0*K); xgreen = (int) (255.0*M); xblue = (int) (255.0*value);} 
      if(lexen==5) { xred = (int) (255.0*value); xgreen = (int) (255.0*M); xblue = (int) (255.0*N);}
    }
  out.red = xred;
  out.green = xgreen;
  out.blue = xblue;
  return out;
}

double  pi = 3.1415926535;

color argcolor(complex q)
{
  return hsv2rgb((atan(q.im/q.re)+pi/2.0)/pi ,1.0,1.0);
}

/* color argcolor(complex q) */
/* { */
/*   return hsv2rgb(atan2(q.im,q.re)/pi ,1.0,1.0); */
/* } */

static const complex origin = { 0.0, 0.0 };
static const complex ai = { 0.0, 1.0 };

complex jtheta3(complex z, complex q)
{
  complex out;
  complex ponent;
  complex qpart;
  complex zpart;
  complex sum;
  sum = origin;
  int n;
  for(n=-40;n<40;n++)
    {
      ponent.re = n*n;
      ponent.im = 0.0;
      qpart = powc(q,ponent);
      zpart = expc(z);
      zpart = powc(zpart, rmult(2*n,ai));
      sum = add(sum,mult(qpart,zpart));
    }
  return sum;
}


/* we have to do this again, because there's a term
   of 2ni for each factor of the theta function, so we
   can't just have a switch for it in the first one */

complex pjtheta3(complex z, complex q)
{
  complex out;
  complex ponent;
  complex qpart;
  complex zpart;
  complex sum;
  complex blorf;
  sum = origin;
  int n;
  for(n=-40;n<40;n++)
    {
      ponent.re = n*n;
      ponent.im = 0.0;
      qpart = powc(q,ponent);
      zpart = expc(z);
      /* the next line explicitly calculates the derivative */
      blorf = rmult(2*n,ai);
      zpart = powc(zpart, rmult(2*n,ai));
      zpart = mult(blorf,zpart);
      sum = add(sum,mult(qpart,zpart));
    }
  return sum;
}

/* here's the actual newton method def -- I"m going to be using
   thirty iterations at the moment, because that seems to get fine
   results with mpmath */
 
complex newt(complex z, complex q)
{
  int ix; 
  /* these don't really need to be defined as global variables, so
     I'll define them locally */
  complex current, next;
  current = z;
  for(ix = 0; ix < 30 ; ix++)
    {
      next = cdiff(current,cdiv(jtheta3(current,q),pjtheta3(current,q)));
      current = next;
    }
  return current;
}

int main(int argc, char *argv[])
{

  static const int displaysize = 40;


  color display[displaysize][displaysize];
  int i,j;


  color gray;
  gray.red = 128;
  gray.green = 128;
  gray.blue = 128;
  for(i=0;i<displaysize;i++)
    {
      for(j=0;j<displaysize;j++)
	{
	  display[i][j] = gray;
	}
    }

  complex testq; 
  testq.re = 0.001;
  testq.im = -.3019;
#if 0
  complex testz;
  testz.re = 0.212;
  testz.im= 0.110;
  printf("theta3(testz,testq)= %.20f %.20f i\n", jtheta3(testz,testq).re, jtheta3(testz,testq).im);
#endif

  int a, b;
  complex z;

  for(a=0;a<displaysize;a++)
    {
      for(b=0;b<displaysize;b++)
	{
	  z.im = 4.3*((double) (a-256))/256.0;
	  z.re = 4.3*((double) (b-256))/256.0;
	  /* I'm going to keep testq constant now */
	  display[a][b] = argcolor(newt(z,testq));
	}
    }

  printf("P3\n%d %d\n255\n", displaysize, displaysize);
  for(a=0;a<displaysize;a++)
    {
      for(b=0;b<displaysize;b++)
	{
	  printf("%d %d %d\n", display[a][b].red, display[a][b].green, display[a][b].blue);
	}
    }

  return 0;
}
