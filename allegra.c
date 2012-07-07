#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>


/* complex number typedefs and assorted effluvia */

typedef struct
{
  double re;
  double im;
} cx;

static cx add(cx m, cx n)
{
  cx out;
  out.re = m.re + n.re;
  out.im = m.im + n.im;
  return out;
}

static cx cdiff(cx m, cx n)
{
  cx out;
  out.re = m.re - n.re;
  out.im = m.im - n.im;
  return out;
}

static cx mult(cx m, cx n)
{
  cx out;
  out.re = m.re*n.re - m.im *n.im;
  out.im = m.im*n.re + m.re * n.im;
  return out;
}

static cx rmult(double u, cx m)
{
  cx out;
  out.re = u*m.re;
  out.im = u*m.im;
  return out;
}

static cx jcon(cx m)
{
  cx out;
  out.re = m.re;
  out.im = -m.im;
  return out;
}

static double norm2(cx m)
{
  double out;
  out = (m.re*m.re + m.im*m.im);
  return out;
}


static cx recip(cx m)
{
  cx out;
  out = rmult(1/norm2(m),jcon(m));
  return out;
}

static cx cdiv(cx m, cx v)
{
  cx out;
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

#if 0
static cx expc(cx m)
{
  cx out;

  out.re = exp(m.re) * cos(m.im);
  out.im = exp(m.re) * sin(m.im);
  return out;
}
#else
static cx expc(cx m)
{
  double complex retval = cexp(m.re + m.im * I);
  cx out;

  out.re = creal(retval); out.im = cimag(retval);
  return out;
}
#endif

#if 1
static cx powc(cx ag, cx bg)
{  
  cx out;
  cx mesp, frim;
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
#else
static cx powc(cx ag, cx bg)
{
  cx out;
  double complex retval = cpow(ag.re + ag.im + I, bg.re + bg.im + I);
  out.re = creal(retval); out.im = cimag(retval);
  return out;
}
#endif

typedef struct  {
  int red;
  int green;
  int blue;
} color;

static color hsv2rgb(double hue, double saturation, double value)
{
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

static double  pi = 3.1415926535;

static color argcolor(cx q)
{
  return hsv2rgb((atan(q.im/q.re)+pi/2.0)/pi ,1.0,1.0);
}

/* color argcolor(cx q) */
/* { */
/*   return hsv2rgb(atan2(q.im,q.re)/pi ,1.0,1.0); */
/* } */

static const cx origin = { 0.0, 0.0 };
static const cx ai = { 0.0, 1.0 };

/* here's the actual newton method def -- I"m going to be using
   thirty iterations at the moment, because that seems to get fine
   results with mpmath */
 
static cx newt(cx z, cx q)
{
  /* these don't really need to be defined as global variables, so
     I'll define them locally */
  cx current;
  current = z;

  /* precalculate some stuff.  blorf could be hoisted out of newt() 
     but it's already 2 inner loops up from where it was */
  cx qpart[80];
  cx blorf[80];
  int n;
  for (n = -40; n < 40; n++) {
    cx ponent;
    ponent.re = n*n;
    ponent.im = 0.0;
    qpart[n+40] = powc(q,ponent);
    blorf[n+40] = rmult(2*n,ai);
  }

  int ix; 
  for(ix = 0; ix < 30 ; ix++)
    {
      const cx zpart = expc(current);
      cx sum = origin;
      cx psum = origin;
      int n;
      for(n=-40;n<40;n++)
	{
	  /* calculate jtheta3 */
	  cx zpart_n = powc(zpart, blorf[n+40]);
	  sum = add(sum,mult(qpart[n+40],zpart_n));
	  /* calculate pjtheta3 from that */
	  zpart_n = mult(blorf[n+40],zpart_n);
	  psum = add(psum,mult(qpart[n+40],zpart_n));
	}
      current = cdiff(current,cdiv(sum, psum));
    }
  return current;
}

#define MAXDISPLAY 2048
static color display[MAXDISPLAY][MAXDISPLAY];
static const int maxdisplay = MAXDISPLAY;
#undef MAXDISPLAY

int main(int argc, char *argv[])
{
  if (argc != 2) {
    fprintf(stderr, "%s displaysize\n", argv[0]);
    exit(1);
  }

  int displaysize = atoi(argv[1]);
  if (displaysize % 2 || displaysize < 2 || displaysize > maxdisplay) {
    fprintf(stderr, "%s: displaysize needs to be even and between 2 and %d\n",
	    argv[0], maxdisplay);
    exit(1);
  }

  int halfdisplay = displaysize / 2;


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

  cx testq; 
  testq.re = 0.001;
  testq.im = -.3019;
#if 0
  cx testz;
  testz.re = 0.212;
  testz.im= 0.110;
  printf("theta3(testz,testq)= %.20f %.20f i\n", jtheta3(testz,testq).re, jtheta3(testz,testq).im);
#endif

  int a, b;
  cx z;

  for(a=0;a<displaysize;a++)
    {
      for(b=0;b<displaysize;b++)
	{
	  z.im = 4.3*((double) (a-halfdisplay))/(double)halfdisplay;
	  z.re = 4.3*((double) (b-halfdisplay))/(double)halfdisplay;
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
