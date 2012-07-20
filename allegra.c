/* compile with cc -g -O3 -Wall -Wextra -ffast-math -std=c99 allegra.c -lm -o allegra */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#undef CALLCOUNT
#define RANGE 40
#define ITERS 30

/* complex number typedefs and assorted effluvia */

typedef struct
{
  double re;
  double im;
} cx;

static const cx origin = { 0.0, 0.0 };
static const cx ai = { 0.0, 1.0 };
static const cx minusi = {0.0, -1.0};
static const cx rone = {1.0, 0.0};
static const cx rtwo = {2.0, 0.0};

#if 1
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
#else
static cx add(cx m, cx n)
{
  complex o = (m.re + m.im * I) + (n.re + n.im * I);
  cx out;
  out.re = creal(o); out.im = cimag(o);
  return out;
}

static cx cdiff(cx m, cx n)
{
  complex o = (m.re + m.im * I) - (n.re + n.im * I);
  cx out;
  out.re = creal(o); out.im = cimag(o);
  return out;
}

static cx mult(cx m, cx n)
{
  complex o = (m.re + m.im * I) * (n.re + n.im * I);
  cx out;
  out.re = creal(o); out.im = cimag(o);
  return out;
}

static cx rmult(double u, cx m)
{
  complex o = (m.re + m.im * I) * u;
  cx out;
  out.re = creal(o); out.im = cimag(o);
  return out;
}

static cx cdiv(cx m, cx n)
{
  complex o = (m.re + m.im * I) / (n.re + n.im * I);
  cx out;
  out.re = creal(o); out.im = cimag(o);
  return out;
}
#endif


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

#ifdef CALLCOUNT
static int exp_count, exp_re_count, exp_im_count;
static int pow_count, pow_a_re_count, pow_a_im_count, pow_b_re_count, pow_b_im_count;
#endif


#if 0
/* m.im == 0 is a common case for this function, optimize for it. */
static cx expc(cx m)
{
#ifdef CALLCOUNT
  exp_count++;
  if (m.re == 0) exp_re_count++;
  if (m.im == 0) exp_im_count++;
#endif

  cx out;

  if (m.im == 0) {
    out.re = exp(m.re);
    out.im = 0;
  } else {
    out.re = exp(m.re) * cos(m.im);
    out.im = exp(m.re) * sin(m.im);
  }
  return out;
}
#else
static cx expc(cx m)
{
  complex retval = cexp(m.re + m.im * I);
  cx out;

  out.re = creal(retval); out.im = cimag(retval);
  return out;
}
#endif

#if 0
/* bg.re == 0 is a common case for this function, but optimizing for that doesn't
   help significantly. */
static cx powc(cx ag, cx bg)
{  
#ifdef CALLCOUNT
  pow_count++;
  if (ag.re == 0) pow_a_re_count++;
  if (ag.im == 0) pow_a_im_count++;
  if (bg.re == 0) pow_b_re_count++;
  if (bg.im == 0) pow_b_im_count++;
#endif

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
  complex retval = cpow(ag.re + ag.im * I, bg.re + bg.im * I);
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
  //return hsv2rgb(atan2(q.im,q.re)/pi ,1.0,1.0);
  //return hsv2rgb((atan2(q.im,q.re)/pi+1)/2.0 ,1.0,1.0);
}

/* here's the actual newton method def -- I"m going to be using
   thirty iterations at the moment, because that seems to get fine
   results with mpmath */
 
#ifdef CALLCOUNT
static int newtzeros = 0;
static int newtnonzeros = 0;
static int newtfullcount = 0;
static int newtlatestexit = 0;
#endif

static cx newt(cx z, cx q)
{
  /* these don't really need to be defined as global variables, so
     I'll define them locally */
  cx current;
  current = z;

  /* precalculate some stuff.  blorf could be hoisted out of newt() 
     but it's already 2 inner loops up from where it was */
  cx qpart[2 * RANGE + 1];
  cx blorf[2 * RANGE + 1];
  int n;
  for (n = -RANGE; n <= RANGE; n++) {
    cx ponent;
    ponent.re = n*n;
    ponent.im = 0.0;
    qpart[n+RANGE] = powc(q,ponent);
    blorf[n+RANGE] = rmult(2*n,ai);
  }

  int ix; 
#ifdef CALLCOUNT
  int zerocount = 0;
#endif
  for(ix = 0; ix < ITERS; ix++)
    {
      const cx zpart = expc(current);
      cx sum = origin;
      cx psum = origin;
      int n;
      for(n=-RANGE;n<=RANGE;n++)
	{
	  /* if qpart is 0, none of the rest of this matters */
	  /* XXX this could be moved up to the precalculation loop to control the iterations we do here */
	  if (qpart[n+RANGE].re == 0 && qpart[n+RANGE].im == 0) continue;
	  /* calculate jtheta3 */
	  cx zpart_n = powc(zpart, blorf[n+RANGE]);
	  sum = add(sum,mult(qpart[n+RANGE],zpart_n));
	  /* calculate pjtheta3 from that */
	  zpart_n = mult(blorf[n+RANGE],zpart_n);
	  psum = add(psum,mult(qpart[n+RANGE],zpart_n));
	}
      cx next = cdiff(current,cdiv(sum, psum));
      // If the iteration has converged, bail.
      // If next has NANs, these comparisons will be true, and we'll
      // fall out here, which is what we want anyway
      if (current.re == next.re && current.im == next.im) {
#ifdef CALLCOUNT
	if (ix > newtlatestexit) newtlatestexit = ix;
	if (!zerocount) {
	  // fprintf(stderr, "iter %d\n", ix);
	  newtzeros++;
	  zerocount++;
	}
	if (memcmp(&current, &next, sizeof current)) {
	  fprintf(stderr, "%a %a %a %a\n", current.re, current.im, next.re, next.im);
	}
#endif
	// next may have NANs, use it instead of current, to be similar to old code
	current = next;
	break;
      }
#ifdef CALLCOUNT
      else if (zerocount) newtnonzeros++;
      if (ix == ITERS) newtfullcount++;
#endif
      current = next;
    }
  return current;
}

int main(int argc, char *argv[])
{
  if (argc != 2 && argc != 3) {
    fprintf(stderr, "%s displaysize [width]\n", argv[0]);
    exit(1);
  }

  int displaysize = atoi(argv[1]);
  if (displaysize % 2 || displaysize < 2) {
    fprintf(stderr, "%s: displaysize needs to be even and greater than 0\n",
	    argv[0]);
    exit(1);
  }
  int width = displaysize;
  if (argc == 3) {
    width = atoi(argv[2]);
    if (width % 2 || width < 2) {
      fprintf(stderr, "%s: width needs to be even and greater than 0\n",
	      argv[0]);
      exit(1);
    }
  }



  double halfdisplay = displaysize / 2.0, halfwidth = width / 2.0;

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

  printf("P3\n%d %d\n255\n", width, displaysize);
  for(a=0;a<displaysize;a++)
    {
      for(b=0;b<width;b++)
	{
	  z.im = 4.3 * (a-halfdisplay) / halfdisplay;
	  z.re = 4.3 * (b-halfwidth) / halfwidth;
	  /* I'm going to keep testq constant now */
	  // display[a][b] = argcolor(newt(z,testq));
	  color c = argcolor(newt(z, testq));
	  printf("%d %d %d\n", c.red, c.green, c.blue);
	}
    }

#ifdef CALLCOUNT
  fprintf(stderr, "exp %d re %d im %d\n"
	  "pow %d a re %d a im %d b re %d b im %d\n",
	  exp_count, exp_re_count, exp_im_count,
	  pow_count, pow_a_re_count, pow_a_im_count, pow_b_re_count, pow_b_im_count);
  fprintf(stderr, "full iter count %d zeros %d following nonzeros %d latest exit %d\n", newtfullcount, newtzeros, newtnonzeros, newtlatestexit);
#endif
  return 0;
}
