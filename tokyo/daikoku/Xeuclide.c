//date : 20210325（ver.1.0）
//date      :  20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229,20191230
//auther    : the queer who thinking about cryptographic future
//code name :  一変数多項式演算ライブラリのつもり
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : majer release (ver 1.0)
// Niederreiter Cryotosysytem by EEA decoding

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <assert.h>
#include <execinfo.h>

#include "8192.h"
//#include "4096.h"
#include "global.h"
#include "struct.h"
#include "debug.c"
#include "chash.c"
#include "lu.c"
#include "sha3.c"
#include "inv_mat.c"
//#include "golay.c"

#define TH omp_get_max_threads()

extern unsigned long xor128(void);
extern int mlt(int x, int y);
extern int mltn(int n, int a);
extern void makeS();

//#pragma omp threadprivate(mat)
//シンドロームのコピー
unsigned short sy[K] = {0};

//Goppa多項式
static unsigned short g[K + 1] = {1, 0, 0, 0, 1, 0, 1};
//{1,0,1,1};
MTX BB = {0};
MTX H = {0};

unsigned int AA = 0, B = 0; //, C = 0, A2 = 0;

//有限体の元の逆数
unsigned short
oinv(unsigned short a)
{
  int i;

  if (a == 0)
    return -1;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(fg[a], i)] == 1)
      return (unsigned short)i;
  }

  printf("no return \n");
  //  exit (1);
}

//aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
  int i;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(fg[a], fg[i])] == b)
      break;
  }
  return i;
}

//OP型からベクトル型への変換
vec o2v(OP f)
{
  vec a = {0};
  int i;

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0 && f.t[i].n < DEG)
      a.x[f.t[i].n] = f.t[i].a;
  }

  return a;
}

//ベクトル型からOP型への変換
OP v2o(vec a)
{
  int i, j = 0;
  OP f = {0};

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      f.t[j].n = i;
      f.t[j++].a = a.x[i];
    }
  }

  return f;
}

//停止コマンド
void wait(void)
{
  int a;                                     // 読み込む変数はローカルに取るべし
  printf(" (enter number and hit return) "); // 何か表示させたほうが良いだろう
  fflush(stdout);                            // just in case
  scanf("%d", &a);                           // fgets(line, LINESIZE, stdin); という手も
}

//OP型を正規化する
OP conv(OP f)
{
  vec v = {0};
  OP g = {0};

  v = o2v(f);
  g = v2o(v);

  return g;
}

//多項式を表示する（OP型）
void oprintpol(OP f)
{
  int i, n;

  f = conv(f);
  n = odeg(f);
  printf("n=%d\n", n);
  printf("terms=%d\n", terms(f));
  printf("deg=%d\n", odeg(f));

  for (i = n; i > -1; i--)
  {
    if (f.t[i].a > 0)
      printf("%ux^%u+", f.t[i].a, f.t[i].n);
  }
}

void op_print_raw(const OP f)
{
  puts("op_print_raw:");
  for (int i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
      printf("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
  }
}

bool op_verify(const OP f)
{
  bool end = false;
  unsigned short n_max = 0;
  for (int i = 0; i < DEG; i++)
  {
    if (end && (f.t[i].n != 0 || f.t[i].a != 0))
    {
      op_print_raw(f);
      printf("found data after end: i=%d\n", i);
      print_trace();
      fflush(stdout);
      return false;
    }
    if (f.t[i].a == 0)
    {
      end = true;
      continue;
    }
    if (f.t[i].n + 1 <= n_max)
    {
      op_print_raw(f);
      printf("found invalid order: i=%d\n", i);
      print_trace();
      fflush(stdout);
      return false;
    }
    n_max = f.t[i].n + 1;
  }
  return true;
}

OP norm(OP f)
{
  OP h = {0};
  int i;

  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
    {
      //      h.t[f.t[i].n].n=f.t[i].n;
      h.t[f.t[i].n].a = f.t[i].a;
    }
  }

  // exit(1);

  return h;
}

//20200816:正規化したいところだがうまく行かない
//多項式の足し算
OP oadd(OP f, OP g)
{
  f = conv(f);
  g = conv(g);
  assert(op_verify(f));
  assert(op_verify(g));

  vec a = {0}, b = {0}, c = {0};
  int i, j, k, l = 0;
  OP h = {0}, f2 = {0}, g2 = {0};

  a = o2v(f);
  b = o2v(g);

  //k=deg(o2v(f));
  //l=deg(o2v(g));

  for (i = 0; i < DEG; i++)
  {
    c.x[i] = a.x[i] ^ b.x[i];
  }
  h = v2o(c);
  //h=conv(h);
  assert(op_verify(h));
  return h;
}

//項の順序を降順に揃える
OP sort(OP f)
{
  oterm o = {0};
  int i, j, k;

  k = terms(f);
  for (i = 0; i < k + 1; ++i)
  {
    for (j = i + 1; j < k + 1; ++j)
    {
      if (f.t[i].n > f.t[j].n)
      {
        o = f.t[i];
        f.t[i] = f.t[j];
        f.t[j] = o;
      }
    }
  }

  return f;
}

//リーディングタームを抽出(OP型）
oterm oLT(OP f)
{
  int i, k, j;
  oterm s = {0};

  k = terms(f);
  s = f.t[0];
  for (i = 0; i < k + 1; i++)
  {
    //printf("a=%d %d\n",f.t[i].a,f.t[i].n);
    if (f.t[i].a > 0)
    {
      printf("in LT=%d %d\n", s.a, s.n);
      for (j = i; j < k + 1; j++)
      {
        if (s.n < f.t[j].n)
        {
          s.n = f.t[j].n;
          s.a = f.t[j].a;
        }

        //  else{
        // t=s;
        // }
      }
    }
  }
  //  exit(1);

  return s;
}

//多項式を足し算する（OP型）
OP add(OP f, OP g)
{
  //  vec a={0},b={0},c={0};
  unsigned long long int i, j, n1 = 0, n2 = 0, m1 = 0, count = 0;
  OP h = {0};
  oterm o1 = {0}, o2 = {
                      0};

  n1 = terms(f);
  printf("n1=%d\n", n1);
  n2 = terms(g);
  printf("n2=%d\n", n2);
  if (n1 > n2)
  {
  }

  oprintpol(f);
  printf(" fff==============\n");
  oprintpol(g);
  printf(" ggg==============\n");
  o1 = oLT(f);
  o2 = oLT(g);
  printf("LTadd==%d %d\n", o1.n, o2.n);
  m1 = n1 + n2;
  printf("m1=%d\n", m1);
  // exit(1);

  for (i = 0; i < n1 + 1; i++)
  {
    for (j = 0; j < n2 + 1; j++)
    {
      if (f.t[i].n == g.t[j].n && g.t[j].a > 0 && f.t[i].a > 0)
      {
        o1 = oLT(f);
        o2 = oLT(g);
        printf("LT==%d %d\n", o1.n, o2.n);
        printf("f.n==%d %d %d %d\n", f.t[i].n, g.t[j].n, i, j);
        f.t[i].a = 0;
        g.t[j].a = 0;
      }
    }
  }
  for (i = 0; i < n2 + 1; i++)
  {
    if (g.t[i].a > 0)
    {
      h.t[count++] = g.t[i];
      g.t[i].a = 0;
    }
  }
  for (i = 0; i < n1 + 1; i++)
  {
    if (f.t[i].a > 0)
    {
      h.t[count++] = f.t[i];
      f.t[i].a = 0;
    }
  }

  h = sort(h);
  /*
     for (i=0; i<count; ++i) {
     for (j=i+1; j<count; ++j) {
     if (h.t[i].n > h.t[j].n) {
     oo =  h.t[i];
     h.t[i] = h.t[j];
     h.t[j] = oo;
     }
     }
     }
   */
  h = conv(h);
  if (odeg(h) > 0)
    oprintpol(h);
  printf(" addh==============\n");
  //   exit(1);

  return h;
}

//多項式を項ずつ掛ける
OP oterml(OP f, oterm t)
{
  f = conv(f);
  assert(op_verify(f));
  int i, k, j;
  OP h = {0};
  vec test;
  unsigned int n;

  //f=conv(f);
  //k = deg (o2v(f));
  j = 0;
  for (i = 0; i < DEG; i++)
  {
    h.t[i].n = f.t[i].n + t.n;
    h.t[i].a = gf[mlt(fg[f.t[i].a], fg[t.a])];
  }

  h = conv(h);
  assert(op_verify(h));
  return h;
}

//多項式の掛け算
OP omul(OP f, OP g)
{
  f = conv(f);
  g = conv(g);
  assert(op_verify(f));
  assert(op_verify(g));
  int i, count = 0, k, l;
  oterm t = {0};
  OP h = {0}, e = {0}, r = {0};
  vec c = {0};

  k = odeg(f);
  l = odeg(g);
  if (l > k)
  {
    k = l;
  }

  for (i = 0; i < k + 1; i++)
  {
    t = g.t[i];
    e = oterml(f, t);
    h = oadd(h, e);
  }
  assert(op_verify(h));
  return h;
}

//リーディグタームを抽出(default)
oterm LT(OP f)
{
  int i, k;
  oterm t = {0};

  //k = deg (o2v (f));
  for (i = 0; i < DEG; i++)
  {
    //printf("a=%d %d\n",f.t[i].a,f.t[i].n);
    if (f.t[i].a > 0)
    {
      t.n = f.t[i].n;
      t.a = f.t[i].a;
    }
  }

  return t;
}

//多項式の最後の項を抽出
oterm LT2(OP f)
{
  int i, k;
  oterm t = {0};

  t.n = f.t[0].n;
  t.a = f.t[0].a;

  return t;
}

//多項式を単行式で割る
oterm LTdiv(OP f, oterm t)
{
  oterm tt = {0}, s = {
                      0};

  tt = LT(f);
  if (tt.n < t.n)
  {
    s.n = 0;
    s.a = 0;
  }
  else if (tt.n == t.n)
  {
    s.n = 0;
    s.a = equ(t.a, tt.a);
  }
  else if (tt.n > t.n)
  {
    s.n = tt.n - t.n;
    s.a = equ(t.a, tt.a);
    //printf("%u\n",s.a);
  }
  else if (t.n == 0 && t.a > 0)
  {
    s.a = gf[mlt(fg[tt.a], oinv(t.a))];
    s.n = tt.n;
  }

  return s;
}

//モニック多項式にする
OP coeff(OP f, unsigned short d)
{
  int i, j, k;
  vec a, b;

  f = conv(f);
  k = odeg((f)) + 1;
  for (i = 0; i < k; i++)
    f.t[i].a = gf[mlt(fg[f.t[i].a], oinv(d))];

  return f;
}

//多項式を表示する(default)
void printpol(vec a)
{
  int i, n;

  n = deg(a);

  //printf ("baka\n");
  assert(("baka\n", n >= 0));

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      printf("%u", a.x[i]);
      //if (i > 0)
      printf("x^%d", i);
      //if (i > 0)
      printf("+");
    }
  }
  //  printf("\n");

  return;
}

//多項式の剰余を取る
OP omod(OP f, OP g)
{
  int i = 0, j, n, k;
  OP h = {0}, e = {
                  0};
  oterm a, b = {0}, c = {0};

  n = LT(g).n;

  //  assert (("baka^\n", LT (f).n != 0));

  //  assert (("baka(A)\n", LT (g).n != 0));

  if (LT(f).n < LT(g).n)
  {
    //    exit(1);
    return f;
  }

  //printf ("in omod\n");
  //exit(1);

  k = LT(g).n;
  b = LT(g);

  assert(("double baka\n", b.a != 0 && b.n != 0));
  while (LT(f).n > 0 && LT(g).n > 0)
  {

    c = LTdiv(f, b);
    h = oterml(g, c);
    f = oadd(f, h);
    if (odeg((f)) == 0 || odeg((g)) == 0)
    {
      //      printf("blake1\n");
      break;
    }

    if (c.n == 0 || b.n == 0)
      break;
  }

  return f;
}

/*
//多項式の剰余を取る
OP
omod (OP f, OP g)
{
  int i = 0, j, n, k;
  OP h = { 0 }, tmp,e = {
    0
  };
  oterm a, b = { 0 }, c = {
    0
  };
  n = LT (g).n;
//printpol(o2v(f));
//printf(" ====@f\n");
//printpol(o2v(g));
//printf(" ====@g\n");
  assert (("baka^\n", LT (f).n != 0));
  assert (("baka(A)\n", LT (g).n != 0));
		 
  if (LT (f).n < LT (g).n)
    {
      printf("us\n");
      //wait();
      //    exit(1);
      return f;
    }
  //printf ("in omod\n");
  //exit(1);
  k = LT (g).n;
  b = LT (g);
  assert (("double baka\n", b.a != 0 && b.n != 0));
  while ((LT (f).n > 0 && LT (g).n > 0))
    {
      //printf("b.a=%d b.n=%d\n",b.a,b.n);
      c = LTdiv (f, b);
      h = oterml (g, c);
//      printpol(o2v(h));
//      printf(" ===@h\n");
    //  printpol(o2v(f));
  //    printf(" ===@fx\n");
      
      f = oadd (f, h);
      if (deg (o2v(f)) == 0 || deg (o2v(g)) == 0)
	{
//    printf("f.a=%d f.n=%d\n",LT(f).a, deg(o2v(f)));
//    printf("g.a=%d g.n=%d\n",LT(g).a, deg(o2v(g)));
//	  printf ("blake1\n");
	  break;
	}
      if (c.n == 0 || b.n == 0)
	      break;
//  b=LT(h);
    }
  return f;
}
*/

//多項式の商を取る
OP odiv(OP f, OP g)
{

  f = conv(f);
  g = conv(g);
  assert(op_verify(f));
  assert(op_verify(g));
  int i = 0, j, n, k;
  OP h = {0}, e = {0}, tt = {0};
  oterm a, b = {0}, c = {0};

  if (LT(f).n == 0 && LT(g).a == 0)
  {
    printf("baka^\n");
    //return f;
    exit(1);
  }
  if (LT(g).a == 0)
  {
    print_trace();
    exit(1);
  }
  if (LT(g).n == 0 && LT(g).a > 1)
    return coeff(f, LT(g).a);

  k = odeg(g);
  b = LT(g);
  if (b.a == 1 && b.n == 0)
    return f;
  if (b.a == 0 && b.n == 0)
  {
    printf("baka in odiv\n");
    exit(1);
  }
  if (odeg((f)) < odeg((g)))
  {
    return f;
    //  a=LT(f);
  }

  i = 0;
  while (LT(f).n > 0 && LT(g).n > 0)
  {
    c = LTdiv(f, b);
    assert(c.n < DEG);
    tt.t[i] = c;
    i++;

    h = oterml(g, c);

    f = oadd(f, h);
    if (odeg((f)) == 0 || odeg((g)) == 0)
    {
      //printf ("blake2\n");
      break;
    }

    if (c.n == 0)
      break;
  }

  // tt は逆順に入ってるので入れ替える
  OP ret = {0};
  int tt_terms = terms(tt);
  for (i = 0; i < tt_terms; i++)
  {
    ret.t[i] = tt.t[tt_terms - i - 1];
  }
  ret = conv(ret);
  assert(op_verify(ret));
  return ret;
}

//多項式のべき乗
OP opow(OP f, int n)
{
  int i;
  OP g = {0};

  g = f;

  for (i = 1; i < n; i++)
    g = omul(g, f);

  return g;
}



OP opowmod1(OP f, OP mod, int n) {
    vec v={0};
    OP ret;

     v.x[0]= 1;
    ret=v2o(v);
    while (n > 0) {
        if (n & 1) ret = omod(omul(ret , f),mod) ;  // n の最下位bitが 1 ならば x^(2^i) をかける
        f = omod(omul(f , f),mod);
        n >>= 1;  // n を1bit 左にずらす
    }
    return ret;
}


//多項式のべき乗余
OP opowmod(OP f, OP mod, int n)
{
  int i, j = 0,l;


  //printpol(o2v(mod));
  //printf(" =mod %d\n",LT(mod).n);
  //繰り返し２乗法
  for (i = 1; i <  n+1; i++)
  {
    f = omod(omul(f, f),mod);
  }

  return f;
}


//多項式の代入値
unsigned short
trace(OP f, unsigned short x)
{
  int i, d;
  unsigned short u = 0;

  d = deg(o2v(f));

  for (i = 0; i < d + 1; i++)
  {
    u ^= gf[mlt(fg[f.t[i].a], mltn(f.t[i].n, fg[x]))];
  }

  return u;
}

// invert of polynomial
OP inv(OP a, OP n)
{
  OP d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, tmp = {0};
  oterm b = {0};
  vec vv = {0}, xx = {
                    0};

  if (odeg((a)) > odeg((n)))
  {
    tmp = a;
    a = n;
    n = tmp;
    printf("baka_i\n");
    //exit (1);
  }
  if (LT(a).a == 0)
  {
    printf("in inv a ga 0\n");
    exit(1);
  }

  tt = n;

  d = n;
  x.t[0].a = 0;
  x.t[0].n = 0;
  s.t[0].a = 1;
  s.t[0].n = 0;
  while (odeg((a)) > 0)
  {
    if (odeg((a)) > 0)
      r = omod(d, a);
    if (LT(a).a == 0)
      break;
    if (LT(a).a > 0)
      q = odiv(d, a);

    d = a;
    a = r;
    t = oadd(x, omul(q, s));
    ////printpol (o2v (a));
    //printf ("\nin roop a==================%d\n", odeg ((a)));
    //printf ("\n");

    x = s;
    s = t;
  }
  // exit(1);
  //  if(LT(a).a>0){
  d = a;
  a = r;
  ////printpol (o2v (a));
  //printf ("\nin roop a|==================%d\n", odeg ((a)));
  //printf ("\n");

  x = s;
  s = t;

  ////printpol (o2v (d));
  //printf ("\nout1================\n");
  gcd = d; // $\gcd(a, n)$
  printpol(o2v(gcd));
  printf(" =========gcd\n");
  //exit(1);
  //printf ("\n");
  ////printpol (o2v (n));
  //printf ("\n");
  //printf ("out2===============\n");

  printf("before odiv\n");
  //w=tt;

  b = LT(w);
  ////printpol (o2v (w));
  //printf ("\nw=======%d %d\n", b.a, b.n);
  //w=tt;
  v = oadd(x, n);
  ////printpol (o2v (v));
  //printf ("\n");
  /*
     if (LT (v).a == 0)
     {
     printf ("v=============0\n");
     }
     printf ("d==============\n");
   */
  //  } //end of a>0
  w = tt;
  ////printpol (o2v (v));
  //printf ("\n");
  //printf ("ss==============\n");
  //       exit(1);
  // if(odeg((w))>0)
  if (LT(v).n > 0 && LT(w).n > 0)
  {
    u = omod(v, w);
  }
  else
  {
    //printpol (o2v (v));
    printf(" v===========\n");
    //printpol (o2v (x));
    printf(" x==0?\n");
    //printpol (o2v (n));
    printf(" n==0?\n");

    exit(1);
  }
  //caution !!
  if (LT(u).a > 0 && LT(d).a > 0)
  {
    u = odiv(u, d);
  }

  if (LT(u).a == 0 || LT(d).a == 0)
  {
    printf("inv div u or d==0\n");
    // exit(1);
  }
  //u=coeff(u,d.t[0].a);
  ////printpol (o2v (u));
  //printf ("\nu==================\n");
  if (LT(u).a == 0)
  {
    printf("no return at u==0\n");
    exit(1);
  }

  return u;
}

unsigned short
v2a(oterm a)
{
  int i, j;

  if (a.a == 0)
    return 0;

  //printf("aa=%d\n",a.a);
  for (j = 0; j < M; j++)
  {
    if (gf[j] == a.a && a.a > 0)
    {
      //printf("j==%d\n",j);
      return j - 1;
    }
  }
}

void printsage(vec a)
{
  int i, j, k;
  oterm b;

  printf("poly=");
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      b.a = a.x[i];
      b.n = i;
      j = v2a(b);
      //printf("%d,==ba\n",b.a);
      //printf ("X**%d+", i); //for GF2
      printf("B('a^%d')*X**%d+", j, i); //for GF(2^m)
    }
  }
}

OP gcd(OP a, OP b)
{
  OP r = {0}, h = {0}, tmp = {0};

  h.t[0].a = 1;
  h.t[0].n = 0;

  if (odeg(a) < odeg(b))
  {
    tmp = a;
    a = b;
    b = tmp;
  }
  /*
  printpol(o2v(a));
  printf(" ========f\n");
  printpol(o2v(b));
  printf(" ========g\n");
*/
  /* 自然数 a > b を確認・入替 */
  if (odeg(a) < odeg(b))
  {
    tmp = a;
    a = b;
    b = tmp;
  }

  r = omod(a, b);
  while (odeg(r) > 0)
  {
    a = b;
    b = r;
    r = omod(a, b);
    if (LT(r).a == 0)
      return b;
  }

  if (LT(r).a == 0)
  {
    return b;
  }
  else
  {
    //if(LT(r).a>0)
    return h;
  }
}

/*
OP gcd(OP a, OP b){
  OP r={0},h={0},tmp={0};
  h.t[0].a=1;
  h.t[0].n=0;
  if(odeg(a)<odeg(b)){
    tmp = a;
    a = b;
    b = tmp;
  }
  printsage(o2v(a));
  printf(" ========a\n");
  printsage(o2v(b));
  printf(" ========b\n");
  
  // 自然数 a > b を確認・入替 
  r = omod(a , b);
  while(odeg(r)>0){
    a = b;
    b = r;
    r = omod(a , b);
    if(LT(r).a==0){
      printf("e-!\n");
      printpol(o2v(a));
      printf(" =a\n");
      printpol(o2v(b));
      printf(" =b\n");
      printpol(o2v(r));
      printf(" =r\n");
      //wait();
      
      return b;
    }
  }
  
  if(LT(r).a==0){
    return b;
  }else if(LT(r).a==1){
    return h;
  }
}
*/

// gcd for patterson
OP zgcd(OP a, OP n)
{
  OP d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, rt = {0};
  oterm b = {0};
  vec vv = {0}, xx = {
                    0};

  if (odeg(a) > odeg(n))
  {
    rt = a;
    a = n;
    n = rt;
    printf("big is good\n");
    //exit (1);
  }
  if (LT(a).a == 0)
  {
    printf(" a ga 0\n");
    exit(1);
  }

  tt = n;

  d = n;
  x.t[0].a = 0;
  x.t[0].n = 0;
  s.t[0].a = 1;
  s.t[0].n = 0;
  while (LT(a).n > T)
  {

    r = omod(d, a);
    q = odiv(d, a);

    d = a;
    a = r;
    t = oadd(x, omul(q, s));

    x = s;
    s = t;
  }

  d = a;
  a = r;

  x = s;
  s = t;
  gcd = d; // $\gcd(a, n)$

  printpol(o2v(x));
  printf(" =======x\n");
  printpol(o2v(a));
  printf(" =======a\n");
  printpol(o2v(s));
  printf(" =======s\n");
  printpol(o2v(r));
  printf(" =======r\n");

  return x;
}

// GCD for decode
OP ogcd(OP xx, OP yy)
{
  OP tt;

  while (odeg(yy) > T - 1)
  {
    tt = omod(xx, yy);
    xx = yy;
    yy = tt;
  }

  printpol(o2v(yy));
  printf(" =========yy\n");
  printpol(o2v(tt));
  printf(" =========tt\n");

  return tt;
}

//test gcd
OP agcd(OP xx, OP yy)
{
  OP tt = {0}, tmp;

  if (deg(o2v(xx)) < deg(o2v(yy)))
  {
    tmp = xx;
    xx = yy;
    yy = tmp;
  }
  tt = omod(xx, yy);
  while (LT(tt).n > 0)
  {
    xx = yy;
    yy = tt;
    tt = omod(xx, yy);
  }

  return yy;
}

//error locater for decode
OP vx(OP f, OP g)
{
  OP h = {0}, ww = {
                  0};
  OP v[K] = {0}, vv = {
                     0};
  oterm a, b;
  int i, j;

  v[0].t[0].a = 0;
  v[0].t[1].n = 0;
  v[1].t[0].a = 1;
  v[1].t[1].n = 0;

  i = 0;

  while (1)
  {
    if (odeg((g)) == 0)
      break;
    h = omod(f, g);
    if (LT(g).a == 0)
      break;
    ww = odiv(f, g);
    v[i + 2] = oadd(v[i], omul(ww, v[i + 1]));
    f = g;
    g = h;

    vv = v[i + 2];

    if (odeg((vv)) == T)
      break;
    i++;
  }

  return vv;
}

OP sabun(OP f, OP g)
{
  OP h = {0}, ww = {
                  0};
  OP v[K * 2] = {0}, vv = {
                         0};
  oterm a, b;
  int i, j;

  v[0].t[0].a = 0;
  v[0].t[1].n = 0;
  v[1].t[0].a = 1;
  v[1].t[1].n = 0;

  i = 0;

  while (1)
  {
    if (odeg((g)) == 0)
    {
      printf("is a\n");
      break;
    }
    h = omod(f, g);
    if (LT(g).a == 0)
    {
      printf("is b\n");
      break;
    }
    ww = odiv(f, g);
    v[i + 2] = oadd(v[i], omul(ww, v[i + 1]));
    f = g;
    g = h;

    vv = v[i + 2];
    printpol(o2v(vv));
    printf(" deg%d\n", odeg(vv));
    if (odeg((vv)) == T * 2)
    {
      //printpol(o2v(vv));
      break;
      //exit(1);
    }
    i++;
  }
  //printf("baka\n");
  //exit(1);
  return vv;
}

//最終の項までの距離
int distance(OP f)
{
  int i, j = 0;

  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
      j = i;
  }

  return j;
}

//項の数
int terms(OP f)
{
  int i, count = 0;

  for (i = 0; i < DEG; i++)
    if (f.t[i].a > 0)
      count++;

  return count;
}

//多項式の次数(degのOP型)
int odeg(OP f)
{
  int i, j = 0, k;

  if (f.t[0].a == 0)
    return 0;

  //k=terms(f);
  for (i = 0; i < DEG; i++)
  {
    if (j < f.t[i].n && f.t[i].a > 0)
      j = f.t[i].n;
  }

  return j;
}

//０多項式かどうかのチェック
unsigned char
chk(OP f)
{
  int i, flg = 0;
  vec x = {0};

  x = o2v(f);
  for (i = 0; i < DEG; i++)
  {
    if (x.x[i] > 0)
    {
      flg = 1;
      return 1;
    }
  }
  if (flg == 0)
    return 0;

  exit(1);
}

//拡張ユークリッドアルゴリズム
EX xgcd(OP f, OP g)
{
  OP h = {0}, ww = {0}, *v, *u;
  oterm a, b;
  int i = 0, j, k, flg = 0;
  EX e = {0}, ee = {0};

  v = (OP *)malloc(sizeof(OP) * (DEG));
  u = (OP *)malloc(sizeof(OP) * (DEG));
  memset(v, 0, sizeof(OP) * DEG);
  memset(u, 0, sizeof(OP) * DEG);

  u[0].t[0].a = 1;
  u[0].t[0].n = 0;
  u[1].t[0].a = 0;
  u[1].t[0].n = 0;
  u[2].t[0].a = 1;
  u[2].t[0].n = 0;

  v[0].t[0].a = 0;
  v[0].t[0].n = 0;
  v[1].t[0].a = 1;
  v[1].t[0].n = 0;

  printpol(o2v(f));
  printf(" f===============\n");
  printpol(o2v(g));
  printf(" s===============\n");
  //  exit(1);

  k = 0;
  i = 0;
  while (LT(g).n > 0)
  //for (i = 0; i < T * 2; i++)
  {

    if ((odeg((g)) == 0 && LT(g).a == 0) || odeg((f)) == 0)
    {
      flg = 1;
      printf("v[%d]=%d skipped deg(g)==0!\n", i, odeg((v[i])));
      printpol(o2v(g));
      printf(" g========\n");
      e.d = f;
      e.v = v[i];
      e.u = u[i];

      free(v);
      free(u);
      //wait();
      return e;

      //exit (1);
      //return e;

      // break;
    }

    if (LT(g).n > 0)
      h = omod(f, g);

    if (LT(g).a > 0)
      ww = odiv(f, g);

    v[i + 2] = oadd(v[i], omul(ww, v[i + 1]));
    u[i + 2] = oadd(u[i], omul(ww, u[i + 1]));
    //printf ("i+1=%d %d %d g=%d\n", i + 1, odeg ((v[i])), T - 1, odeg ((g)));
    f = g;
    g = h;

    if (odeg(v[i]) > T - 2)
    {
      //printf("vaka\n");
      //wait();
      e.d = f;
      e.v = v[i];
      e.u = u[i];

      free(v);
      free(u);
      //wait();
      return e;

      //exit(1);
    }
    else if (deg(o2v(f)) == T - 1)
    {
      //  wait();
      printf("i=%d\n", i);
      //wait();
      e.d = f;
      e.v = v[i];
      e.u = u[i];

      free(v);
      free(u);

      return e;
    }
  }

  //printf ("i=%d\n", i);
  //wait();
  //oprintpol ((v[i]));
  printf("deg(v)=%d\n", odeg((v[i])));
  printf(" v=============\n");
  printf("deg(u)=%d\n", odeg((u[i])));
  //printpol (o2v (u[i]));
  printf(" u=============\n");
  printf("deg(f)=%d\n", odeg((f)));
  printf(" f=============\n");
  //exit(1);
  //  if(deg(v[i])==T-1){
  e.d = f;
  e.v = v[i];
  e.u = u[i];

  free(v);
  free(u);
  printf("end of fnc\n");
  //  wait ();

  return e;
}

OP init_pol(OP f)
{
  int i;

  for (i = 0; i < DEG; i++)
  {
    f.t[i].a = 0;
    f.t[i].n = 0;
  }

  return f;
}

//ランダム多項式の生成
static void
ginit(void)
{
  int j, count = 0, k = 0;
  unsigned short gg[K + 1] = {0};

  printf("in ginit\n");

  g[K] = 1;          //xor128();
  g[0] = rand() % N; //or N
  k = rand() % (K - 2);
  if (k > 0)
  {
    while (count < k)
    {
      printf("in whule\n");
      j = rand() % (K);
      if (j < K && j > 0 && g[j] == 0)
      {
        g[j] = rand() % N; //or N;
        count++;
      }
    }
  }

  for (j = 0; j < K + 1; j++)
    gg[j] = g[K - j];

  memcpy(g, gg, sizeof(g));
}

//ランダム置換の生成（Niederreoter 暗号における置換）
void random_permutation(unsigned short *a)
{
  int i, j, x;

  for (i = 0; i < N; i++)
  {
    a[i] = i;
  }
  for (i = 0; i < N - 2; i++)
  {
    j = (rand() % (N - 1 - i)) + i + 1;

    x = a[j];
    a[j] = a[i];
    a[i] = x;
  }
  if (a[N - 1] == N - 1)
  {
    a[N - 1] = a[N - 2];
    a[N - 2] = N - 1;
  }
}

/*
//配列から置換行列への変換
void P2Mat(unsigned short P[N])
{
  int i;
  for (i = 0; i < N; i++)
    AH.x[i][P[i]] = 1;
}
*/

unsigned short
b2B(unsigned short b[E])
{
  int i;
  unsigned short a = 0;

  for (i = E - 1; i > -1; i--)
    a ^= (b[E - i - 1] << i);

  return a;
}

//多項式の次数(default)
int deg(vec a)
{
  int i, n = 0, flg = 0;

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      n = i;
      flg = 1;
    }
  }
  if (flg == 0)
    return 0;

  return n;
}

//整数からベクトル型への変換
vec i2v(unsigned int n)
{
  vec v = {0};
  int i = 0;

  while (n > 0)
  {
    v.x[i++] = n % 2;
    n = (n >> 1);
  }

  return v;
}

//ベクトル型から整数への変換
unsigned int
v2i(vec v)
{
  unsigned int d = 0, i, e = 0;

  for (i = 0; i < deg(v) + 1; i++)
  {
    e = v.x[i];
    d ^= (e << i);
  }

  return d;
}

//配列からベクトル表現の多項式へ変換する
vec Setvec(int n)
{
  int i;
  vec v = {0};

  for (i = 0; i < n; i++)
  {
    v.x[n - 1 - i] = c[i];
  }

  return v;
}

void printvec(vec v)
{
  int i, j;

  for (i = 0; i < deg(v) + 1; i++)
  {
    //if (v.x[i] > 0)
    printf("%d:%d\n", i, v.x[i]);
  }
}

vec vadd(vec a, vec b)
{
  vec c = {0};
  int i, k;

  if (deg(a) >= deg(b))
  {
    k = deg(a) + 1;
  }
  else
  {

    k = deg(b) + 1;
  }
  for (i = 0; i < k; i++)
  {
    c.x[i] = a.x[i] ^ b.x[i];
  }
  // c.x[i] = a.x[i] ^ b.x[i];
  //  h = v2o (c);

  return c;
}

vec vmul(vec a, vec b)
{
  int i, j, k, l;
  vec c = {0};

  k = deg(a);
  l = deg(b);

  for (i = 0; i < k; i++)
  {
    for (j = 0; j < l; j++)
      if (a.x[i] > 0)
      {
        c.x[i + j] ^= gf[mlt(fg[a.x[i]], fg[b.x[j]])];
      }
  }

  return c;
}

//整数のべき乗
unsigned int
ipow(unsigned int q, unsigned int u)
{
  unsigned int i, m = 1;

  for (i = 0; i < u; i++)
    m *= q;

  printf("in ipow====%d\n", m);

  return m;
}

OP ww[T] = {0};

OP bib(int i, OP d)
{
  int id, j;

  OP t[T] = {0};
  //omp_set_num_threads(omp_get_max_threads());
  id = omp_get_thread_num();
  t[id] = d;

  //#pragma omp parallel for
  for (j = 0; j < T; j++)
  {
    // #pragma omp critical
    if (i != j)
    {
      t[id] = omul(t[id], ww[j]);
    }
  }

  return t[id];
}

//多項式の形式的微分
OP bibun(vec a)
{
  OP w[T * 2] = {0};
  OP l = {0}, t[T] = {0}, d = {0};
  int i, j, n, id;
  vec tmp = {0};

  n = deg(a);
  printf("n=%d\n", n);
  if (n == 0)
  {
    printf("baka8\n");
    //  exit(1);
  }
  memset(ww, 0, sizeof(ww));
  // #pragma omp parallel num_threads(8)
  for (i = 0; i < T; i++)
  {
    ww[i].t[0].a = a.x[i];
    ww[i].t[0].n = 0;
    ww[i].t[1].a = 1;
    ww[i].t[1].n = 1;

    ////printpol(o2v(w[i]));
  }
  //  exit(1);

  tmp.x[0] = 1;
  //
  d = v2o(tmp);

// omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel num_threads(omp_get_max_threads()) //num_threads(TH)
  {
    //#pragma omp parallel for
#pragma omp for schedule(static)
    for (i = 0; i < T; i++)
    {
      t[i] = bib(i, d);
    }
  }

  for (i = 0; i < T; i++)
    l = oadd(l, t[i]);

  return l;
}

//多項式の形式的微分
OP bibun_old(vec a)
{
  OP w[T * 2] = {0};
  OP l = {0}, t = {
                  0};
  int i, j, n, id;
  vec tmp = {0};

  n = deg(a);
  printf("n=%d\n", n);
  if (n == 0)
  {
    printf("baka8\n");
    //  exit(1);
  }

  //
  for (i = 0; i < T; i++)
  {
    w[i].t[0].a = a.x[i];
    w[i].t[0].n = 0;
    w[i].t[1].a = 1;
    w[i].t[1].n = 1;
    ////printpol(o2v(w[i]));
  }
  //exit(1);

  tmp.x[0] = 1;
  //

  //#pragma omp parallel for private(i,j)
  for (i = 0; i < T; i++)
  {
    t = v2o(tmp);
    //
    for (j = 0; j < T; j++)
    {
      // #pragma omp critical
      if (i != j)
      {
        t = omul(t, w[j]);
      }
    }

    //printpol(o2v(t));
    //
    if (deg(o2v(t)) == 0)
    {
      printf("baka9\n");
      // exit(1);
    }
    l = oadd(l, t);
  }

  return l;
}

//chen探索
vec chen(OP f)
{
  vec e = {0};
  int i, count = 0, n, x = 0;
  unsigned short z;

  n = odeg((f));
  //exit(1);
  //#pragma omp parallel for private(i)
  for (x = 0; x < N; x++)
  {
    z = 0;
    //#pragma omp parallel for reduction (^:z)
    for (i = 0; i < n + 1; i++)
    {
      if (f.t[i].a > 0)
        z ^= gf[mlt(mltn(f.t[i].n, fg[x]), fg[f.t[i].a])];
    }
    if (z == 0)
    {
      e.x[count] = x;
      count++;
      printf("x=%d\n", x);
    }
  }
  printpol(e);
  printf(" ==eee!\n");
  //exit(1);

  return e;
}

int oequ(OP f, OP g)
{
  vec v, x;
  int i, flg = 0;

  v = o2v(f);
  x = o2v(g);
  for (i = 0; i < DEG; i++)
  {
    if (v.x[i] != x.x[i])
      return -1;
  }

  return 0;
}

//GF(2^m) then set m in this function.
int ben_or(OP f)
{
  int i, n, flg = 0;
  OP s = {0}, u = {0}, r = {0};
  vec v = {0}, x = {0};
  //if GF(8192) is 2^m and m==13 or if GF(4096) and m==12 if GF(16384) is testing
  int m = E;
  // m=12 as a for GF(4096)=2^12 defined @ gloal.h or here,for example m=4 and GF(16)

  v.x[1] = 1;
  s = v2o(v);
  r = s;
  n = deg(o2v(f));

  if (n == 0)
    return -1;

  i = 0;

  //r(x)^{q^i} square pow mod
  while (i < n / 2 + 1)
  {
    printf("i=%d\n",i);
    flg = 1;
    // irreducible over GH(8192) 2^13
    //for(k=1;k<m+1;k++)
    r = opowmod1(r, f, N);

    // irreducible over GF2
    //r=omod(opow(r,2),f);

    u = oadd(r, s);
    if (deg(o2v(u)) == 0 && LT(u).a == 0)
      return -1;
    if (deg(o2v(u)) == 0 && LT(u).a == 1)
    {
      i++;
      flg = 0;
    }
    if (deg(o2v(u)) > 0)
      u = gcd(f, u);

    if (deg(o2v(u)) > 0)
      return -1;

    if (flg == 1)
      i++;
  }

  return 0;
}

//ユークリッドアルゴリズムによる復号関数
OP decode(OP f, OP s)
{
  int i, j, k, count = 0;
  OP r = {0}, w = {0}, e = {0}, l = {0};
  oterm t1, t2, d1, a, b;
  vec x = {0};
  unsigned short d = 0;
  OP h = {0};
  EX hh = {0};

  printf("in decode\n");
  printpol(o2v(s));
  printf("\nsyn===========\n");
  r = vx(f, s);
  //h=ogcd(f,s);

  if (odeg((r)) == 0)
  {
    printf("baka12\n");
    exit(1);
  }
  k = 0;
  // exit(1);
  x = chen(r);
  // exit(1);

  for (i = 0; i < T; i++)
  {
    printf("x[%d]=1\n", x.x[i]);
    if (x.x[i] == 0)
      k++;
    if (k > 1)
    {
      printf("baka0\n");
      printvec(o2v(f));
      //for (i = 0; i < N; i++)
      //printf("%d,", zz[i]);
      exit(1);
      //return f;
    }
  }
  //exit(1);

  //  printf("\n");

  printf("あっ、でる！\n");
  //  exit(1);

  if (odeg((r)) < T)
  {
    printpol(o2v(r));
    printf("baka5 deg(r)<T\n");
    exit(1);
  }

  w = bibun(x);
  //exit(1);
  //  w=oterml(w,d1);
  printpol(o2v(w));
  printf("@@@@@@@@@\n");
  //exit(1);

  //hh = xgcd(f, s);
  h = ogcd(f, s);
  //printpol(o2v(hh.d));
  printpol(o2v(h));
  //wait();

  //  exit(1);
  t1 = LT(r);

  t2.a = t1.a;
  t2.n = 0;

  if (odeg((w)) == 0)
  {
    printpol(o2v(w));
  }
  l = oterml(w, t2);

  j = deg(x) + 1;
  printf("%d\n", j);

  //    exit(1);

  for (i = 0; i < j; i++)
  {
    //if (x.x[i] > 0)
    {
      //e.t[i].a =
      //  gf[mlt(fg[trace(hh.d, x.x[i])], oinv(trace(l, x.x[i])))];
      e.t[i].a = gf[mlt(fg[trace(h, x.x[i])], oinv(trace(l, x.x[i])))];
      e.t[i].n = x.x[i];
    }
  }
  printpol(o2v(f));
  printf(" f============\n");
  printpol(o2v(l));
  printf(" l============\n");
  //  exit(1);

  for (i = 0; i < T; i++)
    if (gf[trace(h, x.x[i])] == 0)
      printf("h=0");
  //printf("\n");
  for (i = 0; i < T; i++)
    if (gf[oinv(trace(l, x.x[i]))] == 0)
      printf("l=0\n");
  //  printf("\n");

  return e;
}

//ユークリッドアルゴリズムによる復号関数
OP sendrier(OP f, OP s)
{
  int i, j, k, count = 0;
  OP r = {0}, w = {0}, e = {0}, l = {0};
  oterm t1, t2, d1, a, b;
  vec x = {0};
  unsigned short d = 0;
  OP h = {0};
  EX hh = {0};

  printf("in decode\n");
  printpol(o2v(s));
  printf("\nsyn===========\n");
  r = vx(f, s);
  //h=ogcd(f,s);

  if (odeg((r)) == 0)
  {
    printf("baka12\n");
    exit(1);
  }
  k = 0;
  // exit(1);
  x = chen(r);
  // exit(1);

  for (i = 0; i < T * 2; i++)
  {
    printf("x[%d]=1\n", x.x[i]);
    if (x.x[i] == 0)
      k++;
    if (k > 1)
    {
      printf("baka0\n");
      printvec(o2v(f));
      //for (i = 0; i < N; i++)
      //printf("%d,", zz[i]);
      exit(1);
      //return f;
    }
  }
  //exit(1);

  //  printf("\n");

  printf("あっ、でる！\n");
  //  exit(1);

  if (odeg((r)) < T * 2)
  {
    printpol(o2v(r));
    printf("baka5 deg(r)<T\n");
    exit(1);
  }

  w = bibun(x);
  //exit(1);
  //  w=oterml(w,d1);
  printpol(o2v(w));
  printf("@@@@@@@@@\n");
  //exit(1);

  //hh = xgcd(f, s);
  h = ogcd(f, s);
  //printpol(o2v(hh.d));
  printpol(o2v(h));
  //wait();

  //  exit(1);
  t1 = LT(r);

  t2.a = t1.a;
  t2.n = 0;

  if (odeg((w)) == 0)
  {
    printpol(o2v(w));
  }
  l = oterml(w, t2);

  j = deg(x) + 1;
  printf("%d\n", j);

  //    exit(1);

  for (i = 0; i < j; i++)
  {
    //if (x.x[i] > 0)
    {
      //e.t[i].a =
      //  gf[mlt(fg[trace(hh.d, x.x[i])], oinv(trace(l, x.x[i])))];
      e.t[i].a = gf[mlt(fg[trace(h, x.x[i])], oinv(trace(l, x.x[i])))];
      e.t[i].n = x.x[i];
    }
  }
  printpol(o2v(f));
  printf(" f============\n");
  printpol(o2v(l));
  printf(" l============\n");
  //  exit(1);

  for (i = 0; i < T; i++)
    if (gf[trace(h, x.x[i])] == 0)
      printf("h=0");
  //printf("\n");
  for (i = 0; i < T; i++)
    if (gf[oinv(trace(l, x.x[i]))] == 0)
      printf("l=0\n");
  //  printf("\n");

  return e;
}

//配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
  OP g;
  vec a;
  int i;

  memset(c, 0, sizeof(c));
  memcpy(c, f, 2 * n);
  a = Setvec(n);

  g = v2o(a);

  return g;
}

//バイナリ型パリティチェック行列を生成する
void bdet()
{
  int i, j, k, l;
  unsigned char dd[E * K] = {0};
  FILE *ff;

  //ff = fopen("Hb.key", "wb");

  for (i = 0; i < N; i++)
  {
    for (j = 0; j < K; j++)
    {
      l = mat[i][j];
      //#pragma omp parallel for
      for (k = 0; k < E; k++)
      {
        BB.x[i][j * E + k] = l % 2;
        l = (l >> 1);
      }
    }
  }

  for (i = 0; i < N; i++)
  {
    //#pragma omp parallel for
    for (j = 0; j < E * K; j++)
    {
      printf("%d,", BB.x[i][j]);
      //dd[j] = BH[j][i];
    }
    //fwrite(dd, 1, E * K, ff);
    printf("\n");
  }

  //fclose(ff);
}

//バイナリ型パリティチェック行列を生成する
void toBit(MTX L)
{
  int i, j, k, l;
  unsigned char dd[E * K] = {0};
  FILE *ff;

  //ff = fopen("Hb.key", "wb");

  for (i = 0; i < N; i++)
  {
    for (j = 0; j < K; j++)
    {
      l = L.x[i][j];
      printf("l=%d,", l);
      //#pragma omp parallel for
      for (k = 0; k < E; k++)
      {
        BB.x[i][j * E + k] = l % 2;
        l = (l >> 1);
      }
    }
    printf("\n");
  }
  //exit(1);

  for (i = 0; i < N; i++)
  {
    //#pragma omp parallel for
    for (j = 0; j < E * K; j++)
    {
      printf("%d,", BB.x[i][j]);
      //dd[j] = BH[j][i];
    }
    //fwrite(dd, 1, E * K, ff);
    printf("\n");
  }

  //fclose(ff);
}

unsigned short HH[N][K];

void toByte(MTX SH)
{
  vec v = {0};
  int i, j, k, cnt;

  for (i = 0; i < N; i++)
  {
    printf("%d\n",i);
    //#pragma omp parallel for
    for (j = 0; j < K; j++)
    {
      cnt = 0;
      for (k = j * E; k < j * E + E; k++)
        v.x[cnt++] = SH.x[i][k];

      HH[i][j] = v2i(v);
      
      //printf("%d,", HH[i][j]);
      //= BH[j][i];
    }
    //fwrite(dd, 1, E * K, ff);
    //printf("\n");
  }
  printf("end of byte\n");
  //exit(1);
}

//秘密置換を生成する
void Pgen()
{
  unsigned int i, j;
  FILE *fp;

  fp = fopen("P.key", "wb");
    for(i=0;i<N;i++)
    P[i]=i;
    random_shuffle(P,SIZE_OF_ARRAY(P));
//  random_permutation(P);
  for (i = 0; i < N; i++)
    inv_P[P[i]] = i;
  fwrite(P, 2, N, fp);
  fclose(fp);

  //for (i = 0; i < N; i++)
  //printf ("%d,", inv_P[i]);
  //printf ("\n");

  //fp = fopen("inv_P.key", "wb");
  //fwrite(inv_P, 2, N, fp);
  //fclose(fp);
}

//ハッシュ１６進表示
static void
byte_to_hex(uint8_t b, char s[23])
{
  unsigned i = 1;
  s[0] = s[1] = '0';
  s[2] = '\0';
  while (b)
  {
    unsigned t = b & 0x0f;
    if (t < 10)
    {
      s[i] = '0' + t;
    }
    else
    {
      s[i] = 'a' + t - 10;
    }
    i--;
    b >>= 4;
  }
}

//有限体の元の平方を計算する
int isqrt(unsigned short u)
{
  int i, j, k;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(i, i)] == u)
      return i;
  }

  printf("来ちゃいけないところに来ました\n");
  exit(1);
}

//多項式の平方を計算する
OP osqrt(OP f, OP w)
{
  int i, j, k, jj, n, flg = 0;
  OP even = {0}, odd = {0}, h = {0}, r = {0}, ww = {0}, s = {0}, tmp = {0}, t = {0};
  oterm o = {0};
  vec v = {0};

  j = 0;
  jj = 0;
  k = deg(o2v(f));
  for (i = 0; i < k + 1; i++)
  {
    if (f.t[i].n % 2 == 0 && f.t[i].a > 0)
    {
      even.t[j].n = f.t[i].n / 2;
      even.t[j++].a = gf[isqrt(f.t[i].a)];
      printf("a=%d %d\n", f.t[i].a, i);
    }
    if (f.t[i].n % 2 == 1 && f.t[i].a > 0)
    {
      odd.t[jj].n = (f.t[i].n - 1) / 2;
      odd.t[jj++].a = gf[isqrt(f.t[i].a)];
      printf(" odd %d\n", i);
      flg = 1;
    }
  }

  k = deg(o2v(w));
  //printf ("%d\n", k);
  //exit(1);
  j = 0;
  jj = 0;
  for (i = 0; i < k + 1; i++)
  {
    if (w.t[i].n % 2 == 0 && w.t[i].a > 0)
    {
      h.t[j].a = gf[isqrt(w.t[i].a)];
      h.t[j++].n = w.t[i].n / 2;
      printf("h==%d %d\n", (w.t[i].n / 2), i);
    }
    if (w.t[i].n % 2 == 1 && w.t[i].a > 0)
    {
      r.t[jj].a = gf[isqrt(w.t[i].a)];
      r.t[jj++].n = (w.t[i].n - 1) / 2;
      printf("r=====%d %d\n", (w.t[i].n - 1) / 2, i);
    }
  }
  printpol(o2v(r));
  printf(" sqrt(g1)=======\n");

  //  exit(1);
  if (LT(r).n > 0)
  {
    s = inv(r, w);
  }
  else if (LT(r).n == 0)
  {
    printpol(o2v(r));
    printf(" deg(r)======0!\n");
    printpol(o2v(w));
    printf(" goppa======0\n");
    printpol(o2v(f));
    printf(" syn======0\n");
    if (LT(r).a > 0)
    {
      s.t[0].a = gf[isqrt(LT(r).a)];
      s.t[0].n = 0;
      return omod(oadd(even, omul(coeff(h, s.t[0].a), odd)), w);
    }
    return even;
    //s = inv (r, w);
    //wait ();
    //exit (1);
  }
  if (deg(o2v(s)) > 0)
    tmp = omod(omul(s, r), w);
  if (odeg((tmp)) > 0)
  {
    //printpol (o2v (tmp));
    printf(" r is not inv==========\n");
    wait();
    exit(1);
  }
  if (LT(h).n > 0 && odeg(s) > 0)
  {
    ww = omod(omul(h, s), w);
  }
  if (LT(h).n == 0 || odeg(s) == 0)
  {
    printpol(o2v(h));
    printf(" h=========0\n");
    exit(1);
  }

  if (LT(ww).n == 0 && LT(ww).a == 0)
  {
    printpol(o2v(s));
    printf(" s===========\n");
    printsage(o2v(w));
    printf(" w==============\n");
    printpol(o2v(r));
    printf(" r===========\n");
    printpol(o2v(h));
    printf(" h============\n");
    printpol(o2v(ww));
    printf(" ww==============\n");
    printf(" wwが0になりました。error\n");
    wait();
    //return ww;;
    exit(1);
  }

  tmp = omod(omul(ww, ww), w);
  if (LT(tmp).n == 1)
  {
    printpol(o2v(ww));
    printf(" ww succsess!===========\n");
  }
  else
  {
    //printpol (o2v (tmp));
    printf(" mod w^2==========\n");
    //printpol (o2v (ww));
    printf(" ww^2 failed!========\n");
    printpol(o2v(s));
    printf(" g1^-1==============\n");
    //printpol (o2v (w));
    printf(" w==============\n");
    //printpol (o2v (h));
    printf(" g0===========\n");
    //printpol (o2v (r));
    printf(" r===========\n");
    printf("この鍵では逆元が計算できません。error");
    wait();
    //return ww;
    exit(1);
  }

  //    exit(1);
  printpol(o2v(s));
  printf(" g1^-1=========\n");
  printpol(o2v(h));
  printf(" g0=========\n");
  //exit(1);
  printpol(o2v(ww));
  printf(" ww==========\n");
  //  exit(1);
  h = ww;
  if (odeg(omod(omul(h, ww), w)) == 1)
  {
    ww = h;
    h = omod(oadd(even, omul(ww, odd)), w);
    return h;
  }
  else if (LT(ww).a == 0)
  {
    printf("vaka\n");
    exit(1);
  }

  // //printpol(o2v(ww));
  printf(" 来ちゃだめなところに来ました\n");

  exit(1);
}

vec p2()
{
}

EX extgcd(OP a, OP b)
{

  OP s = {0}, sx = {0}, sy = {0}, t = {0}, tx = {0}, ty = {0}, tmp = {0};
  EX c = {0};

  if (odeg(b) > odeg(a))
  {
    tmp = a;
    a = b;
    b = tmp;
  }
  s = a;
  t = b;
  sx.t[0].a = 1;
  sx.t[0].n = 0;
  ty.t[0].a = 1;
  ty.t[0].n = 0;

  //  OP temp={0};
  tmp = omod(s, t);
  if (odeg(tmp) == 0)
  {
    c.d = t;
    c.v = tx;
    c.u = ty;
    printf("ppp\n");
    return c;
  }
  while (odeg(tmp) > 0)
  {
    printpol(o2v((tmp)));
    printf(" ========omod\n");
    OP temp = odiv(s, t);
    OP u = oadd(s, omul(t, temp));
    OP ux = oadd(sx, omul(tx, temp));
    OP uy = oadd(sy, omul(ty, temp));
    /*
       */
    s = t;
    sx = tx;
    sy = ty;
    t = u;
    tx = ux;
    ty = uy;
    tmp = omod(s, t);
  }
  printpol(o2v(tmp));
  printf(" ========omod!\n");

  if (LT(tmp).a == 1)
  {
    c.d.t[0].a = 1;
    c.d.t[0].n = 0;
    //c.d=t;
    c.v = tx;
    c.u = ty;
    printf("bbb\n");
    return c;
  }
  if (LT(tmp).a == 0)
  {

    c.d = t;
    c.v = tx;
    c.u = ty;
    printf("ccc\n");

    return c;
  }
}

//パターソンアルゴリズムでバイナリGoppa符号を復号する
vec patterson(OP w, OP f)
{
  OP g1 = {0}, ll = {0}, op = {0};
  int i, j, k, l, c, n, count = 0, flg2 = 0;
  int flg, o1 = 0;
  OP tt = {0}, ff = {0}, h = {0};
  EX hh = {0}, opp = {0};
  vec v;
  oterm rr = {0};
  OP r2 = {0}, b2 = {
                   0};
  //unsigned short g[K+1]={2,2,12,1,2,8,4,13,5,10,8,2,15,10,7,3,5};

  tt.t[0].n = 1;
  tt.t[0].a = 1;
  o1 = 0;

  ff = inv(f, w);
  //  //printpol (o2v (ff));
  b2 = omod(omul(ff, f), w);
  if (odeg((b2)) > 0)
  {
    printf("逆元が計算できません。error\n");
    //wait ();
    exit(1);
  }
  printf("locater==========\n");
  //exit(1);
  r2 = oadd(ff, tt);
  printpol(o2v(r2));
  printf(" h+x==============\n");
  //wait();
  //  exit(1);
  g1 = osqrt(r2, w);
  printpol(o2v(g1));
  printf(" sqrt(h+x)=======\n");
  b2 = omod(omul(g1, g1), w);
  printpol(o2v(b2));
  printf(" g1*g1%w===========\n");
  if (deg(o2v(b2)) == 0)
  {
    printf("b2 is 0\n");
    exit(1);
  }
  if (LT2(b2).a != LT2(r2).a)
  {
    printpol(o2v(w));
    printf(" w============\n");
    printpol(o2v(r2));
    printf(" r2============\n");
    printpol(o2v(g1));
    printf(" g1============\n");
    printf(" g1は平方ではありません。error");
    //wait ();
    exit(1);
  }
  printpol(o2v(w));
  printf(" ========goppa\n");
  printpol(o2v(g1));
  printf(" g1!=========\n");
  if (LT(g1).n == 0 && LT(g1).a == 0)
  {
    //printpol (o2v (w));
    printf(" badkey=========\n");
    //printpol (o2v (g1));
    printf(" g1============\n");
    printf("平方が０になりました。 error\n");
    //wait ();
    exit(1);
  }
  h = gcd(f, w);
  if (odeg((h)) > 0)
  {
    printsage(o2v(w));
    printf(" =goppa\n");
    printsage(o2v(f));
    printf(" =syn\n");
    printpol(o2v(h));
    printf(" =gcd\n");

    //printf("f-ben_or=%d\n",ben_or(f));
    //printf("w-ben_or=%d\n",ben_or(w));
    printf(" s,wは互いに素じゃありません。\n");
    flg2 = 1;
    //wait ();
    exit(1);
    //goto label;
  }

  //exit(1);
  OP ppo = {0};

  hh = xgcd(w, g1);
  ppo = zgcd(w, g1);
  printpol(o2v(ppo));
  printf(" =========ppo\n");
  printpol(o2v(ppo));
  printf(" =========ppo\n");
  //if(LT(ppo).a!=LT(ppo).a)
  //exit(1);
  flg = 0;
  printpol(o2v(ppo));
  printf(" =========ppo\n");
  ff = omod(omul(ppo, g1), w);
  printpol(o2v(ppo));
  printf(" alpha!=========\n");
  printpol(o2v(ff));
  printf(" beta!=========\n");
  printpol(o2v(w));
  printf(" goppa=========\n");
  printpol(o2v(f));
  printf(" syn=========\n");

  printpol(o2v(ppo));
  printf(" alpha!=========\n");
  hh = xgcd(w, g1);
  ppo = zgcd(w, g1);
  //  h=zgcd(w,g1);
  ff = omod(omul(ppo, g1), w);
  ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo)));
  v = chen(ll);
  if (v.x[K - 1] > 0)
  {
    //wait ();
    return v;
  }

  ppo = zgcd(w, g1);
  ff = omod(omul(ppo, g1), w);
  if (odeg((ff)) != K / 2)
  {

    printf("\nbefore h.d\n");
    ff = omod(omul(ppo, g1), w);
    flg = 1;
    printpol(o2v(ff));
    printf(" ==========beta!\n");
    printpol(o2v(ppo));
    printf(" alpha!=========\n");
    ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo)));
    v = chen(ll);
    if (v.x[K - 1] > 0)
    {
      wait();
      return v;
    }
  }

  if (odeg((ff)) == 1)
  {
    ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo))); //ff;
    printf("deg==1\n");
    exit(1);
  }

  printf("あっ、でる・・・！\n");
  count = 0;
  printpol(o2v(ll));
  printf(" ll=================\n");
  printpol(o2v(f));
  printf(" syn=================\n");
  v = chen(ll);

  if (v.x[K - 1] > 0)
  {
    return v;
  }

  printf("flg2=%d\n", flg2);
  //exit(1);
  //ero2(v);
  return v;
}

//512bitの秘密鍵を暗号化
void encrypt(char buf[], unsigned char sk[64])
{
  const uint8_t *hash = {0};
  sha3_context c = {0};
  int image_size = 512, i;
  FILE *fp;
  //  unsigned short dd=0;

  printf("plain text=");
  for (i = 0; i < 64; i++)
    printf("%u,", sk[i]);
  printf("\n");
  //  puts(buf);
  //printf("\n");
  //exit(1);

  //scanf("%s",buf);
  sha3_Init256(&c);
  sha3_Update(&c, (char *)buf, strlen(buf));
  hash = sha3_Finalize(&c);

  //j=0;

  for (i = 0; i < 64; i++)
  {
    printf("%d", hash[i]);
    //char s[3];
    //byte_to_hex(hash[i],s);

    sk[i] ^= hash[i];
  }
  printf("\nencrypt sk=");
  for (i = 0; i < 64; i++)
    printf("%d,", sk[i]);
  printf("\n");

  fp = fopen("enc.sk", "wb");
  fwrite(sy, 2, K, fp);
  fwrite(sk, 1, 64, fp);
  fclose(fp);
}

void decrypt(OP w)
{
  FILE *fp;
  int i, j;
  unsigned char sk[64] = {0}, err[N] = {
                                  0};
  unsigned short buf[K] = {0}, tmp[K] = {
                                   0};
  OP f = {0}, r = {0};
  vec v = {0};
  const uint8_t *hash = {0};
  sha3_context c = {0};
  int image_size = 512;

  j = 0;
  fp = fopen("enc.sk", "rb");

  fread(tmp, 2, K, fp);
  fread(sk, 1, 64, fp);
  fclose(fp);

  for (i = 0; i < K; i++)
    buf[i] = tmp[K - i - 1];

  printf("in decrypt\n");
  // reconstruct syndrome polynomial
  f = setpol(buf, K);
  // decode error from syndrome
  r = decode(w, f);
  elo2(r);
  //exit(1);

  // elo(r);
  //exit(1);
  //v=o2v(r);

  j = 0;
  if (r.t[1].a > 0 && r.t[0].a == 0)
  {
    err[0] = 1;
    j++;
  }

  printf("j=%d\n", j);
  printf("after j\n");
  for (i = j; i < T; i++)
  {
    if (r.t[i].a > 0 && r.t[i].n < N)
    {
      err[r.t[i].n] = 1;
    }
  }

  char buf0[8192] = {0}, buf1[10] = {
                             0};

  //make error pattarn e for Hash(e)
  for (i = 0; i < N; i++)
  {
    snprintf(buf1, 10, "%d", err[i]);
    strcat(buf0, buf1);
  }
  //puts (buf0);
  printf("vector=%d\n", strlen(buf0));
  //exit(1);
  printf("cipher sk2=");
  for (i = 0; i < 64; i++)
    printf("%u,", sk[i]);
  printf("\n");

  sha3_Init256(&c);
  sha3_Update(&c, (char *)buf0, strlen(buf0));
  hash = sha3_Finalize(&c);

  j = 0;
  printf("hash=");
  for (i = 0; i < 64; i++)
  {
    printf("%d", hash[i]);
    //char s[3];
    //byte_to_hex(hash[i],s);

    // m=c^Hash(e)
    sk[i] ^= hash[i];
  }
  printf("\ndecript sk=");
  for (i = 0; i < 64; i++)
    printf("%u,", sk[i]);
  printf("\n");
  //  exit(1);

  return;
}

OP synd(unsigned short zz[])
{
  unsigned short syn[K] = {0};
  unsigned short s = 0;
  int i, j, t1;
  OP f = {0};

  printf("in synd2\n");

  for (i = 0; i < K; i++)
  {
    syn[i] = 0;
    s = 0;
    //#pragma omp parallel num_threads(16)
    for (j = 0; j < N; j++)
    {
      s ^= gf[mlt(fg[zz[j]], fg[mat[j][i]])];
    }
    syn[i] = s;
    //printf ("syn%d,", syn[i]);
  }
  //printf ("\n");

  f = setpol(syn, K);
  printpol(o2v(f));
  printf(" syn=============\n");
  //  exit(1);

  return f;
}

//64バイト秘密鍵の暗号化と復号のテスト
void test(OP w, unsigned short zz[])
{
  int i;
  vec v = {0};
  const uint8_t *hash;
  sha3_context c;
  //int image_size=512;
  OP f = {0}, r = {0};
  FILE *fp;

  fp = fopen("aes.key", "rb");
  /*
     static char base64[] = {
     'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
     'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
     'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
     'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
     'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
     'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
     'w', 'x', 'y', 'z', '0', '1', '2', '3',
     '4', '5', '6', '7', '8', '9', '+', '/',
     };
   */
  char buf[8192] = {0}, buf1[10] = {
                            0};
  unsigned char sk[64] = {0};
  // unsigned short s[K]={0};
  //fread(sk,1,32,fp);
  for (i = 0; i < 64; i++)
    sk[i] = i + 1;

  for (i = 0; i < N; i++)
  {
    snprintf(buf1, 10, "%d", zz[i]);
    strcat(buf, buf1);
  }
  //puts (buf);
  printf("vector=%u\n", strlen(buf));
  //exit(1);

  printf("sk0=");
  for (i = 0; i < 64; i++)
    printf("%u,", sk[i]);
  printf("\n");
  //exit(1);

  f = synd(zz);
  //r=decode(w,f);
  //elo2(r);
  //exit(1);

  v = o2v(f);
  //printf("v=");

  for (i = 0; i < odeg(f) + 1; i++)
  {
    sy[i] = v.x[i];
    printf("%d,", sy[i]);
  }
  printf("\n");

  //buf=e, sk=m, Hash(e)^m
  encrypt(buf, sk);
  decrypt(w);

  sha3_Init256(&c);
  sha3_Update(&c, (char *)buf, strlen(buf));
  hash = sha3_Finalize(&c);
}

void readkey()
{
  FILE *fp, *fq;
  unsigned short dd[K * N] = {0};
  int i, j;
  unsigned char r[K * E] = {0};
  vec v = {0};
  unsigned short o[K] = {0};

  //鍵をファイルに書き込むためにはkey2を有効にしてください。

  fp = fopen("sk.key", "rb");
  fread(g, 2, K + 1, fp);
  fclose(fp);
  //固定した鍵を使いたい場合はファイルから読み込むようにしてください。
  fq = fopen("Pub.key", "rb");
  fread(dd, 2, K * N, fq);
  //#pragma omp parallel for
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < K; j++)
    {
      HH[i][j] = dd[K * i + j];
      //printf("%d,", HH[i][j]);
    }
    //printf(" ===HH\n");
  }
  printf("\n");
  fclose(fq);
  //exit(1);

  fq = fopen("P.key", "rb");
  fread(P, 2, N, fp);
  fclose(fq);

  fq = fopen("inv_S.key", "rb");
  /*
for(i=0;i<K*E;i++){
  fread(r,1,K*E,fq);
  for(j=0;j<K*E;j++)
  inv_S.x[i][j]=r[j];
}
  fclose(fq);
*/
  memset(inv_S.x, 0, sizeof(inv_S.x));
  for (i = 0; i < K * E; i++)
  {
    fread(o, 2, K, fq);
    for (j = 0; j < K; j++)
    {
      v = i2v(o[j]);
      //printf("%d,", o[j]);
      for (int k = 0; k < E; k++)
      {
        inv_S.x[i][j * E + k] = v.x[k];
        //printf("%d,",inv_S.x[i][j]);
      }
    }
    //printf("\n");
  }
  fclose(fq);
  //exit(1);

  for (i = 0; i < K * E; i++)
  {
    for (j = 0; j < K * E; j++)
      printf("%d,", inv_S.x[i][j]);
    printf("\n");
  }
  //exit(1);
}

//OP sx={0},ty={0};

unsigned short vb[K * 2][N] = {0};
unsigned short gt[K * 2][K * 2] = {0};

void van()
{
  int i, j, k;

  printf("van der\n");

  for (i = 0; i < N; i++)
    vb[0][i] = 1;
  //#pragma omp parallel for private(i, j)
  for (i = 1; i < K; i++)
  {
    for (j = 0; j < N; j++)
      vb[i][j] = gf[mltn(i, fg[j])];
  }
}

void ogt()
{
  int i, j, k;
  OP w = {0};
  unsigned short abc[N][K] = {0};

#pragma omp parallel for private(i, j)
  for (i = 0; i < K; i++)
  {
    for (j = 0; j < K - i; j++)
      gt[i][j + i] = g[j];
  }
}
void van2()
{
  int i, j, k;

  printf("van der\n");

  for (i = 0; i < N; i++)
    vb[0][i] = 1;
  //#pragma omp parallel for private(i, j)
  for (i = 1; i < K * 2; i++)
  {
    for (j = 0; j < N; j++)
      vb[i][j] = gf[mltn(i, fg[j])];
  }
}

void ogt2()
{
  int i, j, k;
  OP w = {0};
  unsigned short abc[N][K] = {0};

#pragma omp parallel for private(i, j)
  for (i = 0; i < K * 2; i++)
  {
    for (j = 0; j < K * 2 - i; j++)
      gt[i][j + i] = g[j];
  }
}

int isquad(OP w)
{
  int i, j, flg = 0;
  vec b = {0};

  b = o2v(w);
  for (i = 0; i < DEG; i++)
  {
    if (b.x[i] > 0 && i % 2 == 1)
      return 0;
  }

  return -1;
}

OP mkpol()
{
  int i, j, k, fail, flg, l, ii = 0;
  OP w = {0};

  do
  {
    fail = 0;
    j = 0;
    k = 0;
    flg = 0;
    l = 0;
    memset(g, 0, sizeof(g));
    //memset(ta, 0, sizeof(ta));
    memset(w.t, 0, sizeof(w));
    ginit();
    ii++;
    if (ii > 300)
    {
      printf("erro=%d\n", ii);
      exit(1);
    }

    for (i = 0; i < K; i++)
    {
      if (g[K - 1] > 0)
        flg = 1;
      if (i % 2 == 1 && g[i] > 0 && i < K)
        k++;
    }

    //偶数項だけにならないようにする
    if ((k > 0 && flg == 0) || (k > 1 && flg == 1))
    //if(k>0)
    {
      w = setpol(g, K + 1);
      j = 1;
      //if(isquad(w)==-1)
      //exit(1);
    }
    // exit(1);

  } while (j == 0);

  return w;
}

unsigned short dd[N][N] = {0};

OP mkg()
{
  int i, j, k, l, ii = 0;
  OP w = {0};
  unsigned short tr[N] = {0};
  unsigned short ta[N] = {0};

aa:

  //printf("\n");

  //既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
  //既約多項式しか使わない。
  l = -1;
  ii = 0;
  // irreducible goppa code
  /*
  while (l == -1)
  {
    w = mkpol();
    l = ben_or(w);
    printf("irr=%d\n", l);
    if (ii > 300)
    {
      printf("too many tryal\n");
      exit(1);
    }
    ii++;
    //
  }
*/

  //separable goppa code
  w = mkpol();
  //w=setpol(g,K+1);
  //多項式の値が0でないことを確認
  for (i = 0; i < N; i++)
  {
    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      //fail = 1;
      goto aa;
    }
  }
  for (i = 0; i < N; i++)
  {
    tr[i] = oinv(ta[i]);
    //printf("%d,", tr[i]);
  }

  printpol(o2v(w));
  printf(" =irreducible\n");
  printsage(o2v(w));
  printf("\n");
  //wait();

  //多項式を固定したい場合コメントアウトする。
  /*
  memset(ta, 0, sizeof(ta));
  w = setpol(g, K + 1);
  printpol(o2v(w));
  //printf(" =poly\n");
  for (i = 0; i < N; i++)
  {    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      fail = 1;
      break;
    }
  }
  oprintpol(w);
  printf("\n");
  printsage(o2v(w));
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");
  */

  van();
  ogt();
  memset(mat, 0, sizeof(mat));

  //wait();

  //#pragma omp parallel for

  printf("\nすげ、オレもうイキそ・・・\n");
  //keygen(g);
  //exit(1);

  for (j = 0; j < N; j++)
  {
    for (i = 0; i < K; i++)
    {
      ma[j][i] = gf[mlt(fg[vb[i][j]], tr[j])];
    }
    //printf("tr[%d]=%d\n",j,tr[j]);
  }

  unsigned short s;
#pragma omp parallel for default(none) private(i, j, k, s) shared(mat, gt, ma, gf, fg)
  for (i = 0; i < K; i++)
  {
    for (j = 0; j < N; j++)
    {
      s = 0;

      for (k = 0; k < K; k++)
        s ^= gf[mlt(fg[gt[k][i]], fg[ma[j][k]])];
      //printf("%d,",s);
      mat[j][i] = s;
    }
  }
  //printf("\n");
  //exit(1);

  for (j = 0; j < N; j++)
  {
    for (i = 0; i < K; i++)
      printf("%d,", mat[j][i]);
    printf("\n");
  }

  //wait();

  return w;
}

OP mkd(OP w)
{
  int i, j, k, l, ii = 0;

  unsigned short tr[N] = {0};
  unsigned short ta[N] = {0};
  vec v = {0};

aa:

  //printf("\n");

  //既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
  //既約多項式しか使わない。
  l = -1;
  ii = 0;

  memset(ta, 0, sizeof(ta));
  //w = setpol(g, K + 1);
  printpol(o2v(w));
  //printf(" =poly\n");

  //多項式の値が0でないことを確認
  for (i = 0; i < N; i++)
  {
    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      //fail = 1;
      exit(1);
    }
  }
  for (i = 0; i < N; i++)
  {
    tr[i] = oinv(ta[i]);
    //printf("%d,", tr[i]);
  }

  //多項式を固定したい場合コメントアウトする。
  oprintpol(w);
  printf("\n");
  printsage(o2v(w));
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");

  van2();
  ogt2();
  memset(bm, 0, sizeof(bm));

  //wait();

  //#pragma omp parallel for

  printf("\nすげ、オレもうイキそ・・・\n");
  //keygen(g);
  //exit(1);

  for (j = 0; j < N; j++)
  {
    for (i = 0; i < K * 2; i++)
    {
      bm[j][i] = gf[mlt(fg[vb[i][j]], tr[j])];
    }
    //printf("tr[%d]=%d\n",j,tr[j]);
  }

  unsigned short s;
#pragma omp parallel for default(none) private(i, j, k, s) shared(bm2, gt, bm, gf, fg)
  for (i = 0; i < K; i++)
  {
    for (j = 0; j < N; j++)
    {
      s = 0;

      for (k = 0; k < K; k++)
        s ^= gf[mlt(fg[gt[k][i]], fg[bm[j][k]])];
      //printf("%d,",s);
      bm2[j][i] = s;
    }
  }
  //printf("\n");
  //exit(1);

  for (j = 0; j < N; j++)
  {
    for (i = 0; i < K * 2; i++)
      printf("%d,", bm2[j][i]);
    printf("\n");
  }
  //exit(1);
  //wait();

  return w;
}

//Niederreiter暗号の公開鍵を作る
OP pubkeygen()
{
  int i, j, k, l;
  unsigned short n[K] = {0};
  FILE *fp;
  unsigned short dd[K] = {0};
  OP w = {0};
  vec v = {0};
  MTX Q={0},O;

  w = mkg();
  printpol(o2v(w));
  printf(" ==goppa polynomial\n");

  v = o2v(w);
  fp = fopen("sk.key", "wb");
  fwrite(g, 2, K + 1, fp);
  fclose(fp);
  oprintpol(w);
  printf("\n");
  printsage(o2v(w));
  printf("\n");
  printf("sagemath で既約性を検査してください！\n");

  bdet();
  //  toByte(BB);
  //exit(1);

  Pgen();
  fp = fopen("P.key", "w");
  //  fwrite(P, 2, N, fp);
  fclose(fp);
  //makeS();
  do{
  memset(Q.x,0,sizeof(Q.x));
  memset(O.x,0,sizeof(O.x));
  memset(S.x,0,sizeof(S.x));
for(i=0;i<(K)*E;i++){
  for(j=0;j<(K)*E;j++)
  S.x[i][j]=xor128()%2;
}
}while(is_reg(S,inv_S.x) == -1);

  fp = fopen("inv_S.key", "wb");

  /*
for(i=0;i<K*E;i++){
  for(j=0;j<K*E;j++)
    n[j]=inv_S.x[i][j];  
    fwrite(n,1,K*E,fp);
} 
*/
  fclose(fp);

  for (i = 0; i < K * E; i++)
  {
    for (j = 0; j < K; j++)
    {
      memset(v.x, 0, sizeof(v.x));
      for (k = 0; k < E; k++)
        v.x[k] = inv_S.x[i][j * E + k];
      n[j] = v2i(v);
      //printf("%d,", n[j]);
    }
    //printf("\n");
    //  fwrite(n, 2, K, fp);
  }

  // SH
  H = mulmat(S, BB, 1);
  // SHP
  for (i = 0; i < K * E; i++)
  {
    for (j = 0; j < N; j++)
      S.x[j][i] = H.x[P[j]][i];
  }
  //binary を unsigned short にパッキング
  toByte(S);
  fp = fopen("Pub.key", "w");
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < K; j++)
      dd[j] = HH[i][j];
    //fwrite(dd, 2, K, fp);
  }
  fclose(fp);

  return w;
}

OP dec(unsigned short ss[])
{
  int i, j, k;
  vec v = {0};
  OP s = {0};
  unsigned ch[K * E] = {0};
  unsigned char h2o[K * E] = {0};

  for (i = 0; i < K; i++)
  {
    v = i2v(ss[i]);
    for (j = 0; j < E; j++)
      ch[i * E + j] = v.x[j];
  }
  for (i = 0; i < K * E; i++)
    printf("%d", ch[i]);
  printf("\n");

  unsigned short uk[K] = {0};

  for (i = 0; i < K * E; i++)
  {
    for (j = 0; j < K * E; j++)
      h2o[i] ^= (ch[j] & inv_S.x[i][j]);
  }
  //for (i = 0; i < K * E; i++)
  //printf("%d,", h2o[i]);
  //printf("\n");

  for (i = 0; i < K; i++)
  {
    memset(v.x, 0, sizeof(v.x));
    for (j = 0; j < E; j++)
      v.x[j] = h2o[i * E + j];
    uk[i] = v2i(v);
  }
  for (i = 0; i < K; i++)
    printf("%d,", uk[i]);
  printf("\n");
  //    exit(1);
  s = setpol(uk, K);

  return s;
}

//鍵生成
void key2(unsigned short g[])
{
  FILE *fp;
  unsigned short dd[K] = {0};
  int i, j, k;

  printf("鍵を生成中です。４分程かかります。\n");
  fp = fopen("H.key", "wb");
  i = 0;

  mkg();

  //exit(1);
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < K; j++)
      dd[j] = mat[i][j];
    fwrite(dd, 2, K, fp);
  }
  fclose(fp);
  fp = fopen("sk.key", "wb");
  fwrite(g, 2, K + 1, fp);
  fclose(fp);
}

int elo(OP r)
{
  int count, i, j, k;

  count = 0;

  for (i = 0; i < T; i++)
  {

    if (i > 0 && r.t[i].n == 0)
    {
      printf("err baka-z\n");
      //return -1;
      exit(1);
    }

    if (r.t[i].a > 0 && i > 0) // == r.t[i].n)
    {
      printf("err=%d %d %s\n", r.t[i].a, r.t[i].n, "お");
      count++;
    }
    if (i == 0)
    {
      printf("\nerr=%d %d %s\n", r.t[i].a, r.t[i].n, "う");
      count++;
    }
    //zz[r.t[i].n]=r.t[i].a;
  }
  if (count != T)
  {
    printf("error pattarn too few %d\n", count);
    return -1;
  }
  //exit(1);

  return count;
}

int elo2(OP r)
{
  int count, i, j, k;

  count = 0;

  unsigned char x[N] = {0}, yy[N] = {0};
  for (i = 0; i < T; i++)
  {

    if (i > 0 && r.t[i].n == 0)
    {
      printf("err baka-z\n");
      //return -1;
      exit(1);
    }

    if (r.t[i].a > 0 && i > 0) // == r.t[i].n)
    {
      x[r.t[i].n] = r.t[i].a;
      //printf("err=%d %d %s\n", r.t[i].a, r.t[i].n, "お");
      count++;
    }
    if (i == 0)
    {
      x[r.t[i].n] = r.t[i].a;
      //printf("\nerr=%d %d %s\n", r.t[i].a, r.t[i].n, "う");
      count++;
    }
    //zz[r.t[i].n]=r.t[i].a;
  }

  for (i = 0; i < N; i++)
    yy[i] = x[P[i]];
  for (i = 0; i < N; i++)
  {
    if (yy[i] > 0)
      printf("err= %d\n", i);
  }

  return count;
}

int ero(vec v)
{
  int i, j, count = 0;

  for (i = 0; i < T * 2; i++)
  {
    if (i == 0)
    {
      printf("error position=%d %d う\n", i, v.x[i]);
      count++;
    }
    if (i > 0 && v.x[i] > 0)
    {
      printf("error position=%d %d お\n", i, v.x[i]);
      count++;
    }
    if (i > 0 && v.x[i] == 0)
    {
      printf("baka %d %d\n", i, v.x[i]);
      printf("v.x[K-1]=%d\n", v.x[K - 1]);
      /*
        printpol(o2v(w));
        printf(" ============goppa\n");
        printsage(o2v(w));
        printf(" ============sage\n");
        printsage(o2v(f));
        printf(" ============syn\n");
        printpol(o2v(f));
        printf(" ==========synd\n");
        printf("{");
        for (k = 0; k < N; k++)
        {
          if (z1[k] > 0)
            printf("%d,", z1[k]);
        }
        printf("};\n");
        //AA++;
        //wait();
        */
      break;
      //
      //exit (1);
    }
    int cnt = 0;
    /*
      for (k = 0; k < N; k++)
      {
        if (z1[k] > 0)
        {
          if (k != v.x[cnt])
          {
            printf("%d,%d\n", k, v.x[cnt]);
            printsage(o2v(w));
            printf(" ========w\n");
            AA++;
            break;
            //exit(1);
          }
          cnt++;
        }
      }
*/
  }

  if (count == T * 2)
  {
    printf("err=%dっ!! \n", count);
    B++;
  }
  if (count < T * 2)
  {
    printf("error is too few\n");

    AA++;
    //memcpy (zz, z1, sizeof (zz));
    /*
      printf("{");
      for (i = 0; i < D; i++)
        printf("%d,", z1[i]);
      printf("};\n");
      printpol(o2v(w));
      printf(" =========goppa\n");
      printsage(o2v(w));
      printf(" =========sage\n");
      printsage(o2v(f));
      printf(" =========syn\n");
      printpol(o2v(f));
      printf(" ==========synd\n");
      */
    printf("へげえええーっ\n");
    //exit(1);
    exit(1);
  }

  return count;
}

int ero2(vec v)
{
  int i, j, count = 0;
  unsigned short ya[N] = {0}, xa[N] = {0};

  for (i = 0; i < T * 2; i++)
  {
    if (i == 0)
    {
      xa[v.x[i]] = 1;
      //printf("error position=%d %d う\n", i, v.x[i]);
      count++;
    }
    if (i > 0 && v.x[i] > 0)
    {
      xa[v.x[i]] = 1;
      //printf("error position=%d %d お\n", i, v.x[i]);
      count++;
    }
    if (i > 0 && v.x[i] == 0)
    {
      printf("baka %d %d\n", i, v.x[i]);
      printf("v.x[K-1]=%d\n", v.x[K - 1]);
      /*
        printpol(o2v(w));
        printf(" ============goppa\n");
        printsage(o2v(w));
        printf(" ============sage\n");
        printsage(o2v(f));
        printf(" ============syn\n");
        printpol(o2v(f));
        printf(" ==========synd\n");
        printf("{");
        for (k = 0; k < N; k++)
        {
          if (z1[k] > 0)
            printf("%d,", z1[k]);
        }
        printf("};\n");
        //AA++;
        //wait();
        */
      break;
      //
      //exit (1);
    }
    int cnt = 0;
    /*
      for (k = 0; k < N; k++)
      {
        if (z1[k] > 0)
        {
          if (k != v.x[cnt])
          {
            printf("%d,%d\n", k, v.x[cnt]);
            printsage(o2v(w));
            printf(" ========w\n");
            AA++;
            break;
            //exit(1);
          }
          cnt++;
        }
      }
*/
  }
  for (i = 0; i < N; i++)
    ya[i] = xa[P[i]];
  for (i = 0; i < N; i++)
  {
    if (ya[i] > 0 && i == 0)
    {
      printf("error position=%d う\n", i);
    }
    else if (ya[i] > 0)
    {
      printf("error position=%d お\n", i);
    }
  }
  //exit(1);

  if (count == T * 2)
  {
    printf("err=%dっ!! \n", count);
    B++;
  }
  if (count < T * 2)
  {
    printf("error is too few\n");

    AA++;
    //memcpy (zz, z1, sizeof (zz));
    /*
      printf("{");
      for (i = 0; i < D; i++)
        printf("%d,", z1[i]);
      printf("};\n");
      printpol(o2v(w));
      printf(" =========goppa\n");
      printsage(o2v(w));
      printf(" =========sage\n");
      printsage(o2v(f));
      printf(" =========syn\n");
      printpol(o2v(f));
      printf(" ==========synd\n");
      */
    printf("へげえええーっ\n");
    //exit(1);
    exit(1);
  }

  return count;
}

void mkerr(unsigned short *z1, int num)
{
  int i, j, l;
  j = 0;

  while (j < num)
  {
    l = xor128() % N;
    //printf ("l=%d\n", l);
    if (0 == z1[l])
    {
      z1[l] = 1;
      printf("l=%d %d\n", l, j);
      j++;
    }
  }

  for (i = 0; i < N; i++)
  {
    if (z1[i] > 0)
      printf("%d=%d\n", i, z1[i]);
  }
}

void fun()
{
  unsigned short i, k;

  OP qq = {0};
  for (i = 0b1000000000001; i < 0b1111111111111 + 1; i++)
  {
    qq = v2o(i2v(i));
    k = ben_or(qq);
    if (k == 0)
    {
      printpol(o2v(qq));
      printf(" =irreducible\n");
    }
  }
}

vec sin2(unsigned short zz[])
{
  int i, j;
  OP s = {0};
  vec v = {0};
  unsigned short ss[K] = {0};

  for (i = 0; i < N; i++)
  {
    if (zz[i] > 0)
    {
      for (j = 0; j < K; j++)
      {
        v.x[j] ^= HH[i][j];
        //printf("%d,", HH[i][j]);
      }
      //printf("\n");
    }
  }

  return v;
}

OP cos2(unsigned short zz[])
{
  unsigned short ss[K * 2] = {0};
  int i, j;
  OP s = {0};
  vec v = {0};

  //unsigned short ss[K] = {0};

  for (i = 0; i < N; i++)
  {
    if (zz[i] > 0)
    {
      for (j = 0; j < K * 2; j++)
      {
        ss[j] ^= gf[mlt(fg[zz[i]], fg[bm2[i][j]])];
        //printf("%d,", HH[i][j]);
      }
      //printf("\n");
    }
    printf("ss==%d\n", ss[i]);
  }

  for (j = 0; j < K; j++)
    printf("%d,", ss[j]);
  printf(" ==ss\n");
  //exit(1);
  s = setpol(ss, K * 2);

  return s;
}

OP kof(unsigned short c, OP f)
{
  int i, j, k;
  vec b = {0}, h = {0};
  OP g = {0};

  b = o2v(f);
  k = deg(b);
  for (i = 0; i < k + 1; i++)
  {
    h.x[i] = gf[mlt(fg[c], fg[b.x[i]])];
  }
  g = v2o(h);

  return g;
}

unsigned short logx(unsigned short u)
{
  unsigned short i;

  return oinv(u);

  printf("baka-von\n");
}

OP rev(OP f)
{
  int i, tmp, j = 0, c[512] = {0}, d[512] = {0}, count = 0;
  vec v = {0};

  j = odeg(f) + 1;
  printf("d=");
  for (i = 0; i < j; i++)
  {
    d[count] = f.t[i].n;
    c[count] = f.t[i].a;
    printf("%d,", d[count]);
    count++;
  }
  printf("\n");
  printf("c=");
  for (i = 0; i < count; i++)
    printf("%d,", c[i]);
  printf("\n");
  for (i = 0; i < count; i++)
    v.x[d[count - i - 1]] = c[i];

  printpol(v);
  printf(" ==rev?\n");
  //exit(1);
  f = v2o(v);

  return f;
}

vec kof2(unsigned short c, vec f)
{
    int i, j, k;
    vec b = {0}, h = {0};
    OP g = {0};

    c = fg[c];
    printf("c=%d\n", c);
    // exit(1);
    b = f; // o2v(f);
    k = deg(b);
    printpol(b);
    printf(" =b debugi\n");
    for (i = 0; i < k + 1; i++)
    {
        h.x[i] = gf[mlt(c, fg[b.x[i]])];
    }
    // g = v2o(h);
    printpol(h);
    printf(" =h in kof2\n");
    return h;
}

vec vadd2(vec a, vec b)
{
    int i;
    vec c = {0};

    // printf("deg=%d %d\n",deg(a),deg(b));

    for (i = 0; i < DEG; i++)
        c.x[i] = a.x[i] ^ b.x[i];

    return c;
}

int mul = 0, mul2 = 0;
vec vmul2(vec a, vec b)
{
    int i, j, k, l;
    vec c = {0};
    if (deg(a) > 128 && deg(b) > 128)
        mul++;
    mul2++;

    k = deg(a);
    l = deg(b);

    for (i = 0; i < k + 1; i++)
    {
        for (j = 0; j < l + 1; j++)
        // if (a.x[i] > 0)
        {
            c.x[i + j] ^= gf[mlt(fg[a.x[i]], fg[b.x[j]])];
        }
    }

    return c;
}

vec bms(unsigned short s[])
{
    int L = 0, m = -1, d[K] = {0}, k = 0, i, e;
    vec f = {0}, g = {0}, h, v;

    f.x[0] = g.x[0] = 1;

    while (k <= (2 * T - 1))
    {
        e = 0;
        for (i = 0; i < L; i++)
            e ^= gf[mlt(fg[f.x[i]], fg[s[k - i]])];
            
        d[k] = gf[mlt(fg[f.x[i]], fg[s[k - i]])] ^ e; // s[k] ^ e;
        if (d[k] > 0)
        {
            h = f;
            memset(v.x, 0, sizeof(v.x));
            v.x[k - m] = 1;

            unsigned short a;
            a = (m < 0) ? 1 : oinv(d[m]);
            f = vadd2(f, vmul2(kof2(gf[mlt(fg[d[k]], a)], g), v));
            if (L <= k / 2)
            {
                L = k + 1 - L;
                m = k;
                g = h;
            }
        }
        k++;
    }

    return f;
}

OP bma(unsigned short s[])
{
  int i, j, k, l, d[K*E] = {0};
  OP lo[K + 1] = {0}, b[K + 1] = {0}, t[K + 1] = {0}, a = {0}, f = {0}, h = {0}, g = {0}, hh = {0};
  vec v = {0}, x = {0}, w = {0};

  //https://www.cayrel.net/?Implementation-of-Goppa-codes
  //unsigned short s[4+1]={0,4,6,3,5};

  //unsigned short s[K+1]={0,13,3,5,4,8,5};
  //unsigned short s[K+1]={0,15,10,8,8,0,12};

  /*
//memset(zz,0,sizeof(zz));
//mkerr(zz,2);
zz[0]=1;
zz[1]=1;
r1=synd(zz);
v=o2v(r1);
memset(s,0,K+2);
s[0]=0;
for(i=0;i<6+1;i++){
s[i+1]=v.x[i];
printf("%d,",s[i]);
}
printf("\n");
*/

  x.x[1] = 1;
  h = v2o(x);
  v.x[0] = 1;
  f = v2o(x);
  lo[0] = v2o(v);
  b[0] = lo[0];

  for (j = 1; j < K + 1; j++)
  {
    v = o2v(lo[j - 1]);
    k = 0;
    printpol(v);
    printf(" ==lo\n");

    l = deg(o2v(lo[j - 1]));
    for (i = 1; i < l + 1; i++)
    {
      k ^= gf[mlt(fg[v.x[i]], fg[s[j - i]])];
      printf("v[%d]=%d\n", i, v.x[i]);
    }
    d[j] = s[j] ^ k;
    printf("d[%d]=%d\n", j, d[j]);
    g = omul(kof(d[j], h), b[j - 1]);
    //  if(2*l>j-1)
    t[j] = oadd(lo[j - 1], g);
    printpol(o2v(t[j]));
    printf("==t[%d]\n", j);
    if (odeg(t[j]) <= odeg(t[j - 1]))
    {
      b[j] = omul(b[j - 1], h);
    }
    else //(deg(o2v(t[j]))>odeg(t[j-1]))
    {
      b[j] = kof(gf[oinv(d[j])], lo[j - 1]);
      l = j - l + 1;
    }
    k = 0;
    if (d[j] > 0)
      lo[j] = t[j];
    if (d[j] == 0)
      lo[j] = lo[j - 1];

    printpol(o2v(b[j]));
    printf(" ==b[%d]\n", j);
  }
  printpol(o2v(lo[j - 1]));
  printf("\n");

  //hh=rev(lo[j-1]);
  //exit(1);
  //r=coeff(r,LT(r).a);
  printpol(o2v(lo[j - 1]));
  printf(" ==coef\n");
  x = chen(lo[j - 1]);
  for (i = 0; i < deg(x) + 1; i++)
  {
    printf("x[%d]=1\n", logx(x.x[i]));

    if (x.x[i] == 0)
      k++;
    if (k > 1)
    {
      printf("baka0\n");
      printvec((x));
      //for (i = 0; i < N; i++)
      //printf("%d,", zz[i]);
      exit(1);
      //return f;
    }
  }

  //return lo[j-1];
}

//言わずもがな
int main(void)
{
  unsigned short z1[N] = {0}; //{1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1};
  OP f = {0}, r = {0}, w = {0};
  vec v = {0};
  unsigned short ss[K] = {0}, zz[N] = {0};

  if (K > N)
  {
    printf("configuration error! K is too big K\n");
    exit(1);
  }

  srand(clock());
  //kabatiansky example
  unsigned short s[K + 1] = {0, 15, 1, 9, 13, 1, 14};
  //Berlekamp-Massey法（UnderConstruction）
  //bms(s);
  //exit(1);

  //公開鍵を生成する(H'=SHP)
 // w = pubkeygen();
    //w=mkpol();
    w=mkg();
  //while(1){

  //エラーベクトルを生成する
  memset(z1, 0, sizeof(z1));
  mkerr(z1, T);
  //exit(1);
  /*
  f=synd(z1);
  v=bms(o2v(f).x);
  chen(v2o(v));
  */
  //exit(1);
  //encryotion test
  //test (w, z1);
  //exit(1);

  //  mkd(w); ??

  memset(zz, 0, sizeof(zz));
  memset(ss, 0, sizeof(ss));

  int j = 0, count = 0;
  //decode開始
  int k = 0;
  while (1)
  {

    memset(z1, 0, sizeof(z1));
    //mkerr(z1, T);
    z1[0]=1;
    z1[1]=2;
    z1[2]=4;
    for (int i = 0; i < N; i++)
    {
      if (z1[i] > 0)
        printf("la=%d %d\n", i, z1[i]);
    }
    //暗号化(v=eH')
    //v = sin2(z1);
    //f=v2o(v);
    //復号(S^-1)
    //f = dec(v.x);
    f=synd(z1);
    r = decode(w, f);
    
    //m=m'P^-1
    count = elo(r);
    for (int i = 0; i < N; i++)
    {
      //検算
      if (z1[i] > 0)
        printf("error position= %d %d\n", i, z1[i]);
    }
    if (count < 0)
    {
      printf("baka-@\n");
      exit(1);
    }
    j++;
    /*
    printf("err=%dっ！！\n", count);
    for (int i = 0; i < N; i++)
      printf("%d,", z1[i]);
    printf("\n");
    */
    exit(1);

    if (j == 10000)
      exit(1);
  }

  return 0;
}

