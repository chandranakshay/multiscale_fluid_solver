#include"perlinNoise.h"

#define B 0x100
#define BM 0xff
#define Number 0x1000
  
static int p[B + B + 2];
static double g3[16][3]={{1.0,1.0,0.0},    {-1.0,1.0,0.0},    {1.0,-1.0,0.0},    {-1.0,-1.0,0.0},
		  { 1.0,0.0,1.0},    {-1.0,0.0,1.0},    {1.0,0.0,-1.0},    {-1.0,0.0,-1.0},
		  { 0.0,1.0,1.0},    {0.0,-1.0,1.0},    {0.0,1.0,-1.0},    {0.0,-1.0,-1.0},
		    {1.0,1.0,0.0},    {0.0,-1.0,1.0},    {-1.0,1.0,0.0},    {0.0,-1.0,-1.0}};

static int start = 1;

static void init(void);

#define s_curve(t) (t * t * t * (t * (t * 6.0 - 15.0) + 10.0))

#define lerp(t, a, b) ( a + t * (b - a) )

#define setup(i,b0,b1,r0,r1)\
	t = vec[i] + Number;\
	b0 = ((int)t) & BM;\
	b1 = (b0+1) & BM;\
	r0 = t - (int)t;\
	r1 = r0 - 1.;

double noise3(double x1, double x2, double x3)
{
	double vec[3];
	vec[0] = x1; vec[1] = x2; vec[2] = x3 ;
	int bx0, bx1, by0, by1, bz0, bz1, b00, b10, b01, b11;
	double rx0, rx1, ry0, ry1, rz0, rz1, *q, sy, sz, a, b, c, d, t, u, v;
	int i, j;

	if (start) {
		start = 0;
		init();
	}

	setup(0, bx0,bx1, rx0,rx1);
	setup(1, by0,by1, ry0,ry1);
	setup(2, bz0,bz1, rz0,rz1);

	i = p[ bx0 ];
	j = p[ bx1 ];

	b00 = p[ i + by0 ];
	b10 = p[ j + by0 ];
	b01 = p[ i + by1 ];
	b11 = p[ j + by1 ];

	t  = s_curve(rx0);
	sy = s_curve(ry0);
	sz = s_curve(rz0);

#define at3(rx,ry,rz) ( rx * q[0] + ry * q[1] + rz * q[2] )

	q = g3[ (b00 + bz0)%16 ] ; u = at3(rx0,ry0,rz0);
	q = g3[( b10 + bz0)%16 ] ; v = at3(rx1,ry0,rz0);
	a = lerp(t, u, v);

	q = g3[ (b01 + bz0)%16 ] ; u = at3(rx0,ry1,rz0);
	q = g3[ (b11 + bz0)%16 ] ; v = at3(rx1,ry1,rz0);
	b = lerp(t, u, v);

	c = lerp(sy, a, b);

	q = g3[ (b00 + bz1)%16 ] ; u = at3(rx0,ry0,rz1);
	q = g3[ (b10 + bz1)%16 ] ; v = at3(rx1,ry0,rz1);
	a = lerp(t, u, v);

	q = g3[ (b01 + bz1)%16 ] ; u = at3(rx0,ry1,rz1);
	q = g3[ (b11 + bz1)%16 ] ; v = at3(rx1,ry1,rz1);
	b = lerp(t, u, v);

	d = lerp(sy, a, b);

	return lerp(sz, c, d);
}

static void normalize3(double v[3])
{
	double s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	v[0] = v[0] / s;
	v[1] = v[1] / s;
	v[2] = v[2] / s;
}

static void init(void)
{
	int i, j, k;
	srand((unsigned)time(0));
	for (i = 0 ; i < B ; i++) p[i] = i;
		    
	while (--i) {
		k = p[i];
		p[i] = p[j = random() % B];
		p[j] = k;
	}
	
	for (i = 0 ; i < 16 ; i++)  normalize3(g3[i]);
	for (i = 0 ; i < B + 2 ; i++) 	p[B + i] = p[i];
}

