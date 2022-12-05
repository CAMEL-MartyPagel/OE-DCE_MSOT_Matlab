__constant sampler_t iSamplerNCL = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP | CLK_FILTER_LINEAR; 
__constant sampler_t iSamplerNEL = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR; 
__constant sampler_t iSamplerNCN = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST; 
__constant sampler_t iSamplerUCN = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP | CLK_FILTER_NEAREST; 
__constant sampler_t iSamplerUEN = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST; 

#define ATTENUATION_VALUE 0.0f
#define FSAMPLE 4.0e7f
#define IFSAMPLE 2.5e-8f
#define __MAX_COMPONENTS_COUNT__ 10
#define MARGIN 1e-5f


#ifdef __SPECIAL_CLAMP__

#define VALUE_CLAMP_LO 120.0f
#define VALUE_CLAMP_HI 2029.0f

__attribute__((always_inline))
float4 CLAMP_TOF_INDEX4(float4 x)
{
	x = select(x, (float4)VALUE_CLAMP_LO, as_int4(x - (float4)120.0f));
	return select(x, (float4)VALUE_CLAMP_HI, as_int4((float4)2029.0f - x));
}

__attribute__((always_inline))
float CLAMP_TOF_INDEX(float x)
{
	x = select(x, VALUE_CLAMP_LO, as_int(x - 120.0f));
	return select(x, VALUE_CLAMP_HI, as_int(2029.0f - x));
}
#else
__attribute__((always_inline))
float4 CLAMP_TOF_INDEX4(float4 x)
{
	return clamp(x, 120.0f,2029.0f);
}

__attribute__((always_inline))
float CLAMP_TOF_INDEX(float x)
{
	return clamp(x, 120.0f,2029.0f);
}

#endif

#ifdef FP_FAST_FMAF
#define FMA(a,b,c) fma(a,b,c)
#else
#define FMA(a,b,c) ((a)*(b)+(c))
#endif


#define ALPHA 0.5f
__attribute__((always_inline))
float WCUBIC01(float x) 
{
	return ((2.0f-ALPHA)*x+(ALPHA-3.0f)) * x * x + 1.0f;
}

__attribute__((always_inline))
float WCUBIC12(float x) 
{
	return (((5.0f-x)*x - 8.0f)*x + 4.0f)*ALPHA;
}

__attribute__((always_inline))
float BICUBIC(__read_only image2d_t input, float2 coordin)
{
	float x, y, x_floor, p, q, u, v;
	int2 icoord;
	x = coordin.x;
	y = coordin.y;
	coordin = floor(coordin)-1.0f;
	u = x_floor = coordin.x;
	v = coordin.y;
	icoord = convert_int2(coordin);
	p = read_imagef(input, iSamplerUEN, icoord).x * WCUBIC12(x-u);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC01(x-u);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC01(u-x);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC12(u-x);
	q = p*WCUBIC12(y-v);
	v = v + 1;	u = x_floor;
	icoord = convert_int2((float2)(u,v));
	p = read_imagef(input, iSamplerUEN, icoord).x * WCUBIC12(x-u);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC01(x-u);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC01(u-x);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC12(u-x);
	q += p*WCUBIC01(y-v);
	v = v + 1;	u = x_floor;
	icoord = convert_int2((float2)(u,v));
	p = read_imagef(input, iSamplerUEN, icoord).x * WCUBIC12(x-u);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC01(x-u);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC01(u-x);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC12(u-x);
	q += p*WCUBIC01(v-y);
	v = v + 1;	u = x_floor;
	icoord = convert_int2((float2)(u,v));
	p = read_imagef(input, iSamplerUEN, icoord).x * WCUBIC12(x-u);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC01(x-u);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC01(u-x);
	u = u + 1;	icoord.x++;
	p += read_imagef(input, iSamplerUEN, icoord).x * WCUBIC12(u-x);
	q += p*WCUBIC12(v-y);
	return q;
}

__kernel
void interp2d(__write_only image2d_t output, __read_only image2d_t input, float2 scale, int4 mn)
{
	float4 q;
	int2 coordout = (int2)(get_global_id(0), get_global_id(1));
	if(all(coordout - mn.xy)==0)
	{
		return;
	}
	float2 coordin = convert_float2(coordout)*scale;
	if(mn.w == 1)
	{
		int2 icoord = convert_int2(coordin+0.5f);
		q = read_imagef(input, iSamplerUEN, icoord);
	}
	else if(mn.w == 2)
	{
		float x, y, p;

		x = coordin.x;
		y = coordin.y;
		coordin = floor(coordin);
		x -= coordin.x; 
		y -= coordin.y;
		int2 icoord = convert_int2(coordin);

		p = read_imagef(input, iSamplerUEN, icoord).x * (1.0f-x);
		icoord.x++;
		p += read_imagef(input, iSamplerUEN, icoord).x * x;
		q = p*(1.0f-y);
		icoord.x--;
		icoord.y++;
		p = read_imagef(input, iSamplerUEN, icoord).x * (1.0f-x);
		icoord.x++;
		p += read_imagef(input, iSamplerUEN, icoord).x * x;
		q += p*y;
	}
	else
	{
		q = BICUBIC(input, coordin);
	}
	q.w =1.0f;
	write_imagef(output, coordout, q);
}


__kernel
void LinearInterp2d(__write_only image2d_t output, __read_only image2d_t input, float2 scale)
{
	int2 coordout = (int2)(get_global_id(0), get_global_id(1));
	if(all(coordout - get_image_dim(output))==0)
	{
		return;
	}
	float2 coordin = convert_float2(coordout)*scale;
	float4 p = read_imagef(input, iSamplerNEL, coordin);
	write_imagef(output, coordout, p);
}

__kernel
void BiCubicInterp2d(__write_only image2d_t output, __read_only image2d_t input, float2 scale)
{
	float4 q;
	int2 coordout = (int2)(get_global_id(0), get_global_id(1));
	if(all(coordout - get_image_dim(output)) == 0)
	{
		return;
	}
	float2 coordin = convert_float2(coordout)*scale;

	q = BICUBIC(input, coordin);
	q.w = 1.0f;
	write_imagef(output, coordout, q);
}

__kernel
void BiCubicOffsetInterp2d(__write_only image2d_t output, __read_only image2d_t input, float2 scale, float2 offset)
{
	float4 q;
	int2 coordout = (int2)(get_global_id(0), get_global_id(1));
	if(all(coordout - get_image_dim(output)) == 0)
	{
		return;
	}
	float2 coordin = (convert_float2(coordout))*scale+offset;

	q = BICUBIC(input, coordin);
	q.w = 1.0f;
	write_imagef(output, coordout, q);
}

__kernel 
void copyBufferRectToImage(__global float *buffer, __write_only image2d_t image, int2 src_origin, int src_row_pitch)
{
	int2 id = (int2)(get_global_id(0), get_global_id(1));
	if(all(id - get_image_dim(image)) == 0) return;
	int2 src_coordinates = src_origin + id;
	int offset = src_coordinates.x + src_coordinates.y * src_row_pitch;
	float pixel = buffer[offset];
	write_imagef(image, id, (float4)pixel);
}

__attribute__((always_inline))
float2 conj(float2 a)
{
	return (float2)(a.x, -a.y);
}

__attribute__((always_inline))
float2 conjTransp(float2 a)
{
	return (float2)(-a.y, a.x);
}

__attribute__((always_inline))
float2 complexMul(float2 a, float2 b)
{
//    return (float2)(a.x * b.x - a.y * b.y, a.y * b.x + a.x * b.y);
//    return (float2)(mad(-a.y, b.y, a.x * b.x), mad(a.y, b.x, a.x * b.y));
#ifdef FP_FAST_FMAF
	return (float2)(fma(-a.y, b.y, a.x * b.x), fma(a.y, b.x, a.x * b.y));
#else
//    return (float2)(a.x * b.x - a.y * b.y, a.y * b.x + a.x * b.y);
	return (float2)(mad(-a.y, b.y, a.x * b.x), mad(a.y, b.x, a.x * b.y));
#endif
}

__attribute__((always_inline))
void fftKernel2(float2 *a, int i0, int i1)
{
	float2 c = a[i0];
	a[i0] = c + a[i1];
	a[i1] = c - a[i1];
}

__attribute__((always_inline))
void fftKernel4(float2 *a, int i0, int i1, int i2, int i3, float dir)
{
	fftKernel2(a, i0, i2);
	fftKernel2(a, i1, i3);
	fftKernel2(a, i0, i1);
	a[i3] = (float2)(dir)*conjTransp(a[i3]);
	fftKernel2(a, i2, i3);
	float2 c = a[i1];
	a[i1] = a[i2];
	a[i2] = c; 
}

__attribute__((always_inline))
void bitreverse8(float2 *a)
{
	float2 c;
	c = a[1];    a[1] = a[4];    a[4] = c;
	c = a[3];    a[3] = a[6];    a[6] = c;
}

__attribute__((always_inline))
void fftKernel8(float2 *a, float dir)
{
	const float2 w1 = (float2)(0x1.6a09e6p-1f,  dir*0x1.6a09e6p-1f);
	const float2 w3 = (float2)(-0x1.6a09e6p-1f, dir*0x1.6a09e6p-1f);
	fftKernel2(a, 0, 4);
	fftKernel2(a, 1, 5);
	fftKernel2(a, 2, 6);
	fftKernel2(a, 3, 7);
	a[5] = complexMul(w1, a[5]);
	a[6] = (float2)(dir)*conjTransp(a[6]);
	a[7] = complexMul(w3, a[7]);
	fftKernel2(a, 0, 2);
	fftKernel2(a, 1, 3);
	fftKernel2(a, 4, 6);
	fftKernel2(a, 5, 7);
	a[3] = (float2)(dir)*conjTransp(a[3]);
	a[7] = (float2)(dir)*conjTransp(a[7]);
	fftKernel2(a, 0, 1);
	fftKernel2(a, 2, 3);
	fftKernel2(a, 4, 5);
	fftKernel2(a, 6, 7);
	bitreverse8(a);
}

__attribute__((always_inline))
void bitreverse4x4(float2 *a)
{
	float2 c;
	c = a[1];  a[1]  = a[4];  a[4]  = c;
	c = a[2];  a[2]  = a[8];  a[8]  = c;
	c = a[3];  a[3]  = a[12]; a[12] = c;
	c = a[6];  a[6]  = a[9];  a[9]  = c;
	c = a[7];  a[7]  = a[13]; a[13] = c;
	c = a[11]; a[11] = a[14]; a[14] = c;
}

__attribute__((always_inline))
void fftKernel16(float2 *a, float dir)
{
	const float w0 = 0x1.d906bcp-1f;
	const float w1 = 0x1.87de2ap-2f;
	const float w2 = 0x1.6a09e6p-1f;
	fftKernel4(a, 0, 4, 8, 12, dir);
	fftKernel4(a, 1, 5, 9, 13, dir);
	fftKernel4(a, 2, 6, 10, 14, dir);
	fftKernel4(a, 3, 7, 11, 15, dir);
	a[5]  = complexMul(a[5], (float2)(w0, dir*w1));
	a[6]  = complexMul(a[6], (float2)(w2, dir*w2));
	a[7]  = complexMul(a[7], (float2)(w1, dir*w0));
	a[9]  = complexMul(a[9], (float2)(w2, dir*w2));
	a[10] = (float2)(dir)*conjTransp(a[10]);
	a[11] = complexMul(a[11], (float2)(-w2, dir*w2));
	a[13] = complexMul(a[13], (float2)(w1, dir*w0));
	a[14] = complexMul(a[14], (float2)(-w2, dir*w2));
	a[15] = complexMul(a[15], (float2)(-w0, dir*-w1));
	fftKernel4(a, 0, 1, 2, 3, dir);
	fftKernel4(a, 4, 5, 6, 7, dir);
	fftKernel4(a, 8, 9, 10, 11, dir);
	fftKernel4(a, 12, 13, 14, 15, dir);
	bitreverse4x4(a);
}

__attribute__((always_inline))
void outUpdate(__global float *out_real, __global float *out_imag, int offset, float2 v)
{
	out_real[offset] = v.x;
	out_imag[offset] = v.y;
}


/* Run kernel fft_4096_0 with global dim = {256*BatchSize}, local dim={128} */

#define SCALE_IFFT4096 (1.0f/4096.0f)



////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/* Run kernel fft_4096_0 with global dim = {256*BatchSize}, local dim={128} */
__kernel __attribute__((reqd_work_group_size (128,1,1)))
void filter_fft_4096_0(__global float *sigMat, __global float2 *out, float pwrCorr)
{
    __local float2 sMem[2064];
    int i, j, indexIn, indexOut, tid, bNum, xNum, k, l;
    float2 w, w0;
    float ang, ang1;
    __local float2 *lMemStore, *lMemLoad;
    __global float *in_real, *in_imag;
    float2 a[16];
    int groupId = get_group_id( 0 );
    bNum = groupId & 1;
    xNum = (groupId >> 1) << 12;
    indexIn = mul24(bNum, 16);
    indexOut = mad24(indexIn, 128, xNum);
    tid = get_local_id( 0 );
    i = tid & 15;
    j = tid >> 4;
    indexIn += mad24(j, 32, i+xNum);
    in_real = sigMat + indexIn;
    in_imag = in_real + 2048;

    a[0] = (float2)(in_real[0], in_imag[0]) * pwrCorr;
    a[1] = (float2)(in_real[256], in_imag[256]) * pwrCorr;
    a[2] = (float2)(in_real[512], in_imag[512]) * pwrCorr;
    a[3] = (float2)(in_real[768], in_imag[768]) * pwrCorr;
    a[4] = (float2)(in_real[1024], in_imag[1024]) * pwrCorr;
    a[5] = (float2)(in_real[1280], in_imag[1280]) * pwrCorr;
    a[6] = (float2)(in_real[1536], in_imag[1536]) * pwrCorr;
    a[7] = (float2)(in_real[1792], in_imag[1792]) * pwrCorr;
    a[8] = 0;
    a[9] = 0;
    a[10] = 0;
    a[11] = 0;
    a[12] = 0;
    a[13] = 0;
    a[14] = 0;
    a[15] = 0;

    fftKernel16(a, -1.0f);
    ang = -(M_PI_F/64.0f)*j;
    w0 = (float2)(native_cos(ang), native_sin(ang));
    w = w0;
    lMemStore = sMem + tid;
    lMemStore[0] = a[0];
    lMemStore[128] = complexMul(a[1], w);    w = complexMul(w0, w);
    lMemStore[256] = complexMul(a[2], w);    w = complexMul(w0, w);
    lMemStore[384] = complexMul(a[3], w);    w = complexMul(w0, w);
    lMemStore[512] = complexMul(a[4], w);    w = complexMul(w0, w);
    lMemStore[640] = complexMul(a[5], w);    w = complexMul(w0, w);
    lMemStore[768] = complexMul(a[6], w);    w = complexMul(w0, w);
    lMemStore[896] = complexMul(a[7], w);    w = complexMul(w0, w);
    lMemStore[1024] = complexMul(a[8], w);    w = complexMul(w0, w);
    lMemStore[1152] = complexMul(a[9], w);    w = complexMul(w0, w);
    lMemStore[1280] = complexMul(a[10], w);  w = complexMul(w0, w);
    lMemStore[1408] = complexMul(a[11], w);  w = complexMul(w0, w);
    lMemStore[1536] = complexMul(a[12], w);  w = complexMul(w0, w);
    lMemStore[1664] = complexMul(a[13], w);  w = complexMul(w0, w);
    lMemStore[1792] = complexMul(a[14], w);  w = complexMul(w0, w);
    lMemStore[1920] = complexMul(a[15], w);
    barrier(CLK_LOCAL_MEM_FENCE);

    indexIn = mad24(j, 256, i);
    lMemLoad = sMem + indexIn;
    a[0] = lMemLoad[0];
    a[1] = lMemLoad[16];
    a[2] = lMemLoad[32];
    a[3] = lMemLoad[48];
    a[4] = lMemLoad[64];
    a[5] = lMemLoad[80];
    a[6] = lMemLoad[96];
    a[7] = lMemLoad[112];
    a[8] = lMemLoad[128];
    a[9] = lMemLoad[144];
    a[10] = lMemLoad[160];
    a[11] = lMemLoad[176];
    a[12] = lMemLoad[192];
    a[13] = lMemLoad[208];
    a[14] = lMemLoad[224];
    a[15] = lMemLoad[240];
    barrier(CLK_LOCAL_MEM_FENCE);

    fftKernel8(a + 0, -1.0f);
    fftKernel8(a + 8, -1.0f);
    l = (bNum << 4) + i;
    k = j << 1;
    ang1 = -(M_PI_F/2048.0f)*l;
    w0 = (float2)(native_cos(ang1 * 16.0f), native_sin(ang1 * 16.0f));
    ang = ang1*k;
    w = (float2)(native_cos(ang), native_sin(ang));
    lMemStore = sMem + mad24(i, 129, j << 1);
    lMemStore[0] = complexMul(a[0], w);    w = complexMul(w0, w);
    lMemStore[16] = complexMul(a[1], w);    w = complexMul(w0, w);
    lMemStore[32] = complexMul(a[2], w);    w = complexMul(w0, w);
    lMemStore[48] = complexMul(a[3], w);    w = complexMul(w0, w);
    lMemStore[64] = complexMul(a[4], w);    w = complexMul(w0, w);
    lMemStore[80] = complexMul(a[5], w);    w = complexMul(w0, w);
    lMemStore[96] = complexMul(a[6], w);    w = complexMul(w0, w);
    lMemStore[112] = complexMul(a[7], w);
    ang += ang1;
    w = (float2)(native_cos(ang), native_sin(ang));
    lMemStore[1] = complexMul(a[8], w);    w = complexMul(w0, w);
    lMemStore[17] = complexMul(a[9], w);    w = complexMul(w0, w);
    lMemStore[33] = complexMul(a[10], w);    w = complexMul(w0, w);
    lMemStore[49] = complexMul(a[11], w);    w = complexMul(w0, w);
    lMemStore[65] = complexMul(a[12], w);    w = complexMul(w0, w);
    lMemStore[81] = complexMul(a[13], w);    w = complexMul(w0, w);
    lMemStore[97] = complexMul(a[14], w);    w = complexMul(w0, w);
    lMemStore[113] = complexMul(a[15], w);
    barrier(CLK_LOCAL_MEM_FENCE);
    lMemLoad = sMem + mad24(tid >> 7, 129, tid & 127);
    indexOut += tid;
    out += indexOut;
    out[   0] = lMemLoad[0];
    out[ 128] = lMemLoad[129];
    out[ 256] = lMemLoad[258];
    out[ 384] = lMemLoad[387];
    out[ 512] = lMemLoad[516];
    out[ 640] = lMemLoad[645];
    out[ 768] = lMemLoad[774];
    out[ 896] = lMemLoad[903];
    out[1024] = lMemLoad[1032];
    out[1152] = lMemLoad[1161];
    out[1280] = lMemLoad[1290];
    out[1408] = lMemLoad[1419];
    out[1536] = lMemLoad[1548];
    out[1664] = lMemLoad[1677];
    out[1792] = lMemLoad[1806];
    out[1920] = lMemLoad[1935];
}

/* Run kernel fft_32_1 with global dim = {512*BatchSize}, local dim={64} */
__kernel __attribute__((reqd_work_group_size (64,1,1)))
void filter_fft_32_1(__global float2 *in, __global float2 *out)
{
    __local float2 sMem[512];
    int i, j, indexIn, indexOut, tid, bNum, xNum;
    float2 w, w0;
    float ang;
    __local float2 *lMemStore, *lMemLoad;

    float2 a[8];
    int lId = get_local_id( 0 );
    int groupId = get_group_id( 0 );
    bNum = groupId & 7;
    xNum = (groupId >> 3) << 12;
    indexIn = mul24(bNum, 16);
    tid = indexIn;
    i = tid >> 7;
    j = tid & 127;
    indexOut = mad24(i, 4096, j+xNum);
    tid = lId;
    i = tid & 15;
    j = tid >> 4;
    indexIn += mad24(j, 128, i+xNum);
    in += indexIn;
    a[0] = in[0];
    a[1] = in[512];
    a[2] = in[1024];
    a[3] = in[1536];
    a[4] = in[2048];
    a[5] = in[2560];
    a[6] = in[3072];
    a[7] = in[3584];
    fftKernel8(a, -1.0f);
    ang = -(M_PI_F/16.0f)*j;
    w0 = (float2)(native_cos(ang), native_sin(ang));
    w = w0;
    lMemStore = sMem + tid;
    lMemStore[0] = a[0];
    lMemStore[64] = complexMul(a[1], w);    w = complexMul(w0, w);
    lMemStore[128] = complexMul(a[2], w);    w = complexMul(w0, w);
    lMemStore[192] = complexMul(a[3], w);    w = complexMul(w0, w);
    lMemStore[256] = complexMul(a[4], w);    w = complexMul(w0, w);
    lMemStore[320] = complexMul(a[5], w);    w = complexMul(w0, w);
    lMemStore[384] = complexMul(a[6], w);    w = complexMul(w0, w);
    lMemStore[448] = complexMul(a[7], w);
    barrier(CLK_LOCAL_MEM_FENCE);
    indexIn = mad24(j, 128, i);
    lMemLoad = sMem + indexIn;
    a[0] = lMemLoad[0];
    a[1] = lMemLoad[16];
    a[2] = lMemLoad[32];
    a[3] = lMemLoad[48];
    a[4] = lMemLoad[64];
    a[5] = lMemLoad[80];
    a[6] = lMemLoad[96];
    a[7] = lMemLoad[112];
    fftKernel4(a, 0, 1, 2, 3, -1.0f);
    fftKernel4(a, 4, 5, 6, 7, -1.0f);
    indexOut += mad24(j, 256, i);
    out += indexOut;

	out[0] 		= a[0];
	out[1024]	= a[1];
	out[2048]	= a[2];
	out[3072]	= a[3];
	out[128]	= a[4];
	out[1152]	= a[5];
	out[2176]	= a[6];
	out[3200]	= a[7];
}

/* Run kernel fft_4096_0 with global dim = {256*BatchSize}, local dim={128} */
__kernel __attribute__((reqd_work_group_size (128,1,1)))
void filter_ifft_4096_0(__global float2 *in, __global float2 *spectra, __global float2 *out)
{
    __local float2 sMem[2064];
    int i, j, indexIn, indexOut, tid, bNum, xNum, k, l;
    float2 w, w0;
    float ang, ang1;
    __local float2 *lMemStore, *lMemLoad;
    float2 a[16];
    int groupId = get_group_id( 0 );
    bNum = groupId & 1;
    xNum = (groupId >> 1) << 12;
    indexIn = mul24(bNum, 16);
    indexOut = mad24(indexIn, 128, xNum);
    tid = get_local_id( 0 );
    i = tid & 15;
    j = tid >> 4;
    indexIn += mad24(j, 32, i);
    in += indexIn+xNum;
	spectra += indexIn;

	w = spectra[   0];    w0 = in[   0];	a[0]  = complexMul(w,w0);
	w = spectra[ 256];    w0 = in[ 256];	a[1]  = complexMul(w,w0);
	w = spectra[ 512];    w0 = in[ 512];	a[2]  = complexMul(w,w0);
	w = spectra[ 768];    w0 = in[ 768];	a[3]  = complexMul(w,w0);
	w = spectra[1024];    w0 = in[1024];	a[4]  = complexMul(w,w0);
	w = spectra[1280];    w0 = in[1280];	a[5]  = complexMul(w,w0);
	w = spectra[1536];    w0 = in[1536];	a[6]  = complexMul(w,w0);
	w = spectra[1792];    w0 = in[1792];	a[7]  = complexMul(w,w0);
	w = spectra[2048];    w0 = in[2048];	a[8]  = complexMul(w,w0);
	w = spectra[2304];    w0 = in[2304];	a[9]  = complexMul(w,w0);
    w = spectra[2560];    w0 = in[2560];	a[10] = complexMul(w,w0);
    w = spectra[2816];    w0 = in[2816];	a[11] = complexMul(w,w0);
    w = spectra[3072];    w0 = in[3072];	a[12] = complexMul(w,w0);
    w = spectra[3328];    w0 = in[3328];	a[13] = complexMul(w,w0);
    w = spectra[3584];    w0 = in[3584];	a[14] = complexMul(w,w0);
    w = spectra[3840];    w0 = in[3840];	a[15] = complexMul(w,w0);
		
    fftKernel16(a, 1.0f);
    ang = (M_PI_F/64.0f)*j;
    w0 = (float2)(native_cos(ang), native_sin(ang));
    w = w0;
    lMemStore = sMem + tid;
    lMemStore[0] = a[0];
    lMemStore[128] = complexMul(a[1], w);    w = complexMul(w0, w);
    lMemStore[256] = complexMul(a[2], w);    w = complexMul(w0, w);
    lMemStore[384] = complexMul(a[3], w);    w = complexMul(w0, w);
    lMemStore[512] = complexMul(a[4], w);    w = complexMul(w0, w);
    lMemStore[640] = complexMul(a[5], w);    w = complexMul(w0, w);
    lMemStore[768] = complexMul(a[6], w);    w = complexMul(w0, w);
    lMemStore[896] = complexMul(a[7], w);    w = complexMul(w0, w);
    lMemStore[1024] = complexMul(a[8], w);    w = complexMul(w0, w);
    lMemStore[1152] = complexMul(a[9], w);    w = complexMul(w0, w);
    lMemStore[1280] = complexMul(a[10], w);  w = complexMul(w0, w);
    lMemStore[1408] = complexMul(a[11], w);  w = complexMul(w0, w);
    lMemStore[1536] = complexMul(a[12], w);  w = complexMul(w0, w);
    lMemStore[1664] = complexMul(a[13], w);  w = complexMul(w0, w);
    lMemStore[1792] = complexMul(a[14], w);  w = complexMul(w0, w);
    lMemStore[1920] = complexMul(a[15], w);
    barrier(CLK_LOCAL_MEM_FENCE);

    indexIn = mad24(j, 256, i);
    lMemLoad = sMem + indexIn;
    a[0] = lMemLoad[0];
    a[1] = lMemLoad[16];
    a[2] = lMemLoad[32];
    a[3] = lMemLoad[48];
    a[4] = lMemLoad[64];
    a[5] = lMemLoad[80];
    a[6] = lMemLoad[96];
    a[7] = lMemLoad[112];
    a[8] = lMemLoad[128];
    a[9] = lMemLoad[144];
    a[10] = lMemLoad[160];
    a[11] = lMemLoad[176];
    a[12] = lMemLoad[192];
    a[13] = lMemLoad[208];
    a[14] = lMemLoad[224];
    a[15] = lMemLoad[240];
    barrier(CLK_LOCAL_MEM_FENCE);

    fftKernel8(a + 0, 1.0f);
    fftKernel8(a + 8, 1.0f);
    l = (bNum << 4) + i;
    k = j << 1;
    ang1 = (M_PI_F/2048.0f)*l;
    w0 = (float2)(native_cos(ang1 * 16.0f), native_sin(ang1 * 16.0f));
    ang = ang1*k;
    w = (float2)(native_cos(ang), native_sin(ang));
    lMemStore = sMem + mad24(i, 129, j << 1);
    lMemStore[0] = complexMul(a[0], w);    w = complexMul(w0, w);
    lMemStore[16] = complexMul(a[1], w);    w = complexMul(w0, w);
    lMemStore[32] = complexMul(a[2], w);    w = complexMul(w0, w);
    lMemStore[48] = complexMul(a[3], w);    w = complexMul(w0, w);
    lMemStore[64] = complexMul(a[4], w);    w = complexMul(w0, w);
    lMemStore[80] = complexMul(a[5], w);    w = complexMul(w0, w);
    lMemStore[96] = complexMul(a[6], w);    w = complexMul(w0, w);
    lMemStore[112] = complexMul(a[7], w);
    ang += ang1;
    w = (float2)(native_cos(ang), native_sin(ang));
    lMemStore[1] = complexMul(a[8], w);    w = complexMul(w0, w);
    lMemStore[17] = complexMul(a[9], w);    w = complexMul(w0, w);
    lMemStore[33] = complexMul(a[10], w);    w = complexMul(w0, w);
    lMemStore[49] = complexMul(a[11], w);    w = complexMul(w0, w);
    lMemStore[65] = complexMul(a[12], w);    w = complexMul(w0, w);
    lMemStore[81] = complexMul(a[13], w);    w = complexMul(w0, w);
    lMemStore[97] = complexMul(a[14], w);    w = complexMul(w0, w);
    lMemStore[113] = complexMul(a[15], w);
    barrier(CLK_LOCAL_MEM_FENCE);
	
    lMemLoad = sMem + mad24(tid >> 7, 129, tid & 127);
    indexOut += tid;
    out += indexOut;
    out[   0] = lMemLoad[0];
    out[ 128] = lMemLoad[129];
    out[ 256] = lMemLoad[258];
    out[ 384] = lMemLoad[387];
    out[ 512] = lMemLoad[516];
    out[ 640] = lMemLoad[645];
    out[ 768] = lMemLoad[774];
    out[ 896] = lMemLoad[903];
    out[1024] = lMemLoad[1032];
    out[1152] = lMemLoad[1161];
    out[1280] = lMemLoad[1290];
    out[1408] = lMemLoad[1419];
    out[1536] = lMemLoad[1548];
    out[1664] = lMemLoad[1677];
    out[1792] = lMemLoad[1806];
    out[1920] = lMemLoad[1935];
}

/* Run kernel fft_32_1 with global dim = {512*BatchSize}, local dim={64} */
__kernel __attribute__((reqd_work_group_size (64,1,1)))
void filter_ifft_32_1(__global float2 *in, __global float *sigMat)
{
    __local float2 sMem[512];
    int i, j, indexIn, indexOut, tid, bNum, xNum;
    float2 w, w0;
    float ang;
    __local float2 *lMemStore, *lMemLoad;

    float2 a[8];
    int lId = get_local_id( 0 );
    int groupId = get_group_id( 0 );
    bNum = groupId & 7;
    xNum = (groupId >> 3) << 12;
    indexIn = mul24(bNum, 16);
    tid = indexIn;
    i = tid >> 7;
    j = tid & 127;
    indexOut = mad24(i, 4096, j+xNum);
    tid = lId;
    i = tid & 15;
    j = tid >> 4;
    indexIn += mad24(j, 128, i+xNum);
    in += indexIn;
    a[0] = in[0];
    a[1] = in[512];
    a[2] = in[1024];
    a[3] = in[1536];
    a[4] = in[2048];
    a[5] = in[2560];
    a[6] = in[3072];
    a[7] = in[3584];
    fftKernel8(a, 1.0f);
    ang = (M_PI_F/16.0f)*j;
    w0 = (float2)(native_cos(ang), native_sin(ang));
    w = w0;
    lMemStore = sMem + tid;
    lMemStore[0] = a[0];
    lMemStore[64] = complexMul(a[1], w);    w = complexMul(w0, w);
    lMemStore[128] = complexMul(a[2], w);    w = complexMul(w0, w);
    lMemStore[192] = complexMul(a[3], w);    w = complexMul(w0, w);
    lMemStore[256] = complexMul(a[4], w);    w = complexMul(w0, w);
    lMemStore[320] = complexMul(a[5], w);    w = complexMul(w0, w);
    lMemStore[384] = complexMul(a[6], w);    w = complexMul(w0, w);
    lMemStore[448] = complexMul(a[7], w);
    barrier(CLK_LOCAL_MEM_FENCE);
    indexIn = mad24(j, 128, i);
    lMemLoad = sMem + indexIn;
    a[0] = lMemLoad[0];
    a[1] = lMemLoad[16];
    a[2] = lMemLoad[32];
    a[3] = lMemLoad[48];
    a[4] = lMemLoad[64];
    a[5] = lMemLoad[80];
    a[6] = lMemLoad[96];
    a[7] = lMemLoad[112];
    fftKernel4(a, 0, 1, 2, 3, 1.0f);
    fftKernel4(a, 4, 5, 6, 7, 1.0f);
    indexOut += mad24(j, 256, i);
	__global float* out_real = sigMat + indexOut;
	__global float* out_imag = out_real + 2048;

	outUpdate(out_real, out_imag, 0, a[0] * SCALE_IFFT4096);
	outUpdate(out_real, out_imag, 128, a[4] * SCALE_IFFT4096);
	outUpdate(out_real, out_imag, 1024, a[1] * SCALE_IFFT4096);
	outUpdate(out_real, out_imag, 1152, a[5] * SCALE_IFFT4096);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

__kernel __attribute__((reqd_work_group_size (128,1,1)))
void ri_fft_4096_0(__global float *in_real, __global float2 *out)
{
    __local float2 sMem[2064];
    int i, j, indexIn, indexOut, tid, bNum, xNum, k, l;
    float2 w, w0;
    float ang, ang1;
    __local float2 *lMemStore, *lMemLoad;
    float2 a[16];
    int groupId = get_group_id( 0 );
    bNum = groupId & 1;
    xNum = (groupId >> 1) << 12;
    indexIn = mul24(bNum, 16);
	indexOut = mad24(indexIn, 128, xNum);
	tid = get_local_id(0);
	i = tid & 15;
	j = tid >> 4;
	indexIn += mad24(j, 32, i+xNum);
	in_real += indexIn;
	
    a[0] = (float2)(in_real[0],0);
    a[1] = (float2)(in_real[256],0);
    a[2] = (float2)(in_real[512],0);
    a[3] = (float2)(in_real[768],0);
    a[4] = (float2)(in_real[1024],0);
    a[5] = (float2)(in_real[1280],0);
    a[6] = (float2)(in_real[1536],0);
    a[7] = (float2)(in_real[1792],0);
	a[8] = (float2)(in_real[2048], 0);
	a[9] = (float2)(in_real[2304], 0);
	a[10] = (float2)(in_real[2560], 0);
	a[11] = (float2)(in_real[2816], 0);
	a[12] = (float2)(in_real[3072], 0);
	a[13] = (float2)(in_real[3328], 0);
	a[14] = (float2)(in_real[3584], 0);
	a[15] = (float2)(in_real[3840], 0);

    fftKernel16(a, -1.0f);
    ang = -(M_PI_F/64.0f)*j;
    w0 = (float2)(native_cos(ang), native_sin(ang));
    w = w0;
	lMemStore = sMem + tid;
	lMemStore[0] = a[0];
	lMemStore[ 128] = complexMul(a[ 1], w); w = complexMul(w0, w);
	lMemStore[ 256] = complexMul(a[ 2], w); w = complexMul(w0, w);
	lMemStore[ 384] = complexMul(a[ 3], w); w = complexMul(w0, w);
	lMemStore[ 512] = complexMul(a[ 4], w); w = complexMul(w0, w);
	lMemStore[ 640] = complexMul(a[ 5], w); w = complexMul(w0, w);
	lMemStore[ 768] = complexMul(a[ 6], w); w = complexMul(w0, w);
	lMemStore[ 896] = complexMul(a[ 7], w); w = complexMul(w0, w);
	lMemStore[1024] = complexMul(a[ 8], w);	w = complexMul(w0, w);
	lMemStore[1152] = complexMul(a[ 9], w);	w = complexMul(w0, w);
	lMemStore[1280] = complexMul(a[10], w); w = complexMul(w0, w);
	lMemStore[1408] = complexMul(a[11], w); w = complexMul(w0, w);
	lMemStore[1536] = complexMul(a[12], w); w = complexMul(w0, w);
	lMemStore[1664] = complexMul(a[13], w); w = complexMul(w0, w);
	lMemStore[1792] = complexMul(a[14], w); w = complexMul(w0, w);
	lMemStore[1920] = complexMul(a[15], w);
	barrier(CLK_LOCAL_MEM_FENCE);
	
	indexIn = mad24(j, 256, i);
	lMemLoad = sMem + indexIn;
	a[0] = lMemLoad[0];
	a[1] = lMemLoad[16];
	a[2] = lMemLoad[32];
	a[3] = lMemLoad[48];
	a[4] = lMemLoad[64];
	a[5] = lMemLoad[80];
	a[6] = lMemLoad[96];
	a[7] = lMemLoad[112];
	a[8] = lMemLoad[128];
	a[9] = lMemLoad[144];
	a[10] = lMemLoad[160];
	a[11] = lMemLoad[176];
	a[12] = lMemLoad[192];
	a[13] = lMemLoad[208];
	a[14] = lMemLoad[224];
	a[15] = lMemLoad[240];
	barrier(CLK_LOCAL_MEM_FENCE);
	
    fftKernel8(a + 0, -1.0f);
    fftKernel8(a + 8, -1.0f);
    l = (bNum << 4) + i;
    k = j << 1;
    ang1 = -(M_PI_F/2048.0f)*l;
    w0 = (float2)(native_cos(ang1 * 16.0f), native_sin(ang1 * 16.0f));
    ang = ang1*k;
    w = (float2)(native_cos(ang), native_sin(ang));
	lMemStore = sMem + mad24(i, 129, j << 1);
	lMemStore[0] = complexMul(a[0], w);    w = complexMul(w0, w);
	lMemStore[16] = complexMul(a[1], w);    w = complexMul(w0, w);
	lMemStore[32] = complexMul(a[2], w);    w = complexMul(w0, w);
	lMemStore[48] = complexMul(a[3], w);    w = complexMul(w0, w);
	lMemStore[64] = complexMul(a[4], w);    w = complexMul(w0, w);
	lMemStore[80] = complexMul(a[5], w);    w = complexMul(w0, w);
	lMemStore[96] = complexMul(a[6], w);    w = complexMul(w0, w);
	lMemStore[112] = complexMul(a[7], w);
    ang += ang1;
    w = (float2)(native_cos(ang), native_sin(ang));
	lMemStore[1] = complexMul(a[8], w);    w = complexMul(w0, w);
	lMemStore[17] = complexMul(a[9], w);    w = complexMul(w0, w);
	lMemStore[33] = complexMul(a[10], w);    w = complexMul(w0, w);
	lMemStore[49] = complexMul(a[11], w);    w = complexMul(w0, w);
	lMemStore[65] = complexMul(a[12], w);    w = complexMul(w0, w);
	lMemStore[81] = complexMul(a[13], w);    w = complexMul(w0, w);
	lMemStore[97] = complexMul(a[14], w);    w = complexMul(w0, w);
	lMemStore[113] = complexMul(a[15], w);
	barrier(CLK_LOCAL_MEM_FENCE);
	lMemLoad = sMem + mad24(tid >> 7, 129, tid & 127);
	indexOut += tid;
	out += indexOut;
	out[0] = lMemLoad[0];
	out[128] = lMemLoad[129];
	out[256] = lMemLoad[258];
	out[384] = lMemLoad[387];
	out[512] = lMemLoad[516];
	out[640] = lMemLoad[645];
	out[768] = lMemLoad[774];
	out[896] = lMemLoad[903];
	out[1024] = lMemLoad[1032];
	out[1152] = lMemLoad[1161];
	out[1280] = lMemLoad[1290];
	out[1408] = lMemLoad[1419];
	out[1536] = lMemLoad[1548];
	out[1664] = lMemLoad[1677];
	out[1792] = lMemLoad[1806];
	out[1920] = lMemLoad[1935];
}

/* Run kernel fft_32_1 with global dim = {512*BatchSize}, local dim={64} */
__kernel __attribute__((reqd_work_group_size (64,1,1)))
void ri_fft_32_1(__global float2 *in, __global float *out)
{
    __local float2 sMem[512];
    int i, j, indexIn, indexOut, tid, bNum, xNum;
    float2 w, w0;
    float ang;
    __local float2 *lMemStore, *lMemLoad;

    float2 a[8];
    int lId = get_local_id( 0 );
    int groupId = get_group_id( 0 );
    bNum = groupId & 7;
    xNum = (groupId >> 3) << 12;
    indexIn = mul24(bNum, 16);
    tid = indexIn;
    i = tid >> 7;
    j = tid & 127;
    indexOut = mad24(i, 4096, j+xNum);
    tid = lId;
    i = tid & 15;
    j = tid >> 4;
    indexIn += mad24(j, 128, i+xNum);
    in += indexIn;
    a[0] = in[0];
    a[1] = in[512];
    a[2] = in[1024];
    a[3] = in[1536];
    a[4] = in[2048];
    a[5] = in[2560];
    a[6] = in[3072];
    a[7] = in[3584];
    fftKernel8(a, -1.0f);
    ang = -(M_PI_F/16.0f)*j;
    w0 = (float2)(native_cos(ang), native_sin(ang));
    w = w0;
	lMemStore = sMem + tid;
	lMemStore[0] = a[0];
	lMemStore[64] = complexMul(a[1], w);    w = complexMul(w0, w);
	lMemStore[128] = complexMul(a[2], w);    w = complexMul(w0, w);
	lMemStore[192] = complexMul(a[3], w);    w = complexMul(w0, w);
	lMemStore[256] = complexMul(a[4], w);    w = complexMul(w0, w);
	lMemStore[320] = complexMul(a[5], w);    w = complexMul(w0, w);
	lMemStore[384] = complexMul(a[6], w);    w = complexMul(w0, w);
	lMemStore[448] = complexMul(a[7], w);
    barrier(CLK_LOCAL_MEM_FENCE);
	indexIn = mad24(j, 128, i);
	lMemLoad = sMem + indexIn;
	a[0] = lMemLoad[0];
	a[1] = lMemLoad[16];
	a[2] = lMemLoad[32];
	a[3] = lMemLoad[48];
	a[4] = lMemLoad[64];
	a[5] = lMemLoad[80];
	a[6] = lMemLoad[96];
	a[7] = lMemLoad[112];
    fftKernel4(a, 0, 1, 2, 3, -1.0f);
    fftKernel4(a, 4, 5, 6, 7, -1.0f);
    indexOut += mad24(j, 256, i);
    out += indexOut;
	
	__global float *out_imag = out + 4096;
	outUpdate(out, out_imag, 0, a[0]);
	outUpdate(out, out_imag, 1024, a[1]);
	outUpdate(out, out_imag, 2048, a[2]);
	outUpdate(out, out_imag, 3072, a[3]);
	outUpdate(out, out_imag, 128, a[4]);
	outUpdate(out, out_imag, 1152, a[5]);
	outUpdate(out, out_imag, 2176, a[6]);
	outUpdate(out, out_imag, 3200, a[7]);
}

__kernel
 void doRecon2Dxyz(  __global float *recon, 
				__global float *pSigMat, 
				__global float* pSensor, 
				__global float* pEDT, 
			  int n, float dx, float m2px,
			  float dist, float ua, int4 MN, __global float* pSensMap) 
{
	union _rem {
		int4	i4;
		int		i[4];
	}urem;

	float4 index;
	float4 x, y, z;

	int i = get_global_id(0);
	int j = get_global_id(1);
	if((i < n) && (j < n))
	{
		int rows = MN.x;
		int cols = MN.y;
		int k = n;
		n = hadd(n,-2);
		float acc = convert_float(i - n)*dx;
		float4 deltax = (float4)acc;
		acc = convert_float(n -j)*dx;
		float4 deltay = (float4)acc;
		i = mad24(i,k,j);
		float4 dist4 = (float4)(dist);
		acc = 0.0f;
		int4 rows0123 = (int4) ( 0, rows, 2*rows, 3*rows );
		int4 rows4x = (int4) ( 4*rows );
		float rows_1 = convert_float(rows-1);
		int _cols = (cols & ~3);
		int colsPadded = cols + ((_cols-cols)&3);
		__global float4 *sx = (__global float4 *)pSensor;
		__global float4 *sy = (__global float4 *)(pSensor+colsPadded);
		__global float4 *sz = (__global float4 *)(pSensor+2*colsPadded);
		for(k = 0; k < _cols; k+=4 ) 
		{
			x = *sx++ - deltax;
			y = *sy++ - deltay;
			z = *sz++;
			index = sqrt( x*x+y*y+z*z ) * m2px;
			index = clamp( index - dist4, 0.0f, rows_1);
			int4 rounded = rows0123+convert_int4(index);
			rows0123 += rows4x;
			acc += pSigMat[ rounded.s0 ] ;
			acc += pSigMat[ rounded.s1 ] ;
			acc += pSigMat[ rounded.s2 ] ;
			acc += pSigMat[ rounded.s3 ] ;
		}
		if(k < cols)
		{
			x = *sx - deltax;
			y = *sy - deltay;
			z = *sz;
			index = sqrt( x*x+y*y+z*z ) * m2px;
			index = clamp( index - dist4, 0.0f, rows_1);
			urem.i4 = rows0123 + convert_int4(index);
			for(;k < cols; k++) {
				acc += pSigMat[ urem.i[k&3]];
			}
		}
		recon[i] = -acc * exp(ua*pEDT[i]) * pSensMap[i];
	}
}





__kernel void doRecon3D2dB(  __global float *recon, 
	__global float *pSigMat,
	__global float *pSensor,
	__global float* pEDT, 
	__global float *MipXY, int2 limitsX, int2 limitsY, int2 limitsZ,
	float4 R0, float4 R1, float4 R2, float4 C, float m2px, float ua,
	int4 n, int4 MN)
{
	__global float4 *s_x, *s_y, *s_z;
	float temp[512];
	int rows, i;
	int4 id, rounded, rows0123;
	id.x = get_global_id(0);
	id.y = get_global_id(1);
	if (all(id.xy - n.xy))
	{
		float v, maxv = -MAXFLOAT;
		rows = MN.x;
		id.z = 0;
		id.w = 1;
		rows0123 = (int4)(0, rows, 2 * rows, 3 * rows);
		rows <<= 2;
		s_x = (__global float4 *)pSensor;
		s_y = (__global float4 *)(pSensor + MN.y);
		s_z = (__global float4 *)(pSensor + 2 * MN.y);
		i = mad24(id.y, n.x, id.x);
		recon += i;
		pEDT += i;
		MipXY += i;

		float4 index = convert_float4(id);
		float4 x = (float4)dot(index, R0);
		float4 y = (float4)dot(index, R1);
		for (int j = 0; j < n.z; j++)
		{
			temp[j] = 0;
		}
		for (int k = 0; k < MN.y; k += 4)
		{
			float4 dx = *s_x++ - x;
			float4 dy = *s_y++ - y;
			float4 zs = *s_z++;
			float4 dxy2 = dx*dx + dy*dy;

			for (int iz = 0; iz < n.z; iz++)
			{
				id.z = iz;
				index = convert_float4(id);
				float4 z = (float4)dot(index, R2);
				float4 dz = zs - z;
				index = sqrt(dxy2 + dz*dz) * m2px;
				index = CLAMP_TOF_INDEX4(index);
				rounded = rows0123 + convert_int4(index);

				float acc = pSigMat[rounded.s0] * index.s0;
				acc += pSigMat[rounded.s1] * index.s1;
				acc += pSigMat[rounded.s2] * index.s2;
				acc += pSigMat[rounded.s3] * index.s3;
				temp[iz] -= acc;
			}
			pSigMat += rows;
		}
		i = n.y*n.x;


		for (int j = 0; j < n.z; j++)
		{
			float coef = exp(ua * *pEDT) * IFSAMPLE;
			temp[j] *= coef;
			*recon = temp[j];
			recon += i;
			pEDT += i;
		}
		for (int j = limitsZ.lo; j < limitsZ.hi; j++)
		{
			maxv = fmax(maxv, temp[j]);
		}

		bool cond = (id.x < limitsX.lo) || (id.x >= limitsX.hi) || (id.y < limitsY.lo) || (id.y >= limitsY.hi);
		v = cond ? ATTENUATION_VALUE : 1.0f;
		*MipXY = maxv*v;
	}
}

__kernel void doRecon3D2dF(__global float *recon,
	__global float *pSigMat,
	__global float *pSensor,
	__global float *pEDT,
	__global float *MipXY, int2 limitsX, int2 limitsY, int2 limitsZ,
	float4 R0, float4 R1, float4 R2, float4 C, float m2px, float ua,
	int4 n, int4 MN)
{
	float acc;
	__global float4 *s_x, *s_y, *s_z;
	float temp[512];
	int rows, i;
	int4 id, rounded, rows0123;
	id.x = get_global_id(0);
	id.y = get_global_id(1);
	if (all(id.xy - n.xy))
	{
		float v, maxv = -MAXFLOAT;
		float4 d1, d2, r1;
		float4 c = (float4)C.x;
		rows = MN.x;
		id.z = 0;
		id.w = 1;
		rows0123 = (int4)(0, rows, 2 * rows, 3 * rows);
		rows <<= 2;
		s_x = (__global float4 *)pSensor;
		s_y = (__global float4 *)(pSensor + MN.y);
		s_z = (__global float4 *)(pSensor + 2 * MN.y);
		i = mad24(id.y, n.x, id.x);
		recon += i;
		pEDT += i;
		MipXY += i;
		float4 index = convert_float4(id);
		float4 x = (float4)dot(index, R0);
		float4 y = (float4)dot(index, R1);
		int k = convert_int(C.w - dot(index, R2)) / R2.z;
		int z0 = (MN.w == 0) ? max(k, 0) : 0;
		int z_END = (MN.w == 0) ? n.z : min(k, n.z);
		for (int j = 0; j < n.z; j++)
		{
			temp[j] = 0;
		}
		R2.w -= C.w;
		for (k = 0; k < MN.y; k += 4)
		{
			float4 dx = *s_x++ - x;
			float4 dy = *s_y++ - y;
			float4 z1 = *s_z++ - C.w;
			float4 z1z1 = z1*z1;
			float4 r = sqrt(dx*dx + dy*dy);

			r1 = r*0.5f;
			d1 = sqrt(r1*r1 + z1z1);
			for (int iz = z0; iz < z_END; iz++)
			{
				id.z = iz;
				index = convert_float4(id);
				float4 z = (float4)dot(index, R2);
				if (z.w > 0) continue;
				z *= z;
				d2 = sqrt((r1 - r)*(r1 - r) + z);

				r1 = r*d1 / (d1 + d2*c);
				d1 = sqrt(r1*r1 + z1z1);
				d2 = sqrt(((r1 - r)*(r1 - r) + z));

				r1 = r*d1 / (d1 + d2*c);
				d1 = sqrt(r1*r1 + z1z1);
				d2 = sqrt(((r1 - r)*(r1 - r) + z));

				index = (d1*c + d2) * m2px;
				index = CLAMP_TOF_INDEX4(index);
				rounded = rows0123 + convert_int4(index);

				acc = pSigMat[rounded.s0] * index.s0;
				acc += pSigMat[rounded.s1] * index.s1;
				acc += pSigMat[rounded.s2] * index.s2;
				acc += pSigMat[rounded.s3] * index.s3;
				temp[iz] -= acc;
			}
			pSigMat += rows;
		}
		i = n.y*n.x;
		for (int j = 0; j < z0; j++)
		{
			*recon = 0;
			recon += i;
		}
		pEDT += z0*i;
		for (int j = z0; j < z_END; j++)
		{
			float coef = exp(ua * *pEDT) * IFSAMPLE;
			temp[j] *= coef;
			*recon = temp[j];
			recon += i;
			pEDT += i;
		}
		for (int j = limitsZ.lo; j < limitsZ.hi; j++)
		{
			maxv = fmax(maxv, temp[j]);
		}
		bool cond = (id.x < limitsX.lo) || (id.x >= limitsX.hi) || (id.y < limitsY.lo) || (id.y >= limitsY.hi);
		v = cond ? ATTENUATION_VALUE : 1.0f;
		*MipXY = maxv*v;
	}
}

__kernel void doRecon3D2dC(  __global float *recon, 
	__global float *pSigMat,
	__global float *pSensor,
	__global float* pEDT, 
	__global float *MipXY, int2 limitsX, int2 limitsY, int2 limitsZ,
	float4 R0, float4 R1, float4 R2, float4 C, float m2px, float ua,
	int4 n, int4 MN)
{
	float acc;
	__global float4 *s_x, *s_y, *s_z;
	float temp[512];
	int rows, i;
	int4 id, rounded, rows0123;
	id.x = get_global_id(0);
	id.y = get_global_id(1);
	if (all(id.xy - n.xy))
	{
		float v, maxv = -MAXFLOAT;
		float4 d1, d2;
		float4 c = (float4)(C.x);
		float4 R = (float4)(C.z);
		float4 RR = R*R;
		rows = MN.x;
		id.z = 0;
		id.w = 1;
		rows0123 = (int4)(0, rows, 2 * rows, 3 * rows);
		rows <<= 2;
		s_y = (__global float4 *)pSensor;
		s_x = (__global float4 *)(pSensor + MN.y);
		s_z = (__global float4 *)(pSensor + 2 * MN.y);
		i = mad24(id.y, n.x, id.x);
		recon += i;
		pEDT += i;
		MipXY += i;
		float4 index = convert_float4(id);
		float4 Sx = (float4)dot(index, R0);
		float4 Sy = (float4)dot(index, R1);
		float4 x2py2 = Sx*Sx + Sy*Sy;
		float z = C.w + sqrt(RR.x-x2py2.x);
		int k = convert_int(z - dot(index, R2)) / R2.z;
		int z0 = (MN.w == 0) ? max(k, 0) : 0;
		int z_END = (MN.w == 0) ? n.z : min(k, n.z);
		for (int j = 0; j < n.z; j++)
		{
			temp[j] = 0;
		}
		for (k = 0; k < MN.y; k += 4)
		{
			float4 Dx = *s_x++;
			float4 Dy = *s_y++;
			float4 Dz = *s_z++ - (float4)(C.w);
			float4 dx = Dx - Sx;
			float4 dy = Dy - Sy;
			float4 dx2pdy2 = dx*dx + dy*dy;
			float4 CD2 = Dx*Dx + Dy*Dy + Dz*Dz;
			float4 CD = sqrt(CD2);
			float4 iCDx2 = 0.5f / CD;
			float4 Px = R;
			float4 Py = 0;
			d1 = sqrt(RR + (CD - 2 * Px)*CD);
			for (int iz = z0; iz < z_END; iz++)
			{
				id.z = iz;
				index = convert_float4(id);
				float4 Sz = (float4)(dot(index, R2) - C.w);
				float4 dz = Dz - Sz;
				float4 CS2 = x2py2 + Sz*Sz;
				float4 SD2 = dx2pdy2 + dz*dz;
				float4 xs = (CD2 + CS2 - SD2)*iCDx2;
				float4 xs2 = xs*xs;
				float4 ys2 = CS2 - xs2;
				float4 ys = sqrt(ys2);

				d2 = sqrt(RR + xs2 + ys2 - 2 * (Px*xs + Py*ys));

				Py = Px*ys*d1 / (d1*xs + d2*c*CD);
				Px = sqrt(RR - Py*Py);
				d1 = sqrt(RR + (CD - 2 * Px)*CD);
				d2 = sqrt(RR + xs2 + ys2 - 2 * (Px*xs + Py*ys));

				index = (d1*c + d2) * m2px;
				index = CLAMP_TOF_INDEX4(index);
				rounded = rows0123 + convert_int4(index);

				acc = pSigMat[rounded.s0]  * index.s0;
				acc += pSigMat[rounded.s1] * index.s1;
				acc += pSigMat[rounded.s2] * index.s2;
				acc += pSigMat[rounded.s3] * index.s3;
				temp[iz] -= acc;
			}
			pSigMat += rows;
		}
		i = n.y*n.x;
		for (int j = 0; j < z0; j++)
		{
			*recon = 0;
			recon += i;
		}
		pEDT += z0*i;
		for (int j = z0; j < z_END; j++)
		{
			float coef = exp(ua * *pEDT) * IFSAMPLE;
			temp[j] *= coef;
			*recon = temp[j];
			recon += i;
			pEDT += i;
		}
		for (int j = limitsZ.lo; j < limitsZ.hi; j++)
		{
			maxv = fmax(maxv, temp[j]);
		}
		bool cond = (id.x < limitsX.lo) || (id.x >= limitsX.hi) || (id.y < limitsY.lo) || (id.y >= limitsY.hi);
		v = cond ? ATTENUATION_VALUE : 1.0f;
		*MipXY = maxv*v;
	}
}

__kernel void doRecon3D2dO(__global float *recon,
	__global float *pSigMat,
	__global float *pSensor,
	__global float* pEDT, 
	__global float *MipXY, int2 limitsX, int2 limitsY, int2 limitsZ,
	float4 R0, float4 R1, float4 R2, float4 C, float m2px, float ua,
	int4 n, int4 MN)
{
	int4 id;
	id.x = get_global_id(0);
	id.y = get_global_id(1);
	if (all(id.xy - n.xy))
	{
		float v, maxv = -MAXFLOAT;
		int rows = MN.x;
		int4 rows0123 = (int4)(0, rows, 2 * rows, 3 * rows);
		rows <<= 2;
		int offset = mad24(id.y, n.x, id.x);
		recon += offset;
		pEDT += offset;
		MipXY += offset;
		offset = n.x*n.y;
		id.w = 1;
		for (int iz = 0; iz < n.z; iz++)
		{
			id.z = iz;

			float4 index = convert_float4(id);
			float4 x = (float4)dot(index, R0);
			float4 y = (float4)dot(index, R1);
			float4 z = (float4)dot(index, R2);
			float acc = 0.0f;
			__global float4 *sx = (__global float4 *)pSensor;
			__global float4 *sy = (__global float4 *)(pSensor + MN.y);
			__global float4 *sz = (__global float4 *)(pSensor + 2 * MN.y);
			__global float *s = pSigMat;
			for (int k = 0; k < MN.y; k += 4)
			{
				float4 dx = *sx++ - x;
				float4 dy = *sy++ - y;
				float4 dz = *sz++ - z;
				index = sqrt(dx*dx + dy*dy + dz*dz) * m2px;
				index = CLAMP_TOF_INDEX4(index);
				int4 rounded = rows0123 + convert_int4(index);
				acc += s[rounded.s0] * index.s0;
				acc += s[rounded.s1] * index.s1;
				acc += s[rounded.s2] * index.s2;
				acc += s[rounded.s3] * index.s3;
				s += rows;
			}
			acc *= IFSAMPLE;
			*recon = -acc * exp(ua * *pEDT);
			recon += offset;
			pEDT += offset;
			if ((iz >= limitsZ.lo) && (iz < limitsZ.hi))
			{
				maxv = fmax(maxv, -acc);
			}
		}
		bool cond = (id.x < limitsX.lo) || (id.x >= limitsX.hi) || (id.y < limitsY.lo) || (id.y >= limitsY.hi);
		v = cond ? ATTENUATION_VALUE : 1.0f;
		*MipXY = maxv*v;
	}
}

__kernel void doRecon3D2dP(  __global float *recon, 
	__global float *pSigMat, 
	__global float *pSensor, 
	__global float* pEDT, 
	__global float *MipXY, int2 limitsX, int2 limitsY, int2 limitsZ,
	float4 R0, float4 R1, float4 R2, float4 C, float m2px, float ua,
	int4 n, int4 MN)
{
	int4 id;
	__local float cmat[35];
	id.x = get_global_id(0);
	id.y = get_global_id(1);
	id.z = 0;
	if(all(id.xy-n.xy))
	{
		float v, maxv = -MAXFLOAT;
		int rows = MN.x;
		int offset = mad24(id.y, n.x, id.x);
		recon += offset;
		pEDT += offset;
		MipXY += offset;
		offset = n.x*n.y;
		id.w = 1;
		float c = C.x*C.y; 
		float4 idf = convert_float4(id);
		float x = dot(idf,R0);
		float y = dot(idf,R1);
		float xx = x*x;
		float yy = y*y;
		float xyy = yy*x;
		float yyy = yy*y;
		float xy = x*y;
		for(int iz = 0; iz < n.z; iz++)
		{
			id.z = iz;
			idf = convert_float4(id);
			float z = dot(idf,R2);
			float zz = z*z;
			float xz = x*z;
			float yz = y*z;
			float acc = 0.0f;
			__global float *s = pSigMat;
			for(int k = 0; k < MN.y; k++ ) 
			{
				async_work_group_copy(cmat, pSensor + k*35, 35, 0);
				barrier(CLK_LOCAL_MEM_FENCE);
				float index = FMA(FMA(cmat[2], c, cmat[1]), c, cmat[0]) * c;
				index = FMA(FMA(FMA(cmat[5], c, cmat[4]), c, cmat[3]), x, index);
				index = FMA(FMA(cmat[8],  x, FMA(cmat[7], c, cmat[6])), xx, index);
				index = FMA(FMA(FMA(cmat[11], c, cmat[10]), c, cmat[9]), y, index);
				index = FMA(FMA(cmat[14], x, FMA(cmat[13], c, cmat[12])), xy, index);
				index = FMA(FMA(cmat[17], x, FMA(cmat[16], c, cmat[15])), yy, index);
				index = FMA(cmat[18], yyy, index);
				index = FMA(FMA(FMA(cmat[21], c, cmat[20]), c, cmat[19]), z, index);
				index = FMA(FMA(cmat[24], x, FMA(cmat[23], c, cmat[22])), xz, index);
				index = FMA(FMA(cmat[28], y, FMA(cmat[27], x, FMA(cmat[26], c, cmat[25]))), yz, index);
				index = FMA(FMA(cmat[33], z, FMA(cmat[32], y, FMA(cmat[31], x, FMA(cmat[30], c, cmat[29])))), zz, index + cmat[34]);

				index *= 4e7f;
				index = CLAMP_TOF_INDEX(index);
				int rounded = convert_int(index);
				acc += s[ rounded ]*index ;
				s += rows;
			}
			acc *= IFSAMPLE;
			*recon = -acc * exp(ua * *pEDT);
			recon += offset;
			pEDT += offset;
			if((iz >= limitsZ.lo) && (iz < limitsZ.hi))
			{
				maxv = fmax(maxv, -acc);
			}
		}
		bool cond = (id.x < limitsX.lo) || (id.x >= limitsX.hi) || (id.y < limitsY.lo) || (id.y >= limitsY.hi);
		v = cond ? ATTENUATION_VALUE : 1.0f;
		*MipXY = maxv*v;
	}
}

__kernel void doRecon3D2dP20(  __global float *recon, 
	__global float *pSigMat, 
	__global float *pSensor, 
	__global float* pEDT, 
	__global float *MipXY, int2 limitsX, int2 limitsY, int2 limitsZ,
	float4 R0, float4 R1, float4 R2, float4 C, float m2px, float ua,
	int4 n, int4 MN)
{
	int4 id;
	__local float cmat[20];
	id.x = get_global_id(0);
	id.y = get_global_id(1);
	id.z = 0;
	if(all(id.xy-n.xy))
	{
		float v, maxv = -MAXFLOAT;
		int rows = MN.x;
		int offset = mad24(id.y, n.x, id.x);
		recon += offset;
		pEDT += offset;
		MipXY += offset;
		offset = n.x*n.y;
		id.w = 1;
		float c = C.x*C.y; 
		float4 idf = convert_float4(id);
		float x = dot(idf,R0);
		float y = dot(idf,R1);
		float xx = x*x;
		float yy = y*y;
		float xyy = yy*x;
		float yyy = yy*y;
		float xy = x*y;
		for(int iz = 0; iz < n.z; iz++)
		{
			id.z = iz;
			idf = convert_float4(id);
			float z = dot(idf,R2);
			float zz = z*z;
			float xz = x*z;
			float yz = y*z;
			float acc = 0.0f;
			__global float *s = pSigMat;
			for(int k = 0; k < MN.y; k++ ) 
			{
				async_work_group_copy(cmat, pSensor + k*20, 20, 0);
				barrier(CLK_LOCAL_MEM_FENCE);
				float index = (cmat[0] + ((cmat[1] + cmat[2] * c)*c)) * c;
				index += (cmat[3] + (cmat[4] + cmat[5] * c) *c)*x;
				index += (cmat[6] + cmat[7] * c + cmat[8] * x) *xx;
				index += (cmat[9] + (cmat[10] + cmat[11] * c) *c)*y;
				index += (cmat[12] + cmat[13] * c + cmat[14] * x)*xy;
				index += (cmat[15] + cmat[16] * c + cmat[17] * x) *yy;
				index += cmat[18] * yyy + cmat[19];

				index *= 4e7f;
				index = CLAMP_TOF_INDEX(index);
				int rounded = convert_int(index);
				acc += s[ rounded ] * index;
				s += rows;
			}
			acc *= IFSAMPLE;
			*recon = -acc * exp(ua * *pEDT);
			recon += offset;
			pEDT += offset;
			if((iz >= limitsZ.lo) && (iz < limitsZ.hi))
			{
				maxv = fmax(maxv, -acc);
			}
		}
		bool cond = (id.x < limitsX.lo) || (id.x >= limitsX.hi) || (id.y < limitsY.lo) || (id.y >= limitsY.hi);
		v = cond ? ATTENUATION_VALUE : 1.0f;
		*MipXY = maxv*v;
	}
}

__kernel void doMipXY(__global float *src, __global float *dst,
	int4 n, int2 limitsX, int2 limitsY, int2 limitsZ)
{
	float maxv, v;
	int idx = get_global_id(0);
	int idy = get_global_id(1);
	if ((idx < n.x) && (idy < n.y))
	{
		int frame = get_global_id(2);
		int NxNy = n.x*n.y;
		maxv = -MAXFLOAT;
		int offset = mad24(idy, n.x, idx);
		src += mad24(mad24(frame, n.z, limitsZ.lo), NxNy, offset);
		dst += mad24(NxNy, frame, offset);
		for (int idz = limitsZ.lo; idz < limitsZ.hi; idz++)
		{
			maxv = fmax(maxv, *src);
			src += NxNy;
		}
		bool cond = (idx < limitsX.lo) || (idx >= limitsX.hi) || (idy < limitsY.lo) || (idy >= limitsY.hi);
		v = cond ? ATTENUATION_VALUE : 1.0f;
		dst[0] = maxv*v;
	}
}

__kernel void doMipYZ(__global float *src, __global float *dst,
	int4 n, int2 limitsX, int2 limitsY, int2 limitsZ)
{
	float maxv, v;
	int idy = get_global_id(0);
	int idz = get_global_id(1);
	if ((idy < n.y) && (idz < n.z))
	{
		int frame = get_global_id(2);
		maxv = -MAXFLOAT;
		int offset = mad24(mad24(frame, n.z, idz), n.y, idy);
		src += mad24(offset, n.x, limitsX.lo);
		dst += offset;
		for (int i = limitsX.lo; i < limitsX.hi; i++)
		{
			maxv = fmax(maxv, *src);
			src++;
		}
		bool cond = (idy < limitsY.lo) || (idy >= limitsY.hi) || (idz < limitsZ.lo) || (idz >= limitsZ.hi);
		v = cond ? ATTENUATION_VALUE : 1.0f;
		dst[0] = maxv*v;
	}
}

__kernel void doMipZX(__global float *src, __global float *dst,
	int4 n, int2 limitsX, int2 limitsY, int2 limitsZ)
{
	float maxv, v;
	int idz = get_global_id(0);
	int idx = get_global_id(1);
	if ((idz < n.z) && (idx < n.x))
	{
		int frame = get_global_id(2);
		maxv = -MAXFLOAT;
		src += mad24(mad24(mad24(frame, n.z, idz), n.y, limitsY.lo), n.x, idx);
		dst += mad24(mad24(frame, n.x, idx), n.z, idz);
		for (int idy = limitsY.lo; idy < limitsY.hi; idy++)
		{
			maxv = fmax(maxv, *src);
			src += n.x;
		}
		bool cond = (idx >= limitsX.lo) && (idx < limitsX.hi) && (idz >= limitsZ.lo) && (idz < limitsZ.hi);
		v = cond ? 1.0f : ATTENUATION_VALUE;
		dst[0] = maxv*v;
	}
}



__kernel void doBMipXY(__global unsigned char *src, __global unsigned char *dst,
	int4 n, int2 limitsX, int2 limitsY, int2 limitsZ)
{
	unsigned char maxv, v;
	int idx = get_global_id(0);
	int idy = get_global_id(1);
	if ((idx < n.x) && (idy < n.y))
	{
		int frame = get_global_id(2);
		int NxNy = n.x*n.y;
		maxv = 0;
		int offset = mad24(idy, n.x, idx);
		src += mad24(mad24(frame, n.z, limitsZ.lo), NxNy, offset);
		dst += mad24(NxNy, frame, offset);
		for (int idz = limitsZ.lo; idz < limitsZ.hi; idz++)
		{
			if (*src != 0)
			{
				maxv = 255;
				break;
			}
			src += NxNy;
		}
		bool cond = (idx < limitsX.lo) || (idx >= limitsX.hi) || (idy < limitsY.lo) || (idy >= limitsY.hi);
		v = cond ? 0 : maxv;
		dst[0] = v;
	}
}

__kernel void doBMipYZ(__global unsigned char *src, __global unsigned char *dst,
	int4 n, int2 limitsX, int2 limitsY, int2 limitsZ)
{
	unsigned char maxv, v;
	int idy = get_global_id(0);
	int idz = get_global_id(1);
	if ((idy < n.y) && (idz < n.z))
	{
		int frame = get_global_id(2);
		maxv = 0;
		int offset = mad24(mad24(frame, n.z, idz), n.y, idy);
		src += mad24(offset, n.x, limitsX.lo);
		dst += offset;
		for (int i = limitsX.lo; i < limitsX.hi; i++)
		{
			if (*src != 0)
			{
				maxv = 255;
				break;
			}
			src++;
		}
		bool cond = (idy < limitsY.lo) || (idy >= limitsY.hi) || (idz < limitsZ.lo) || (idz >= limitsZ.hi);
		v = cond ? 0 : maxv;
		dst[0] = v;
	}
}

__kernel void doBMipZX(__global unsigned char *src, __global unsigned char *dst,
	int4 n, int2 limitsX, int2 limitsY, int2 limitsZ)
{
	unsigned char maxv, v;
	int idz = get_global_id(0);
	int idx = get_global_id(1);
	if ((idz < n.z) && (idx < n.x))
	{
		int frame = get_global_id(2);
		maxv = 0;
		src += mad24(mad24(mad24(frame, n.z, idz), n.y, limitsY.lo), n.x, idx);
		dst += mad24(mad24(frame, n.x, idx), n.z, idz);
		for (int idy = limitsY.lo; idy < limitsY.hi; idy++)
		{
			if (*src != 0)
			{
				maxv = 255;
				break;
			}
			src += n.x;
		}
		bool cond = (idx >= limitsX.lo) && (idx < limitsX.hi) && (idz >= limitsZ.lo) && (idz < limitsZ.hi);
		v = cond ? maxv : 0;
		dst[0] = v;
	}
}








#ifdef __DEVICE_IS_CPU__

__attribute__((always_inline))
float CRSSpMV(__global float *V, __global int *Idx, __global float *x, int2 Range) 
{
	float acc = 0.0f;
	for (int j = Range.x; j < Range.y; j++)
	{
		acc = FMA(V[j], x[Idx[j]], acc);
	}
	return acc;
}

__kernel void CRS_SpMV(__global float *V, __global int *Idx, __global int2 *Ptr, __global float *x, __global float *y, int rows, int index) 
{
	int i = get_global_id(0);
	if(i < rows)
	{
		y[i+index*rows] = CRSSpMV(V, Idx, x, Ptr[i]);
	}
}



#else

#define SPMV1LS (3*64/2)

__kernel void CRS_SpMV(__global float *V, __global int *Idx, __global int2 *Ptr, __global float *x, __global float *y, int rows, int index) 
{
	__local int2 lRange;
	__local float lSum[SPMV1LS];
	unsigned i = (unsigned)get_group_id(0);
	int lid = get_local_id(0);
	
	if (lid == 0)
	{
		lRange = Ptr[i];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	int2 Range = lRange;
	
	float acc = 0.0f;
	for (int j = Range.lo + lid; j < Range.hi; j+= get_local_size(0))
	{
		int k = Idx[j];
		acc = FMA(V[j], x[k], acc);
	}

	lSum[lid] = acc;
	barrier(CLK_LOCAL_MEM_FENCE);

	for(int j = get_local_size(0)/2; j > 0; j >>= 1)
	{
		if (lid < j)
		{
			lSum[lid] += lSum[lid + j];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	if (lid == 0) 
	{
		y[i+index*rows] = lSum[0];
	}
}

#endif


__attribute__((always_inline))
void SPMV8MAC(float8 *acc, float4 v, __global float8 *x)
{
	(*acc) = FMA((float8)(v.x), x[0], (*acc));
	(*acc) = FMA((float8)(v.y), x[1], (*acc));
	(*acc) = FMA((float8)(v.z), x[2], (*acc));
	(*acc) = FMA((float8)(v.w), x[3], (*acc));
}


#ifdef __DEVICE_IS_CPU__

__attribute__((always_inline))
void hSpMV_8(__global float4 *A, __global int *Idx, int2 Range, __global float8 *x, 
	__global float8 *y, float8 a, unsigned l, unsigned m, unsigned i, unsigned ii) 
{
	float8 acc1 = a*y[i];
	float8 acc2 = a*y[ii];

	for (unsigned j = (unsigned)Range.x; j < (unsigned)Range.y; j++)
	{
		unsigned k = (unsigned)Idx[j];
		float4 v = A[j];
		SPMV8MAC(&acc1, v, &x[k]);
		k = m-(k/l)*l+k%l;
		SPMV8MAC(&acc2, v, &x[k]);
	}
	y[i] = acc1;
	y[ii] = acc2;
}

__kernel void vector_hSpMV8(__global float4 *A, __global int *Idx, __global int2 *Ptr, __global float8 *x, 
	__global float8 *y, float8 a, unsigned n, unsigned las, unsigned lt) 
{
	unsigned i = (unsigned)get_global_id(0);
	unsigned rows = n*n/2;
	if(i < rows)
	{
		unsigned m = (las-1)*lt;
		int2 Range = Ptr[i];
		unsigned ii = (n-1-i/n)*n+i%n;
		hSpMV_8(A, Idx, Range, x, y, a, lt, m, i, ii);
	}
}

__kernel void vector_hSpMV8T(__global float4 *A, __global int *Idx, __global int2 *Ptr, __global float8 *x, 
	__global float8 *y, float8 a, unsigned n, unsigned las, unsigned lt) 
{
	unsigned i = (unsigned)get_global_id(0);
	int rows = las*lt/2;
	if(i < rows)
	{
		unsigned m = (n-1)*n;
		int2 Range = Ptr[i];
		unsigned ii = (las-1-i/lt)*lt+i%lt;
		hSpMV_8(A, Idx, Range, x, y, a, n, m, i, ii);
	}
}

__kernel void CRS_SpMV8(__global float *V, __global int *Idx, __global int2 *Ptr, __global float8 *x, __global float8 *y, int rows, int index) 
{
	int i = get_global_id(0);
	if(i < rows)
	{
		float8 acc = 0.0f;
		int2 Range = Ptr[i];
		for (int j = Range.x; j < Range.y; j++)
		{
			acc += V[j] * x[Idx[j]];
		}
		y[i+index*rows] = acc;
	}
}
__kernel void Accumulate8(__global float8 *src, __global float8 *dst, int rows)
{
	int i = get_global_id(0);
	if(i < rows)
	{
		float8 acc = 0;

		for(int j = 0; j < 16; j++)
		{
			acc += src[i + j*rows];
		}
		dst[i] = acc;
	}
}

#else

__attribute__((always_inline))
void hSpMV_8(__global float4 *A, __global int *Idx, int2 Range, __global float8 *x, 
	__global float8 *y, float8 a, unsigned l, unsigned m, unsigned i, unsigned ii,
	volatile __local float8* l1Sum,
	volatile __local float8* l2Sum,
	int lid
	) 
{
	float8 acc1 = 0;
	float8 acc2 = 0;
	int ll = l*2;
	for (int j = Range.lo + lid; j < Range.hi; j+= get_local_size(0))
	{
		unsigned k = (unsigned)Idx[j];
		float4 v = A[j];
		SPMV8MAC(&acc1, v, &x[k]);
		k += m-(k/l)*ll;
		SPMV8MAC(&acc2, v, &x[k]);
	}
	l1Sum[lid] = acc1;
	l2Sum[lid] = acc2;
	barrier(CLK_LOCAL_MEM_FENCE);

	for(int j = get_local_size(0)/2; j > 0; j >>= 1)
	{
		if (lid < j)
		{
			l1Sum[lid] += l1Sum[lid + j];
			l2Sum[lid] += l2Sum[lid + j];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	
	// write result for this block to global mem
	if (lid == 0) 
	{
		y[i] = FMA(a, y[i], l1Sum[0]);
		y[ii] = FMA(a, y[ii], l2Sum[0]);
	}
}

#define SPMV8LS (3*64/2)
__kernel void vector_hSpMV8(__global float4 *A, __global int *Idx, __global int2 *Ptr, __global float8 *x, 
	__global float8 *y, float8 a, unsigned n, unsigned las, unsigned lt) 
{
	__local int2 lRange;
	__local float8 l1Sum[SPMV8LS];
	__local float8 l2Sum[SPMV8LS];
	unsigned i = (unsigned)get_group_id(0);
	int lid = get_local_id(0);
	
	if (lid == 0)
	{
		lRange = Ptr[i];
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	unsigned m = (las-1)*lt;
	unsigned ii = (n-1-i/n)*n+i%n;
	int2 Range = lRange;
	hSpMV_8(A, Idx, Range, x, y, a, lt, m, i, ii, l1Sum, l2Sum, lid);
}

__kernel void vector_hSpMV8T(__global float4 *A, __global int *Idx, __global int2 *Ptr, __global float8 *x, 
	__global float8 *y, float8 a, unsigned n, unsigned las, unsigned lt) 
{
	__local int2 lRange;
	__local float8 l1Sum[SPMV8LS];
	__local float8 l2Sum[SPMV8LS];
	unsigned i = (unsigned)get_group_id(0);
	int lid = get_local_id(0);

	if (lid == 0)
	{
		lRange = Ptr[i];
	}
	barrier(CLK_LOCAL_MEM_FENCE);

	unsigned m = (n-1)*n;
	unsigned ii = (las-1-i/lt)*lt+i%lt;
	int2 Range = lRange;
	hSpMV_8(A, Idx, Range, x, y, a, n, m, i, ii, l1Sum, l2Sum, lid);
}


__kernel void CRS_SpMV8(__global float *V, __global int *Idx, __global int2 *Ptr, __global float4 *x, __global float4 *y, int rows, int index) 
{
	__local int2 lRange;
	__local float4 l1Sum[SPMV8LS];
	__local float4 l2Sum[SPMV8LS];
	unsigned i = (unsigned)get_group_id(0);
	int lid = get_local_id(0);

	if (lid == 0)
	{
		lRange = Ptr[i];
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	float4 acc1, acc2;
	acc1 = acc2 = 0.0f;
	int2 Range = lRange;
	for (int j = Range.x + lid; j < Range.y; j+=get_local_size(0))
	{
		int ind = 2*Idx[j];
		acc1 = FMA((float4)V[j], x[ind], acc1);
		acc2 = FMA((float4)V[j], x[ind+1], acc2);
	}
	l1Sum[lid] = acc1;
	l2Sum[lid] = acc2;
	barrier(CLK_LOCAL_MEM_FENCE);

	for(int j = get_local_size(0)/2; j > 0; j >>= 1)
	{
		if (lid < j)
		{
			l1Sum[lid] += l1Sum[lid + j];
			l2Sum[lid] += l2Sum[lid + j];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	
	if (lid == 0) 
	{
		int ind = 2*(i+index*rows);
		y[ind] = l1Sum[0];
		y[ind+1] = l2Sum[0];
	}
}

__kernel void Accumulate8(__global float4 *src, __global float4 *dst, int rows)
{
	int i = get_global_id(0);
	if(i < rows)
	{
		float4 acc = 0;

		for(int j = 0; j < 16; j++)
		{
			acc += src[(i + j*rows)*2];
		}
		dst[i*2] = acc;
		dst++;src++;
		acc = 0;
		for(int j = 0; j < 16; j++)
		{
			acc += src[(i + j*rows)*2];
		}
		dst[i*2] = acc;
	}
}

#endif

__constant float2 DB6[12] = {
(float2)( 0.11154074335008017000f,-0.00107730108499557990f),
(float2)( 0.49462389039838539000f,-0.00477725751101065140f),
(float2)( 0.75113390802157753000f, 0.00055384220099380160f),
(float2)( 0.31525035170924320000f, 0.03158203931803115600f),
(float2)(-0.22626469396516913000f, 0.02752286553001628800f),
(float2)(-0.12976686756709563000f,-0.09750160558707936200f),
(float2)( 0.09750160558707936200f,-0.12976686756709563000f),
(float2)( 0.02752286553001628800f, 0.22626469396516913000f),
(float2)(-0.03158203931803115600f, 0.31525035170924320000f),
(float2)( 0.00055384220099380160f,-0.75113390802157753000f),
(float2)( 0.00477725751101065140f, 0.49462389039838539000f),
(float2)(-0.00107730108499557990f,-0.11154074335008017000f)

};

#define ROWEXT (Rows+11)
#define COLEXT (Cols+11)
#define ROWEXTH (ROWEXT/2)
#define COLEXTH (COLEXT/2)

__kernel 
void waveletX_dec2D_ROWS (__global float *source, int Rows, int Cols, __global float *target) 
{ 
	int r = get_global_id (0);
	float2 sum;
	int cj,ex, offset;
	if (r < Rows) 
	{ 
		ex = COLEXTH*Rows;
		offset = r;
		__global float *src = source+r;
		for (int c = 0; c < COLEXT-2; c+=2) 
		{ 
			sum = 0.0f; 
			for (int j = 0; j < 12; j++) 
			{ 
				cj = c+j-10;	
				cj = (cj >= 0) ? cj : ~cj; 
				cj = (cj < Cols) ? cj : 2*Cols-cj-1; 
				sum = FMA(src[cj*Rows], DB6[j], sum); 
			} 
			target[offset] = sum.x;
			target[offset+ex] = sum.y;
			offset += Rows;
		}
	} 
} 

__kernel 
void waveletX_dec2D_COLS (__global float *source, int Rows, int Cols, __global float *target) 
{ 
	int c = get_global_id (0);
	float2 sum;
	int rj,offset,ex;
	if (c < 2*COLEXTH) 
	{ 
		offset = (c % COLEXTH)*ROWEXTH;
		rj = 2*(c / COLEXTH);
		ex = COLEXTH*ROWEXTH;
		offset += rj*ex;
		__global float *src = source+c*Rows;
		for (int r = 0; r < ROWEXT-2; r+=2) 
		{ 
			sum = 0.0f; 
			for (int j = 0; j < 12; j++) 
			{ 
				rj = r+j-10;	
				rj = (rj >= 0) ? rj : ~rj; 
				rj = (rj < Rows) ? rj : 2*Rows-rj-1; 
				sum = FMA(src[rj], DB6[j], sum); 
			} 
			target[offset] = sum.x;
			target[offset+ex] = sum.y;
			offset++;
		} 
	} 
} 


#ifdef __V_INTEL__
__kernel 
void waveletX_dec2D8_ROWS (__global float8 *source, int Rows, int Cols, __global float8 *target) 
{ 
	int r = get_global_id (0);
	float8 sumx, sumy;
	int cj,ex, offset;
	if (r < Rows) 
	{ 
		ex = COLEXTH*Rows;
		offset = r;
		__global float8 *src = source+r;
		for (int c = 0; c < COLEXT-2; c+=2) 
		{ 
			sumx = sumy = 0.0f; 
			for (int j = 0; j < 12; j++) 
			{ 
				cj = c+j-10;	
				cj = (cj >= 0) ? cj : ~cj; 
				cj = (cj < Cols) ? cj : 2*Cols-cj-1; 
				float8 v = src[cj*Rows];
				sumx = FMA(v, (float8)DB6[j].x, sumx); 
				sumy = FMA(v, (float8)DB6[j].y, sumy); 
			} 
			target[offset] = sumx;
			target[offset+ex] = sumy;
			offset += Rows;
		}
	} 
} 

__kernel 
void waveletX_dec2D8_COLS (__global float8 *source, int Rows, int Cols, __global float8 *target) 
{ 
	int c = get_global_id (0);
	float8 sumx,sumy;
	int rj,offset,ex;
	if (c < 2*COLEXTH) 
	{ 
		offset = (c % COLEXTH)*ROWEXTH;
		rj = 2*(c / COLEXTH);
		ex = COLEXTH*ROWEXTH;
		offset += rj*ex;
		__global float8 *src = source+c*Rows;
		for (int r = 0; r < ROWEXT-2; r+=2) 
		{ 
			sumx = sumy = 0.0f; 
			for (int j = 0; j < 12; j++) 
			{ 
				rj = r+j-10;	
				rj = (rj >= 0) ? rj : ~rj; 
				rj = (rj < Rows) ? rj : 2*Rows-rj-1; 
				float8 v = src[rj];
				sumx = FMA(v, (float8)DB6[j].x, sumx); 
				sumy = FMA(v, (float8)DB6[j].y, sumy); 
			} 
			target[offset] = sumx;
			target[offset+ex] = sumy;
			offset++;
		} 
	} 
} 
#else

__kernel 
void waveletX_dec2D8_ROWS (__global float4 *source, int Rows, int Cols, __global float4 *target) 
{ 
	int r = get_global_id (0);
	float4 sumx, sumy, sumz, sumw;
	int cj,ex, offset;
	if (r < Rows) 
	{ 
		ex = COLEXTH*Rows;
		offset = r;
		__global float4 *src = source+r*2;
		for (int c = 0; c < COLEXT-2; c+=2) 
		{ 
			sumx = sumy = sumz  = sumw = 0.0f; 
			for (int j = 0; j < 12; j++) 
			{ 
				cj = c+j-10;	
				cj = (cj >= 0) ? cj : ~cj; 
				cj = (cj < Cols) ? cj : 2*Cols-cj-1; 
				float4 v = src[cj*Rows*2];
				sumx = FMA(v, (float4)DB6[j].x, sumx); 
				sumy = FMA(v, (float4)DB6[j].y, sumy); 
				v = src[cj*Rows*2+1];
				sumz = FMA(v, (float4)DB6[j].x, sumz); 
				sumw = FMA(v, (float4)DB6[j].y, sumw); 
			} 
			target[offset*2] = sumx;
			target[offset*2+1] = sumz;
			target[(offset+ex)*2] = sumy;
			target[(offset+ex)*2+1] = sumw;
			offset += Rows;
		}
	} 
} 

__kernel 
void waveletX_dec2D8_COLS (__global float4 *source, int Rows, int Cols, __global float4 *target) 
{ 
	int c = get_global_id (0);
	float4 sumx, sumy, sumz, sumw;
	int rj,offset,ex;
	if (c < 2*COLEXTH) 
	{ 
		offset = (c % COLEXTH)*ROWEXTH;
		rj = 2*(c / COLEXTH);
		ex = COLEXTH*ROWEXTH;
		offset += rj*ex;
		__global float4 *src = source+c*Rows*2;
		for (int r = 0; r < ROWEXT-2; r+=2) 
		{ 
			sumx = sumy = sumz  = sumw = 0.0f; 
			for (int j = 0; j < 12; j++) 
			{ 
				rj = r+j-10;	
				rj = (rj >= 0) ? rj : ~rj; 
				rj = (rj < Rows) ? rj : 2*Rows-rj-1; 
				float4 v = src[rj*2];
				sumx = FMA(v, (float4)DB6[j].x, sumx); 
				sumy = FMA(v, (float4)DB6[j].y, sumy); 
				v = src[rj*2+1];
				sumz = FMA(v, (float4)DB6[j].x, sumz); 
				sumw = FMA(v, (float4)DB6[j].y, sumw); 
			} 
			target[offset*2] = sumx;
			target[offset*2+1] = sumz;
			target[(offset+ex)*2] = sumy;
			target[(offset+ex)*2+1] = sumw;
			offset++;
		} 
	} 
} 
#endif


__attribute__((always_inline))
float2 rc_range(float4 a, float4 b)
{
	float4 aa,c,d;
	aa.w = 1.0f;
	aa.x = (fabs(a.x) < FLT_EPSILON) ? ((a.x >= 0) ? FLT_EPSILON : -FLT_EPSILON) : a.x;
	aa.y = (fabs(a.y) < FLT_EPSILON) ? ((a.y >= 0) ? FLT_EPSILON : -FLT_EPSILON) : a.y;
	aa.z = (fabs(a.z) < FLT_EPSILON) ? ((a.z >= 0) ? FLT_EPSILON : -FLT_EPSILON) : a.z;
	if(aa.x < 0) {c.x = -b.x; d.x = 1.0f-b.x; } else { c.x = 1.0f-b.x; d.x = -b.x; }
	if(aa.y < 0) {c.y = -b.y; d.y = 1.0f-b.y; } else { c.y = 1.0f-b.y; d.y = -b.y; }
	if(aa.z < 0) {c.z = -b.z; d.z = 1.0f-b.z; } else { c.z = 1.0f-b.z; d.z = -b.z; }
	c /= aa;
	d /= aa;
	return (float2)(fmax(d.x, fmax(d.y,d.z)), fmin(c.x, fmin(c.y,c.z)));
}



#define FVOLUME_COORDINATES(X,Y,Z) (float4)(1.0f-Y, X,1.0f-Z,0)
#define IVOLUME_COORDINATES(X,Y,Z) (int4)(dims.x-1-Y, X,dims.y-1-Z,0)

#define RAYCASTHEAD \
	int2 coord = (int2)(get_global_id(0), get_global_id(1));\
	if(coord.x >= dims.x || coord.y >= dims.x) return; \
	float x = native_divide((float)(get_global_id(0)), (float)(dims.x-1));\
	float y = native_divide((float)(get_global_id(1)), (float)(dims.x-1));\
	float4 nearp = (float4)(x, y, geo.y, 1.0f);\
	float4 pos;\
	pos.x = dot(nearp, iViewMatX);\
	pos.y = dot(nearp, iViewMatY);\
	pos.z = dot(nearp, iViewMatZ);\
	pos.w = 1.0f;\
	float4 ray = fast_normalize(pos - camera);\
	float threshold = geo.z;\
	float2 range = rc_range(ray, pos)

#define RC_SAMPLING_LOOP_B2F for(float t = range.hi; t >= range.lo; t -= geo.x)
#define RC_SAMPLING_LOOP_F2B for(float t = range.lo; t <= range.hi; t += geo.x)


#define TABLE_LOOKUP(_op_) \
		fIndex += 0.5f;\
		fIndex *= _op_;\
		float a = read_imagef(alphaTable, iSamplerNEL, (float2)(I,fIndex)).x;\
		float4 rgb = read_imagef(colorTable, iSamplerNEL, (float2)(I,fIndex))



#define RC_BASIC_ARGS __read_only image3d_t input, __write_only image2d_t output, float4 iViewMatX, float4 iViewMatY, float4 iViewMatZ, float4 camera, int2 dims, float4 geo, __read_only image2d_t colorTable, float4 cw

__kernel void RC_MIP_MC(RC_BASIC_ARGS)
{
	float2 v[4];
	RAYCASTHEAD;

	float color = 0.0f;
	float findex = 0;
	int ind = 0;
	RC_SAMPLING_LOOP_B2F
	{
		float4 tpos = pos + t*ray;
		float4 I4 = read_imagef(input, iSamplerNCL, FVOLUME_COORDINATES(tpos.x, tpos.y, tpos.z));

		v[0] = (float2)(I4.s0 * cw.s0, I4.s0);
		v[1] = (float2)(I4.s1 * cw.s1, I4.s1);
		v[2] = (float2)(I4.s2 * cw.s2, I4.s2);
		v[3] = (float2)(I4.s3 * cw.s3, I4.s3);
		int index = 0;
		index = (v[index].x < v[1].x) ? 1 : index;
		index = (v[index].x < v[2].x) ? 2 : index;
		index = (v[index].x < v[3].x) ? 3 : index;

		float I = v[index].y;

		if (I < threshold) continue;

		if (I>color) { color = I; ind = index; }
	}
	color = clamp(color, 0.0f, 1.0f);
	if (color < threshold)
	{
		color = 0.0f;
		ind = 0;
	}
	findex = (float)ind;
	findex += 0.5f;
	findex /= ((float)get_image_height(colorTable));
	float4 rgb = read_imagef(colorTable, iSamplerNEL, (float2)(color, findex));
	rgb.w = 1.0f;
	write_imagef(output, coord, rgb);
}

__attribute__((always_inline))
float4 Blend_Normal(float4 A, float4 B) { return B; }

__attribute__((always_inline))
float4 Blend_Lighten(float4 A, float4 B) { return max(A, B); }


#define RC_COLOR_type(Method) \
__kernel void RC_COLOR_##Method(RC_BASIC_ARGS, __read_only image2d_t alphaTable)\
{\
	float4 v[4];\
	RAYCASTHEAD;\
	int components = get_image_height(alphaTable);\
	float4 color = 0.0f;\
	float Iacc = 0.0f;\
	float icomps = 1.0f / ((float)components);\
	RC_SAMPLING_LOOP_B2F\
	{\
		float4 tpos = pos + t*ray;\
		float4 I4 = read_imagef(input, iSamplerNCN, FVOLUME_COORDINATES(tpos.x, tpos.y, tpos.z));\
		float4 I4a = read_imagef(input, iSamplerNCL, FVOLUME_COORDINATES(tpos.x, tpos.y, tpos.z));\
		\
		v[0] = (float4)(I4.s0 * cw.s0, I4.s0, I4a.s0, 0);\
		v[1] = (float4)(I4.s1 * cw.s1, I4.s1, I4a.s1, 0);\
		v[2] = (float4)(I4.s2 * cw.s2, I4.s2, I4a.s2, 0);\
		v[3] = (float4)(I4.s3 * cw.s3, I4.s3, I4a.s3, 0);\
		int index = 0;\
		index = (v[index].x < v[1].x) ? 1 : index;\
		index = (v[index].x < v[2].x) ? 2 : index;\
		index = (v[index].x < v[3].x) ? 3 : index;\
		\
		float I = v[index].z;\
		\
		if (I < threshold) continue;\
		\
		float fIndex = (float)(index)+0.5f;\
		fIndex *= icomps;\
		float a = read_imagef(alphaTable, iSamplerNEL, (float2)(I, fIndex)).x;\
		float4 rgb = read_imagef(colorTable, iSamplerNEL, (float2)(I, fIndex));\
		rgb = Blend_##Method(color,rgb);\
		color = mix(color, rgb, a);\
		Iacc = mix(Iacc, I, a);\
	}\
	color.w = 1.0f;\
	write_imagef(output, coord, color);\
}

RC_COLOR_type(Normal)

RC_COLOR_type(Lighten)



#define zStart zConf.x
#define zStop zConf.y
#define zDelta zConf.z

#define SHEARCOREHEAD \
	int xyz[20];\
	int2 coord = (int2)(get_global_id(0), get_global_id(1));\
	int2 IbufferDims = get_image_dim(Ibuffer);\
	if(all(coord-IbufferDims) == 0) return;\
	IbufferDims -= 2;\
	int2 coordnb = coord-1;\
	int2 test = coordnb * (coordnb - IbufferDims);\
	int *pxyz = xyz + zConf.w

#define SW_SAMPLING_LOOP if((test.x <=0) && (test.y <=0)) for(int z = zStart; z != zStop; z += zDelta)

#define SHEARCOREBODY \
	float4 v[4];\
	float4 weights, Ix, Iy, Iz, Iw, I4N, I4L; \
	if (any(pend[z] - coord) == 1) continue;\
	int4 limit = pstart[z];\
	if (any(coord - limit.xy) == 1) continue;\
	limit.zw += coord;\
	xyz[8] = xyz[9] = xyz[10] = xyz[11] = z;\
	xyz[0] = xyz[1] = xyz[12] = xyz[13] = limit.z;\
	xyz[4] = xyz[6] = xyz[16] = xyz[18] = limit.w;\
	limit.zw++;\
	xyz[2] = xyz[3] = xyz[14] = xyz[15] = limit.z;\
	xyz[5] = xyz[7] = xyz[17] = xyz[19] = limit.w;\
	v[0] = read_imagef(vol, iSamplerUCN, IVOLUME_COORDINATES(pxyz[0], pxyz[4], pxyz[8]));\
	v[1] = read_imagef(vol, iSamplerUCN, IVOLUME_COORDINATES(pxyz[1], pxyz[5], pxyz[9]));\
	v[2] = read_imagef(vol, iSamplerUCN, IVOLUME_COORDINATES(pxyz[2], pxyz[6], pxyz[10]));\
	v[3] = read_imagef(vol, iSamplerUCN, IVOLUME_COORDINATES(pxyz[3], pxyz[7], pxyz[11]));\
	Ix = (float4)(v[0].x, v[1].x, v[2].x, v[3].x);\
	Iy = (float4)(v[0].y, v[1].y, v[2].y, v[3].y);\
	Iz = (float4)(v[0].z, v[1].z, v[2].z, v[3].z);\
	Iw = (float4)(v[0].w, v[1].w, v[2].w, v[3].w);\
	weights = percL[z];\
	I4L.s0 = dot(Ix, weights);\
	I4L.s1 = dot(Iy, weights);\
	I4L.s2 = dot(Iz, weights);\
	I4L.s3 = dot(Iw, weights);\
	weights = percN[z];\
	I4N.s0 = dot(Ix, weights);\
	I4N.s1 = dot(Iy, weights);\
	I4N.s2 = dot(Iz, weights);\
	I4N.s3 = dot(Iw, weights);\
	v[0] = (float4)(I4N.s0 * cw.s0, I4L.s0,I4N.s0,0);\
	v[1] = (float4)(I4N.s1 * cw.s1, I4L.s1,I4N.s1,0);\
	v[2] = (float4)(I4N.s2 * cw.s2, I4L.s2,I4N.s2,0);\
	v[3] = (float4)(I4N.s3 * cw.s3, I4L.s3,I4N.s3,0);\
	int ind = 0;\
	ind = (v[0].x < v[1].x) ? 1 : 0;\
	ind = (v[ind].x < v[2].x) ? 2 : ind;\
	ind = (v[ind].x < v[3].x) ? 3 : ind;\
	float I = v[ind].y;\
	float fIndex = (float)(ind);\
	if (I < thrsh) continue


	
#define SW_BASIC_ARGS \
	__read_only image3d_t vol, \
	__write_only image2d_t Ibuffer, \
	__global float4 *percL, \
	__global float4 *percN, \
	__global int4 *pstart, \
	__global int2 *pend, \
	int4 zConf, \
	float thrsh, \
	int2 dims, \
	__read_only image2d_t colorTable, \
	float4 cw

__kernel void ShearWarp_PrepareVars(
	__global float4 *percL, 
	__global int4 *pstart, 
	__global int2 *pend,
	int4 __iSIZE, float4 iMshear,
	__global float4 *percN)
{
	int z = get_global_id(0);
	if(z < __iSIZE.z)
	{
		float zd = (float)(z);
		float d = iMshear.x + iMshear.z * zd;
		float dfloor = floor(d);
		int xi = (int)(dfloor);
		float xCom = d - dfloor;
		d = iMshear.y + iMshear.w * zd;
		dfloor = floor(d);
		int yi = (int)(dfloor);
		float yCom = d - dfloor;
		int pxend = __iSIZE.x - xi;
		int pyend = __iSIZE.y - yi;
		float4 aPerc = (float4)((1.0f - xCom) * (1.0f - yCom), (1.0f - xCom) * yCom, xCom * (1.0f - yCom), xCom * yCom);
		int4 aStart = (int4)((xi > 0) ? 0 : -xi, (yi > 0) ? 0 : -yi, xi, yi);
		int2 aEnd = (int2)((pxend > __iSIZE.w) ? __iSIZE.w - 2 : pxend - 2, (pyend > __iSIZE.w) ? __iSIZE.w - 2 : pyend - 2);
		percL[z] = aPerc;
		pstart[z] = aStart;
		pend[z] = aEnd;

		xCom = round(xCom);
		yCom = round(yCom);
		aPerc = (float4)((1.0f - xCom) * (1.0f - yCom), (1.0f - xCom) * yCom, xCom * (1.0f - yCom), xCom * yCom);
		percN[z] = aPerc;
	}
}

__kernel void ShearCore_MIP_MC(SW_BASIC_ARGS)
{
	SHEARCOREHEAD;

	float mi = 0.0f;
	float index = 0;
	SW_SAMPLING_LOOP
	{
		SHEARCOREBODY;

		if (I>mi) { mi = I; index = fIndex; }
	}
	mi = clamp(mi, 0.0f, 1.0f);
	index += 0.5f;
	index /= ((float)get_image_height(colorTable));
	float4 rgb = read_imagef(colorTable, iSamplerNEL, (float2)(mi, index));
	rgb.w = 1.0f;
	write_imagef(Ibuffer, coord, rgb);
}

__kernel void ShearCore_COLOR_MC(
	SW_BASIC_ARGS,
	__read_only image2d_t alphaTable)
{
	SHEARCOREHEAD;

	float4 cv = 0;
	float Iacc = 0.0f;
	float icomps = 1.0f / ((float)get_image_height(alphaTable));
	SW_SAMPLING_LOOP
	{
		SHEARCOREBODY;

		fIndex += 0.5f; 
		fIndex *= icomps;
		float a = read_imagef(alphaTable, iSamplerNEL, (float2)(I, fIndex)).x; 
		float4 rgb = read_imagef(colorTable, iSamplerNEL, (float2)(I, fIndex));

		cv = mix(cv, rgb, a);
		Iacc = mix(Iacc, I, a);
	}
	cv.w = 1.0f;
	write_imagef(Ibuffer, coord, cv);
}


__kernel void WarpCore(__write_only image2d_t output, __read_only image2d_t Ibuffer, float4 cs, float4 m0314)
{
	int2 p = (int2)(get_global_id(0), get_global_id(1));
	int2 Iout_size = get_image_dim(output);
	if (all(p - Iout_size) == 1)
	{
		float2 xydt;
		float2 xyd = convert_float2(2 * p - Iout_size) / 2.0f;
		xydt.x = dot(xyd, m0314.xy);
		xydt.y = dot(xyd, m0314.zw);
		xydt += cs.xy;
		xydt *= cs.zw;

		float4 data = read_imagef(Ibuffer, iSamplerNEL, xydt);
		write_imagef(output, p, data);
	}
}

__kernel void kf_MinMax(__global float2  *gminmax, __global float  *src, int2 count)
{
#ifdef __DEVICE_IS_CPU__

	int idx  = get_global_id(0) * count.x;
	int iend = idx + count.x;
	if(iend > count.y) iend = count.y;
	float2 pminmax = (float2)src[idx];
	idx++;
	for (; idx < iend; idx++) 

#else

	int idx = get_global_id(0);
	float2 pminmax = (float2)src[idx];
	for (int i = 0; (i < count.x) && (idx < count.y); i++)
	{
		barrier(CLK_GLOBAL_MEM_FENCE);
		float d = src[idx];
		if (pminmax.x > d) { pminmax.x = d; }
		if (pminmax.y < d) { pminmax.y = d; }
		idx += get_local_size(0);
	}
	if (idx < count.y)

#endif

	{
		float d = src[idx];
		if (pminmax.x > d) { pminmax.x = d; }
		if (pminmax.y < d) { pminmax.y = d; }
	}
	int lid = get_local_id(0);
	barrier(CLK_GLOBAL_MEM_FENCE);
	if (lid == 0) gminmax[0] = pminmax;
	for (int n = 1; n < get_local_size(0); n++)
	{
		barrier(CLK_GLOBAL_MEM_FENCE);
		if (lid == n)
		{
			if (pminmax.x < gminmax[0].x) { gminmax[0].x = pminmax.x; }
			if (pminmax.y > gminmax[0].y) { gminmax[0].y = pminmax.y; }
		}
	}
	barrier(CLK_GLOBAL_MEM_FENCE);
}

__kernel void kf_Normalize(__global float* d, __global float* s, float a, float b, int size)
{
	int i = get_global_id(0);
	if(i < size)
	{
		float v = FMA(s[i],a,b);
		d[i] = (v < 0) ? 0 : v;
	}
}

__attribute__((always_inline))
bool isEdge3D(int4 dims, int i)
{
	if (dims.w == 1) return false;
	int xy = dims.x * dims.y;
	int z = i / xy;
	int y = (i - z*xy) / dims.x;
	int x = (i - z*xy - y*dims.x);
	int xe = x*((x + 1) % dims.x) == 0 ? 0 : 1;
	int ye = y*((y + 1) % dims.y) == 0 ? 0 : 1;
	int ze = z*((z + 1) % dims.z) == 0 ? 0 : 1;
	return (xe + ye)*(xe + ze)*(ye + ze) == 0;
}


__kernel void k_ComposeVolume(__global uchar4 *volume, __global float *volumes, float4 a, float4 b, int4 dims, int components)
{
	int i = get_global_id(0);
	int total = dims.x * dims.y * dims.z;
	if (i < total)
	{
		bool edge = isEdge3D(dims, i);
		if (!edge)
		{
			__global float *cv = volumes + i*components;
			float4 v = 0;

			int c = components;
			if (0 < c--)
			{
				v.x = FMA(*cv, a.s0, b.s0);
			}
			if (0 < c--)
			{
				cv++;
				v.y = FMA(*cv, a.s1, b.s1);
			}
			if (0 < c--)
			{
				cv++;
				v.z = FMA(*cv, a.s2, b.s2);
			}
			if (0 < c--)
			{
				cv++;
				v.w = FMA(*cv, a.s3, b.s3);
			}
			volume[i] = convert_uchar4_sat(v*255.0f);
		}
		else
		{
			volume[i] = (uchar4)(255, 0, 0, 0);
		}
	}
}


__kernel void kf_Var(__global float2 *g, __global float  *src, int2 count)
{
	float2 s = 0;

#ifdef __DEVICE_IS_CPU__

	int idx = get_global_id(0) * count.x;
	int iend = idx + count.x;
	if (iend > count.y) iend = count.y;
	for (; idx < iend; idx++)

#else

	int idx = get_global_id(0);
	for (int i = 0; (i < count.x) && (idx < count.y); i++)
	{
		float x = src[idx];
		barrier(CLK_GLOBAL_MEM_FENCE);
		s += (float2)(x, x*x);
		idx += get_local_size(0);
	}
	if (idx < count.y)

#endif

	{
		float x = src[idx];
		s += (float2)(x, x*x);
	}
	int lid = get_local_id(0);
	barrier(CLK_GLOBAL_MEM_FENCE);
	if (lid == 0) g[0] = s;
	for (int n = 1; n < get_local_size(0); n++)
	{
		barrier(CLK_GLOBAL_MEM_FENCE);
		if (lid == n) { g[0] += s; }
	}
	barrier(CLK_GLOBAL_MEM_FENCE);
}

__kernel void kf_CoVar(__global float *g, __global float  *srcx, __global float  *srcy, int2 count)
{
	float3 s = 0;

#ifdef __DEVICE_IS_CPU__

	int idx = get_global_id(0) * count.x;
	int iend = idx + count.x;
	if (iend > count.y) iend = count.y;
	for (; idx < iend; idx++)

#else

	int idx = get_global_id(0);
	for (int i = 0; (i < count.x) && (idx < count.y); i++)
	{
		float x = srcx[idx];
		float y = srcy[idx];
		barrier(CLK_GLOBAL_MEM_FENCE);
		s += (float3)(x, y, x*y);
		idx += get_local_size(0);
	}
	if (idx < count.y)

#endif

	{
		float x = srcx[idx];
		float y = srcy[idx];
		s += (float3)(x, y, x*y);
	}
	int lid = get_local_id(0);
	barrier(CLK_GLOBAL_MEM_FENCE);
	if (lid == 0)
	{
		g[0] = s.s0;
		g[1] = s.s1;
		g[2] = s.s2;
	}
	for (int n = 1; n < get_local_size(0); n++)
	{
		barrier(CLK_GLOBAL_MEM_FENCE);
		if (lid == n) 
		{ 
			g[0] += s.s0;
			g[1] += s.s1;
			g[2] += s.s2;
		}
	}
	barrier(CLK_GLOBAL_MEM_FENCE);
}

__kernel void kf_CorrCoefHelper(__global float  *g, __global float  *srcx, __global float  *srcy, int2 count)
{
	float3 sy = 0;

#ifdef __DEVICE_IS_CPU__

	int idx = get_global_id(0) * count.x;
	int iend = idx + count.x;
	if (iend > count.y) iend = count.y;
	for (; idx < iend; idx++)

#else

	int idx = get_global_id(0);
	for (int i = 0; (i < count.x) && (idx < count.y); i++)
	{
		float x = srcx[idx];
		float y = srcy[idx];
		barrier(CLK_GLOBAL_MEM_FENCE);
		sy += (float3)(y, y*y, x*y);
		idx += get_local_size(0);
	}
	if (idx < count.y)

#endif

	{
		float x = srcx[idx];
		float y = srcy[idx];
		sy += (float3)(y, y*y, x*y);
	}
	int lid = get_local_id(0);
	barrier(CLK_GLOBAL_MEM_FENCE);
	if (lid == 0)
	{
		g[2] = sy.s0;
		g[3] = sy.s1;
		g[4] = sy.s2;
	}
	for (int n = 1; n < get_local_size(0); n++)
	{
		barrier(CLK_GLOBAL_MEM_FENCE);
		if (lid == n)
		{
			g[2] += sy.s0;
			g[3] += sy.s1;
			g[4] += sy.s2;
		}
	}
	barrier(CLK_GLOBAL_MEM_FENCE);
}

__kernel void kf_CorrCoef(__global float  *g, __global float  *srcx, __global float  *srcy, int2 count)
{
	float2 sx = 0;
	float3 sy = 0;

#ifdef __DEVICE_IS_CPU__

	int idx = get_global_id(0) * count.x;
	int iend = idx + count.x;
	if (iend > count.y) iend = count.y;
	for (; idx < iend; idx++)

#else

	int idx = get_global_id(0);
	for (int i = 0; (i < count.x) && (idx < count.y); i++)
	{
		float x = srcx[idx];
		float y = srcy[idx];
		barrier(CLK_GLOBAL_MEM_FENCE);
		sx += (float2)(x, x*x);
		sy += (float3)(y, y*y, x*y);
		idx += get_local_size(0);
	}
	if (idx < count.y)

#endif

	{
		float x = srcx[idx];
		float y = srcy[idx];
		sx += (float2)(x, x*x);
		sy += (float3)(y, y*y, x*y);
	}
	int lid = get_local_id(0);
	barrier(CLK_GLOBAL_MEM_FENCE);
	if (lid == 0)
	{
		g[0] = sx.s0;
		g[1] = sx.s1;
		g[2] = sy.s0;
		g[3] = sy.s1;
		g[4] = sy.s2;
	}
	for (int n = 1; n < get_local_size(0); n++)
	{
		barrier(CLK_GLOBAL_MEM_FENCE);
		if (lid == n)
		{
			g[0] += sx.s0;
			g[1] += sx.s1;
			g[2] += sy.s0;
			g[3] += sy.s1;
			g[4] += sy.s2;
		}
	}
	barrier(CLK_GLOBAL_MEM_FENCE);
}
//////////////////////////////////////////////////////////////////////////////////

__kernel void k_Sobel_sq(read_only image2d_t src, write_only image2d_t dst, int sq)
{

	int x = get_global_id(0);
	int y = get_global_id(1);
	float Gx, Gy;
	float v;

	if ((x >= get_image_width(dst)) || (y >= get_image_height(dst)))
		return;

	v = read_imagef(src, iSamplerUCN, (int2) ( x - 1, y - 1 )).x;
	Gx = -v;
	Gy = -v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x, y - 1 )).x;
	Gy -= 2 * v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x + 1, y - 1 )).x;
	Gy -= v;
	Gx += v;

	v = read_imagef(src, iSamplerUCN, (int2) ( x - 1, y )).x;
	Gx -= 2 * v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x + 1, y )).x;
	Gx += 2 * v;

	v = read_imagef(src, iSamplerUCN, (int2) ( x - 1, y + 1 )).x;
	Gx -= v;
	Gy += v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x, y + 1 )).x;
	Gy += 2 * v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x + 1, y + 1 )).x;
	Gy += v;
	Gx += v;

	float GG = (Gx*Gx + Gy*Gy);
	v = (sq == 1) ? GG : sqrt(GG);
	write_imagef(dst, (int2) (x, y), (float4)(v,v,v,1.0f));
}

//////////////////////////////////////////////////////////////////////////////////
__kernel void k_Sobel3D_sq(read_only image3d_t src, __global float *dst, int sq)
{
	int x = get_global_id(0);
	int y = get_global_id(1);
	int z = get_global_id(2);
	float Gx, Gy, Gz, v;

	if ((x >= get_image_width(src)) ||
		(y >= get_image_height(src)) ||
		(z >= get_image_depth(src)))
		return;

	Gz = Gx = read_imagef(src, iSamplerUCN, (int4)(x - 1, y, z - 1, 0)).x; //1
	v = read_imagef(src, iSamplerUCN, (int4)(x + 1, y, z - 1, 0)).x;//2
	Gx -= v;
	Gz += v;
	Gy = read_imagef(src, iSamplerUCN, (int4)(x - 1, y - 1, z, 0)).x;//3
	Gx += Gy;
	Gx += 2 * read_imagef(src, iSamplerUCN, (int4)(x - 1, y, z, 0)).x;
	v = read_imagef(src, iSamplerUCN, (int4)(x - 1, y + 1, z, 0)).x;//4
	Gx += v;
	Gy -= v;
	v = read_imagef(src, iSamplerUCN, (int4)(x + 1, y - 1, z, 0)).x;//5
	Gx -= v;
	Gy += v;
	Gx -= 2 * read_imagef(src, iSamplerUCN, (int4)(x + 1, y, z, 0)).x;
	v = read_imagef(src, iSamplerUCN, (int4)(x + 1, y + 1, z, 0)).x;//6
	Gx -= v;
	Gy -= v;
	v = read_imagef(src, iSamplerUCN, (int4)(x - 1, y, z + 1, 0)).x;//7
	Gx += v;
	Gz -= v;
	v = read_imagef(src, iSamplerUCN, (int4)(x + 1, y, z + 1, 0)).x;//8
	Gx -= v;
	Gz -= v;
	/////////////////////////////////////////////////////////////////
	v = read_imagef(src, iSamplerUCN, (int4)(x, y - 1, z - 1, 0)).x;//9
	Gy += v;
	Gz += v;
	v = read_imagef(src, iSamplerUCN, (int4)(x, y + 1, z - 1, 0)).x;//10
	Gy -= v;
	Gz += v;
	Gy += 2 * read_imagef(src, iSamplerUCN, (int4)(x, y - 1, z, 0)).x;
	Gy -= 2 * read_imagef(src, iSamplerUCN, (int4)(x, y + 1, z, 0)).x;
	v = read_imagef(src, iSamplerUCN, (int4)(x, y - 1, z + 1, 0)).x;//11
	Gy += v;
	Gz -= v;
	v = read_imagef(src, iSamplerUCN, (int4)(x, y + 1, z + 1, 0)).x;//12
	Gy -= v;
	Gz -= v;
	//////////////////////////////////////////////////////////////////
	Gz += 2 * read_imagef(src, iSamplerUCN, (int4)(x, y, z - 1, 0)).x;
	Gz -= 2 * read_imagef(src, iSamplerUCN, (int4)(x, y, z + 1, 0)).x;

	float GG = (Gx*Gx + Gy*Gy + Gz*Gz);
	v = (sq == 1) ? GG : sqrt(GG);
	dst[(z * get_image_height(src) + y) * get_image_width(src) + x] = v;
}
//////////////////////////////////////////////////////////////////////////////////

__kernel void k_Scharr_sq(read_only image2d_t src, write_only image2d_t dst, int sq)
{

	int x = get_global_id(0);
	int y = get_global_id(1);
	float Gx, Gy;
	float v;

	if ((x >= get_image_width(dst)) || (y >= get_image_height(dst)))
		return;

	v = read_imagef(src, iSamplerUCN, (int2) ( x - 1, y - 1 )).x*3.0f;
	Gx = -v;
	Gy = -v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x, y - 1 )).x;
	Gy -= 10.0f*v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x + 1, y - 1 )).x*3.0f;
	Gy -= v;
	Gx += v;

	v = read_imagef(src, iSamplerUCN, (int2) ( x - 1, y )).x;
	Gx -= 10.0f * v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x + 1, y )).x;
	Gx += 10.0f * v;

	v = read_imagef(src, iSamplerUCN, (int2) ( x - 1, y + 1 )).x*3.0f;
	Gx -= v;
	Gy += v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x, y + 1 )).x;
	Gy += 10.0f * v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x + 1, y + 1 )).x*3.0f;
	Gy += v;
	Gx += v;

	float GG = (Gx*Gx + Gy*Gy);
	v = (sq == 1) ? GG : sqrt(GG);
	write_imagef(dst, (int2) (x, y), (float4)(v,v,v,1.0f));
}

//////////////////////////////////////////////////////////////////////////////////

__kernel void k_Prewitt_sq(read_only image2d_t src, write_only image2d_t dst, int sq)
{

	int x = get_global_id(0);
	int y = get_global_id(1);
	float Gx, Gy;
	float v;

	if ((x >= get_image_width(dst)) || (y >= get_image_height(dst)))
		return;

	v = read_imagef(src, iSamplerUCN, (int2) ( x - 1, y - 1 )).x;
	Gx = -v;
	Gy = -v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x, y - 1 )).x;
	Gy -= v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x + 1, y - 1 )).x;
	Gy -= v;
	Gx += v;

	v = read_imagef(src, iSamplerUCN, (int2) ( x - 1, y )).x;
	Gx -= v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x + 1, y )).x;
	Gx += v;

	v = read_imagef(src, iSamplerUCN, (int2) ( x - 1, y + 1 )).x;
	Gx -= v;
	Gy += v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x, y + 1 )).x;
	Gy += v;
	v = read_imagef(src, iSamplerUCN, (int2) ( x + 1, y + 1 )).x;
	Gy += v;
	Gx += v;

	float GG = (Gx*Gx + Gy*Gy);
	v = (sq == 1) ? GG : sqrt(GG);
	write_imagef(dst, (int2) (x, y), (float4)(v,v,v,1.0f));
}

//////////////////////////////////////////////////////////////////////////////////

__kernel void k_EmbossLaplacian_sq(read_only image2d_t src, write_only image2d_t dst, int sq)
{

	int x = get_global_id(0);
	int y = get_global_id(1);
	float G;
	float v;

	if ((x >= get_image_width(dst)) || (y >= get_image_height(dst)))
		return;

	v = read_imagef(src, iSamplerUCN, (int2)(x - 1, y - 1)).x;
	G = -v;
	v = read_imagef(src, iSamplerUCN, (int2)(x + 1, y - 1)).x;
	G -= v;

	v = read_imagef(src, iSamplerUCN, (int2)(x, y)).x;
	G += 4.0f * v;

	v = read_imagef(src, iSamplerUCN, (int2)(x - 1, y + 1)).x;
	G -= v;
	v = read_imagef(src, iSamplerUCN, (int2)(x + 1, y + 1)).x;
	G -= v;

	v = (sq == 1) ? G*G : G;
	write_imagef(dst, (int2)(x, y), (float4)(v,v,v,1.0f));
}

//////////////////////////////////////////////////////////////////////////////////

__kernel
void k_separableH(read_only image2d_t src, write_only image2d_t dst, __global float2 *mask, int maskSize)
{
	int r = get_global_id(0);
	float sum;
	float2 bmask[16];
	if (r < get_image_height(src))
	{
		float scalex = 1.0f / (float)(get_image_width(src) - 1);
		float scaley = 1.0f / (float)(get_image_height(src) - 1);
		float ypos = (float)r * scaley;
		for (int j = 0; j < maskSize; j++)
		{
			bmask[j] = mask[j]*(float2)(1.0f, scalex);
		}

		for (int c = 0; c < get_image_width(src); c++)
		{
			float xpos = (float)c * scalex;
			sum = read_imagef(src, iSamplerNCL, (float2)(xpos, ypos)).x * bmask[0].x;
			for (int j = 1; j < maskSize; j++)
			{
				sum += read_imagef(src, iSamplerNCL, (float2)(xpos - bmask[j].y, ypos)).x * bmask[j].x;
				sum += read_imagef(src, iSamplerNCL, (float2)(xpos + bmask[j].y, ypos)).x * bmask[j].x;
			}
			write_imagef(dst, (int2)(r, c), (float4)(sum,sum,sum,1.0f));
		}
	}
}

__kernel
void k_separableV(read_only image2d_t src, write_only image2d_t dst, __global float2 *mask, int maskSize)
{
	int c = get_global_id(0);
	float sum;
	float2 bmask[16];
	if (c < get_image_width(src))
	{
		float scalex = 1.0f / (float)(get_image_width(src) - 1);
		float scaley = 1.0f / (float)(get_image_height(src) - 1);
		float xpos = (float)c * scalex;
		for (int j = 0; j < maskSize; j++)
		{
			bmask[j] = mask[j] * (float2)(1.0f, scaley);
		}

		for (int r = 0; r < get_image_height(src); r++)
		{
			float ypos = (float)r * scaley;
			sum = read_imagef(src, iSamplerNCL, (float2)(xpos, ypos)).x * bmask[0].x;
			for (int j = 1; j < maskSize; j++)
			{
				sum += read_imagef(src, iSamplerNCL, (float2)(xpos, ypos - bmask[j].y)).x * bmask[j].x;
				sum += read_imagef(src, iSamplerNCL, (float2)(xpos, ypos + bmask[j].y)).x * bmask[j].x;
			}
			write_imagef(dst, (int2)(r, c), (float4)(sum,sum,sum,1.0f));
		}
	}
}

__attribute__((always_inline))
void fft_2048_0_pass(
	__local float2 *sMem,
	float2 *a,
	int lid,
	float dir)
{
	float2 w, w0;
	float ang;
	__local float2 *lMemStore, *lMemLoad;

	fftKernel8(a + 0, dir);
	ang = dir * ( M_PI_F / 1024.0f ) * lid;
	w = w0 = (float2)(native_cos(ang), native_sin(ang));
	lMemStore = sMem + lid;
	lMemStore[0] = a[0];
	lMemStore[258] = complexMul(a[1], w);		w = complexMul(w0, w);
	lMemStore[516] = complexMul(a[2], w);		w = complexMul(w0, w);
	lMemStore[774] = complexMul(a[3], w);		w = complexMul(w0, w);
	lMemStore[1032] = complexMul(a[4], w);		w = complexMul(w0, w);
	lMemStore[1290] = complexMul(a[5], w);		w = complexMul(w0, w);
	lMemStore[1548] = complexMul(a[6], w);		w = complexMul(w0, w);
	lMemStore[1806] = complexMul(a[7], w);
	barrier(CLK_LOCAL_MEM_FENCE);
	lMemLoad = sMem + mad24(lid & 7, 258, lid >> 3);
	a[0] = lMemLoad[0];
	a[1] = lMemLoad[32];
	a[2] = lMemLoad[64];
	a[3] = lMemLoad[96];
	a[4] = lMemLoad[128];
	a[5] = lMemLoad[160];
	a[6] = lMemLoad[192];
	a[7] = lMemLoad[224];
	barrier(CLK_LOCAL_MEM_FENCE);
	fftKernel8(a + 0, dir);
	ang = dir * ( M_PI_F / 128.0f ) * (lid >> 3);
	w = w0 = (float2)(native_cos(ang), native_sin(ang));
	lMemStore[0] = a[0];
	lMemStore[264] = complexMul(a[1], w);		w = complexMul(w0, w);
	lMemStore[528] = complexMul(a[2], w);		w = complexMul(w0, w);
	lMemStore[792] = complexMul(a[3], w);		w = complexMul(w0, w);
	lMemStore[1056] = complexMul(a[4], w);		w = complexMul(w0, w);
	lMemStore[1320] = complexMul(a[5], w);		w = complexMul(w0, w);
	lMemStore[1584] = complexMul(a[6], w);		w = complexMul(w0, w);
	lMemStore[1848] = complexMul(a[7], w);
	barrier(CLK_LOCAL_MEM_FENCE);
	lMemLoad = sMem + mad24((lid & 63) >> 3, 264, mad24(lid >> 6, 8, lid & 7));
	a[0] = lMemLoad[0];
	a[1] = lMemLoad[32];
	a[2] = lMemLoad[64];
	a[3] = lMemLoad[96];
	a[4] = lMemLoad[128];
	a[5] = lMemLoad[160];
	a[6] = lMemLoad[192];
	a[7] = lMemLoad[224];
	barrier(CLK_LOCAL_MEM_FENCE);
	fftKernel8(a + 0, dir);
	ang = dir * ( M_PI_F / 16.0f ) * (lid >> 6);
	w = w0 = (float2)(native_cos(ang), native_sin(ang));
	lMemStore[0] = a[0];
	lMemStore[256] = complexMul(a[1], w);		w = complexMul(w0, w);
	lMemStore[512] = complexMul(a[2], w);		w = complexMul(w0, w);
	lMemStore[768] = complexMul(a[3], w);		w = complexMul(w0, w);
	lMemStore[1024] = complexMul(a[4], w);		w = complexMul(w0, w);
	lMemStore[1280] = complexMul(a[5], w);		w = complexMul(w0, w);
	lMemStore[1536] = complexMul(a[6], w);		w = complexMul(w0, w);
	lMemStore[1792] = complexMul(a[7], w);
	barrier(CLK_LOCAL_MEM_FENCE);
	lMemLoad = sMem + mad24(lid >> 6, 256, lid & 63);
	a[0] = lMemLoad[0];
	a[1] = lMemLoad[64];
	a[2] = lMemLoad[128];
	a[3] = lMemLoad[192];
	a[4] = lMemLoad[1024];
	a[5] = lMemLoad[1088];
	a[6] = lMemLoad[1152];
	a[7] = lMemLoad[1216];
	barrier(CLK_LOCAL_MEM_FENCE);
	fftKernel4(a, 0, 1, 2, 3, dir);
	fftKernel4(a, 4, 5, 6, 7, dir);
}


/* Run kernel with global dim = {256*BatchSize}, local dim={256} */
__kernel __attribute__((reqd_work_group_size (256,1,1)))
void k_fft_conv_2048(
	__global float *in, 
	__global float2 *out,
	__global float *filter, 
	int4 size)
{
	__local float2 sMem[2112];
	float2 a[8],t;
	int lid = get_local_id( 0 );
	int groupId = get_group_id( 0 );
	unsigned offset =  groupId * 2 * size.x + lid;

	__global float *re = in + offset;
	__global float *im = re + size.x;
	int chk = size.x - lid;
	a[0] = (chk <= 0) ? (float2)0 : (float2)(re[0], im[0]);
	a[1] = (chk <= 256) ? (float2)0 : (float2)(re[256], im[256]);
	a[2] = (chk <= 512) ? (float2)0 : (float2)(re[512], im[512]);
	a[3] = (chk <= 768) ? (float2)0 : (float2)(re[768], im[768]);
	a[4] = a[5] = a[6] = a[7] = (float2)0;
	
	fft_2048_0_pass(sMem, a, lid, 1.0f);
	
	filter += lid;
	a[0] *= filter[0];
	t = a[1];
	a[1] = a[4] * filter[256];
	a[4] = a[2] * filter[1024];
	a[2] = t * filter[512];
	t = a[3];
	a[3] = a[5] * filter[768];
	a[5] = a[6] * filter[1280];
	a[6] = t * filter[1536];
	a[7] *= filter[1792];

	fft_2048_0_pass(sMem, a, lid, -1.0f);
	int tmp = size.x * size.y;
	int z = offset / tmp;
	int y = (offset - z * tmp) / size.x;
	int x = offset - z * tmp - y * size.x;
	if(chk > 0)
	{
		offset = ((x * size.z + z) * size.y + y)/2;
		out[offset] = a[0] * 4.8828125000e-004f;
	}
	if(chk > 256)
	{
		tmp = 128 * size.z * size.y;
		offset += tmp;
		out[offset] = a[4] * 4.8828125000e-004f; 
	}
	if(chk > 512)
	{
		offset += tmp;
		out[offset] = a[1] * 4.8828125000e-004f;
	}
	if(chk > 768)
	{
		offset += tmp;
		out[offset] = a[5] * 4.8828125000e-004f;
	}
}




//
//	Multiply m1 with m2
//  m1[r x q] m2[q x p] mres[r x p]
//
//  CL threads dimension x:r, y:q
//
__kernel void doMatrixMultF(__global float * mRes, __global float * m1, __global float * m2,
		int r, int q, int p)
{
    int i = get_global_id(0);
	if (i >= r)
		return;
    
	int j = get_global_id(1);
	if (j >= p)
		return;

    int writeIdx = j + p * i;
	float res = 0;

    int strideM1 = i * q;
    for (int k = 0; k < q; k++)
    {
        res += (m1[strideM1 + k] * m2[j + p * k]);
    }
    mRes[writeIdx] = res;
}
//
//	Multiply m1 with m2
//  m1[r x q] m2[q x p] mres[r x p]
//
//  CL threads dimension x:p
//
__kernel void doMatrixMultLocalF(__global float * mRes, __global float * m1, __global float * m2,
		int r, int q, int p)
{
    int j = get_global_id(0);
	if (j >= p)
		return;
    for(int i = 0; i < r; i++)
	{ 
		int writeIdx = j + p * i;
		float res = 0;
	
		int strideM1 = i * q;
		for (int k = 0; k < q; k++)
		{
			res += (m1[strideM1 + k] * m2[j + p * k]);
		}
		mRes[writeIdx] = res;
	}
}

//
// Given a limit of max zeros counted per column of 
// a mat[rows x columns] set each column over the 
// limit to zero
//
//  CL threads dimension x:columns
//
__kernel void doMatrixColDiscardNegative(__global float * mat,	
	int rows, int columns, int limit)
{
    int i = get_global_id(0);
	if (i >= columns)
		return;

	int count = 0;
    for (int k = 0; k < rows; k++)
	{
		if(mat[i + k * columns] < 0)
		    count++;
	}

	if (count >= limit)
	    for (int k = 0; k < rows; k++)
		{
			mat[i + k * columns] = 0;
		}
}


__kernel void doFramesAveraging(__global float *src, __global float *dst,
	int frames, int frameLength)
{
	int idx = get_global_id(0);
	src += idx;
	if (idx < frameLength)
	{
		float acc = 0;
		for(int i =0; i < frames; i++)
		{
			acc += *src;
			src += frameLength;
		}
		dst[idx] = acc/frames;
	}
}

__kernel void doFramesAveragingLocal(__global float *src, __global float *dst,
	int frames, int frameLength)
{
	__local float lb[8192];
	int idx = get_global_id(0);
	int lid = get_local_id(0);
	src += idx;
	float acc = 0;
	for(int i =0; i < frames; i++)
	{
		if(lid == 0)
		{
			event_t ev = async_work_group_copy (lb,src, get_local_size(0), 0); 
			wait_group_events(1,&ev);
		}
		barrier(CLK_GLOBAL_MEM_FENCE);
		acc += lb[lid];
		src += frameLength;
	}
	dst[idx] = acc/frames;
}


