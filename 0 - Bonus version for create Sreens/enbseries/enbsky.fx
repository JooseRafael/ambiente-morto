
/* VERSION 2.0 */

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//      /XXXXXXX\          __                  |HHHHHHHHH\                                          \HH\     /HH/    // 
//     |XX|   |XX|         XX                  |HH|    \HHH\   /\                                    \HH\   /HH/     //
//     |XX|    |X|        /XX\                 |HH|      \HH\  \/                                     \HH\ /HH/      //
//      \XX\             /XXXX\                |HH|       |HH| __ __  __     ___     ___      __       \HHHHH/       //
//       \XX\           /XX/\XX\               |HH|       |HH| || || /__\   /___\   /___\   __||__      \HHH/        //
//     _   \XX\        /XX/  \XX\              |HH|       |HH| || ||//  \\ //   \\ //   \\ |==  ==|     /HHH\        //
//    |X|    \XX\     /XXXXXXXXXX\             |HH|       |HH| || ||/      ||===|| ||    -    ||       /HHHHH\       //
//    \XX\    |XX|   /XX/      \XX\            |HH|      /HH/  || ||       ||      ||         ||      /HH/ \HH\      // 
//     \XX\   /XX/  /XX/        \XX\  ________ |HH|    /HHH/   || ||       \\    _ \\    _    ||  _  /HH/   \HH\     //
//      \XXXXXXX/  /XX/          \XX\|________||HHHHHHHHH/     || ||        \===//  \===//    \\=// /HH/     \HH\    //
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++            ENBSeries effect file              ++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++      SA_DirectX by Maxim Dubinov(Makarus)     ++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++    Visit http://www.facebook.com/sadirectx    ++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++    https://www.youtube.com/channel/UCrASy-x5DgwHpYiDv41RL2Q    ++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++          Visit http://enbdev.com              ++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++    Copyright (c) 2007-2018 Boris Vorontsov    ++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

float4 tempF1; float4 tempF2; float4 tempF3;
float4 ScreenSize; float ENightDayFactor; float EInteriorFactor; float4 WeatherAndTime; float4 Timer;
float FieldOfView; float GameTime; float4 SunDirection; float4 CustomShaderConstants1[8]; float4 MatrixVP[4];
float4 MatrixInverseVP[4]; float4 MatrixVPRotation[4]; float4 MatrixInverseVPRotation[4]; float4 MatrixView[4];
float4 MatrixInverseView[4]; float4 CameraPosition; float4x4 MatrixWVP; float4x4 MatrixWVPInverse; float4x4 MatrixWorld;
float4x4 MatrixProj; float4 FogParam; float4 FogFarColor;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//Textures
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

texture2D texDepth;
texture2D texNoise;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//Sampler Inputs
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

sampler2D SamplerNoise = sampler_state
{
	Texture   = <texNoise>;
	MinFilter = POINT;
	MagFilter = POINT;
	MipFilter = NONE;
	AddressU  = Wrap;
	AddressV  = Wrap;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerDepth = sampler_state
{
	Texture   = <texDepth>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Wrap;
	AddressV = Wrap;
	AddressW = Wrap;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

struct VS_INPUT
{
	float3	pos : POSITION;
	float4	diff : COLOR0;
};

struct VS_OUTPUT
{
	float4	pos : POSITION;
	float4	diff : COLOR0;
	float4	vposition : TEXCOORD6;
};

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

VS_OUTPUT	VS_Draw(VS_INPUT IN)
{
    VS_OUTPUT OUT;

	float4	pos=float4(IN.pos.x,IN.pos.y,IN.pos.z,1.0);

	float4	tpos=mul(pos, MatrixWVP);

	OUT.diff=IN.diff;
	OUT.pos=tpos;
	OUT.vposition=tpos;

    return OUT;
}

float4 	WorldPos(in float2 coord)
{
 	const float	depth=1.0;  
   
   float4 tvec; 
          tvec.xy = coord.xy*2.0-1.0;   
          tvec.y = -tvec.y;   
          tvec.z = depth;
          tvec.w = 1.0;
   float4 wpos;
          wpos.x = dot(tvec, MatrixInverseVPRotation[0]);
          wpos.y = dot(tvec, MatrixInverseVPRotation[1]);
          wpos.z = dot(tvec, MatrixInverseVPRotation[2]);
          wpos.w = dot(tvec, MatrixInverseVPRotation[3]);
          wpos.xyz/= wpos.w;
   return wpos;
	
   return wpos;
}

float3 CalculateSun(in float2 coord)
{	  
    float3 sv = SunDirection.xyz;
    float3 sv2 = normalize(float3(-0.0833, -0.946, 0.317));
	
   float t = GameTime;	
	
   float x1 = smoothstep(0.0, 4.0, t);
   float x2 = smoothstep(4.0, 5.0, t);
   float x3 = smoothstep(5.0, 6.0, t);
   float x4 = smoothstep(6.0, 7.0, t);
   float xE = smoothstep(8.0, 11.0, t);
   float x5 = smoothstep(16.0, 17.0, t);
   float x6 = smoothstep(18.0, 19.0, t);
   float x7 = smoothstep(19.0, 20.0, t);
   float xG = smoothstep(20.0, 21.0, t);  
   float xZ = smoothstep(21.0, 22.0, t);
   float x8 = smoothstep(22.0, 23.0, t);
   float x9 = smoothstep(23.0, 24.0, t);
   
   float3 t0 = lerp(0.0, 0.0, x1);
          t0 = lerp(t0, float3(1.0, 0.549, 0.0784)*2.0, x2);
          t0 = lerp(t0, float3(1.0, 0.549, 0.0784)*2.7, x3);
          t0 = lerp(t0, float3(1.0, 0.51, 0.235)*2.0, x4);
          t0 = lerp(t0, float3(1.0, 0.863, 0.549), xE);
          t0 = lerp(t0, float3(1.0, 0.863, 0.549), x5);
          t0 = lerp(t0, float3(1.0, 0.549, 0.0784)*2.0, x6);	 
          t0 = lerp(t0, float3(1.0, 0.549, 0.0784)*0.8, x7);
		  t0 = lerp(t0, float3(1.0, 0.549, 0.0784)*0.6, xG);
		  t0 = lerp(t0, float3(1.0, 0.549, 0.0784)*0.4, xZ);
          t0 = lerp(t0, float3(1.0, 0.549, 0.0784)*0.2, x8);
          t0 = lerp(t0, 0.0, x9);	

   float3 t3 = lerp(0.0, 0.0, x1);
          t3 = lerp(t3, float3(1.0, 0.447, 0.0)*2.0, x2);
          t3 = lerp(t3, float3(1.0, 0.447, 0.0)*2.7, x3);
          t3 = lerp(t3, float3(1.0, 0.784, 0.392)*2.0, x4);
          t3 = lerp(t3, float3(1.0, 1.0, 1.0)*2.0, xE);
          t3 = lerp(t3, float3(1.0, 1.0, 1.0)*2.0, x5);
          t3 = lerp(t3, float3(1.0, 0.447, 0.0)*2.0, x6);		 
          t3 = lerp(t3, float3(1.0, 0.447, 0.0)*0.8, x7);
		  t3 = lerp(t3, float3(1.0, 0.447, 0.0)*0.6, xG);
		  t3 = lerp(t3, float3(1.0, 0.447, 0.0)*0.4, xZ);
          t3 = lerp(t3, float3(1.0, 0.447, 0.0)*0.2, x8);
          t3 = lerp(t3, 0.0, x9); 	
   
   float4 wp = WorldPos(coord);
   float c0 = 475.0 * 25.0;
   float c1 = 18.0;
   float c2 = 1.35;
   
   float3 np0 = normalize(wp.xyz);
   float factor = (0.01/12.0) - dot(-sv, np0);
         factor = pow(factor, c0);
   float factor1 = 0.04 - dot(-sv, np0);
         factor1 = pow(factor1, c1);	
   float factor2 = 0.6 - dot(-sv, np0);
         factor2 = pow(factor2, c2);
		 
   float factor3 = 0.04 - dot(-sv2, np0);
         factor3 = pow(factor3, c1);	
   float factor4 = 0.6 - dot(-sv2, np0);
         factor4 = pow(factor4, c2);		 
		 
   float3 f0 = factor/12.0;
   float3 f1 = (factor1*factor1)/5.0;    
   float3 f2 = (factor2*factor2)/10.0;
   
   float3 f1x = (factor3*factor3)/5.0;    
   float3 f2x = (factor4*factor4)/10.0;   
   float3 fnight = (f1x*float3(0.0863, 0.137, 0.176))+(f2x*float3(0.0, 0.0, 0.0));     
		  fnight*= smoothstep(0.0, 0.3, -np0.y);     
   
   float y1 = smoothstep(0.0, 2.0, t);
   float yZ = smoothstep(2.0, 3.0, t);   
   float y2 = smoothstep(4.0, 23.0, t);
   float y3 = smoothstep(23.0, 24.0, t);   
   
  float3 t4 = lerp(0.0, 0.0, y1);
         t4 = lerp(t4, float3(0.0392, 0.0235, 0.0118), yZ);  
         t4 = lerp(t4, float3(0.0392, 0.0235, 0.0118), y2);
         t4 = lerp(t4, 0.0, y3);
		 
  float3 t5 = lerp(1.0, 0.5, y1);
         t5 = lerp(t5, 0.0, yZ);  
         t5 = lerp(t5, 0.0, y2);
         t5 = lerp(t5, 1.0, y3);	 
 
   float3 wp0 = wp + CameraPosition;  
   float3 np1 = normalize(wp0.xyz);
   float3 f3 = (f1*t0)+(f2*t3)+(fnight*t5);
 
float3 SunCurrent;
float3 SunNext;
float4 WeatherX = WeatherAndTime;

if (WeatherX.x==0,1) SunCurrent = f3;
if (WeatherX.y==0,1) SunNext = f3;

if (WeatherX.x==4) SunCurrent = 0.0;
if (WeatherX.x==7) SunCurrent = 0.0;
if (WeatherX.x==8) SunCurrent = 0.0;
if (WeatherX.x==9) SunCurrent = 0.0;
if (WeatherX.x==12) SunCurrent = 0.0;
if (WeatherX.x==15) SunCurrent = 0.0;
if (WeatherX.x==16) SunCurrent = 0.0;

if (WeatherX.y==4) SunNext = 0.0;
if (WeatherX.y==7) SunNext = 0.0;
if (WeatherX.y==8) SunNext = 0.0;
if (WeatherX.y==9) SunNext = 0.0;
if (WeatherX.y==12) SunNext = 0.0;
if (WeatherX.y==15) SunNext = 0.0;
if (WeatherX.y==16) SunNext = 0.0;	  

float3 sss2 = lerp(SunCurrent, SunNext, WeatherX.z);	   
 
   return sss2;
}

float CalculateGameTime0(in float t)
{	  
   float x1 = smoothstep(0.0, 4.0, t);
   float x2 = smoothstep(4.0, 5.0, t);
   float x3 = smoothstep(5.0, 6.0, t);
   float x4 = smoothstep(6.0, 7.0, t);
   float xE = smoothstep(8.0, 11.0, t);
   float x5 = smoothstep(16.0, 17.0, t);
   float x6 = smoothstep(18.0, 19.0, t);
   float x7 = smoothstep(19.0, 20.0, t);
   float xG = smoothstep(20.0, 21.0, t);  
   float xZ = smoothstep(21.0, 22.0, t);
   float x8 = smoothstep(22.0, 23.0, t);
   float x9 = smoothstep(23.0, 24.0, t);
   
   float3 t0 = lerp(0.0, 0.1, x1);
          t0 = lerp(t0, 0.2, x2);
          t0 = lerp(t0, 0.8, x3);
          t0 = lerp(t0, 0.9, x4);
          t0 = lerp(t0, 1.0, xE);
          t0 = lerp(t0, 1.0, x5);
          t0 = lerp(t0, 0.9, x6);	 
          t0 = lerp(t0, 0.5, x7);
		  t0 = lerp(t0, 0.4, xG);
		  t0 = lerp(t0, 0.3, xZ);
          t0 = lerp(t0, 0.2, x8);
          t0 = lerp(t0, 0.0, x9); 		  
   return t0;	  
}

float CalculateGameTime(in float t)
{	  
   float x1 = smoothstep(0.0, 4.0, t);
   float x2 = smoothstep(4.0, 5.0, t);
   float x3 = smoothstep(5.0, 6.0, t);
   float x4 = smoothstep(6.0, 7.0, t);
   float xE = smoothstep(8.0, 11.0, t);
   float x5 = smoothstep(16.0, 17.0, t);
   float x6 = smoothstep(18.0, 19.0, t);
   float x7 = smoothstep(19.0, 20.0, t);
   float xG = smoothstep(20.0, 21.0, t);  
   float xZ = smoothstep(21.0, 22.0, t);
   float x8 = smoothstep(22.0, 23.0, t);
   float x9 = smoothstep(23.0, 24.0, t);
   
   float3 t0 = lerp(0.0, 0.1, x1);
          t0 = lerp(t0, 0.7, x2);
          t0 = lerp(t0, 1.0, x3);
          t0 = lerp(t0, 1.0, x4);
          t0 = lerp(t0, 1.0, xE);
          t0 = lerp(t0, 1.0, x5);
          t0 = lerp(t0, 0.9, x6);		 
          t0 = lerp(t0, 0.8, x7);
		  t0 = lerp(t0, 0.6, xG);
		  t0 = lerp(t0, 0.4, xZ); 
          t0 = lerp(t0, 0.2, x8);
          t0 = lerp(t0, 0.0, x9); 		  
   return t0;	  
}

float CalculateGameTime2(in float t)
{	  
   float x1 = smoothstep(0.0, 4.0, t);
   float x2 = smoothstep(4.0, 5.0, t);
   float x3 = smoothstep(5.0, 6.0, t);
   float x4 = smoothstep(6.0, 7.0, t);
   float xE = smoothstep(8.0, 11.0, t);
   float x5 = smoothstep(16.0, 17.0, t);
   float x6 = smoothstep(18.0, 19.0, t);
   float x7 = smoothstep(19.0, 20.0, t);
   float xG = smoothstep(20.0, 21.0, t);  
   float xZ = smoothstep(21.0, 22.0, t);
   float x8 = smoothstep(22.0, 23.0, t);
   float x9 = smoothstep(23.0, 24.0, t);
   
   float3 t0 = lerp(0.0, 0.0, x1);
          t0 = lerp(t0, 0.7, x2);
          t0 = lerp(t0, 0.9, x3);
          t0 = lerp(t0, 1.0, x4);
          t0 = lerp(t0, 1.0, xE);
          t0 = lerp(t0, 0.5, x5);
          t0 = lerp(t0, 0.5, x6);		 
          t0 = lerp(t0, 1.0, x7);
		  t0 = lerp(t0, 1.0, xG);
		  t0 = lerp(t0, 1.0, xZ);
          t0 = lerp(t0, 0.0, x8);  
          t0 = lerp(t0, 0.0, x9);	  
   return t0;	  
}

float3 ColorHorizon(float3 r0, float f)
{	
   float t = GameTime;
   float tf = CalculateGameTime2(t); 
   float df = pow(saturate(SunDirection.z * tf), 1.0);
	float3 sg = lerp(1.0, float3(0.902, 0.784, 0.882), pow(r0.y, 2.2));
	float3 sg2 = lerp(1.0, float3(1.0, 1.0, 0.529), pow(r0.y, 1.15));
	float3 sd = lerp(1.0, float3(1.0, 1.0, 0.804), pow(r0.y, 2.15));
	float3 sd2 = lerp(1.0, float3(1.0, 1.0, 1.0), pow(r0.y, 1.15));
	float3 a = lerp(sg, 1.0, f);	
	float3 a2 = lerp(sg2, 1.0, f);
	float3 ad = lerp(sd, 1.0, f);	
	float3 ad2 = lerp(sd2, 1.0, f);
	r0 *= lerp(a, 1.0, df);	
	r0 *= lerp(a2, 1.0, df);
	float3 d1 = lerp(1.0, ad, df);
	float3 d2 = lerp(1.0, ad2, df);
    r0 *= d1*d2;				
   return r0;
}

float3 ColorTop()
{	
float3   cc = float3(1.0, 1.0, 1.0);
float3   cn = float3(1.0, 1.0, 1.0);

float x1 = smoothstep(1.0, 2.0,  WeatherAndTime);
   
float4 WeatherX = WeatherAndTime;

if (WeatherX.x==0,1) cc = float3(0.251, 0.51, 1.0);
if (WeatherX.y==0,1) cn = float3(0.251, 0.51, 1.0);

float3 skyW = float3(0.259, 0.259, 0.259);

if (WeatherX.x==4) cc = skyW;
if (WeatherX.x==7) cc = skyW;
if (WeatherX.x==8) cc = skyW;
if (WeatherX.x==9) cc = skyW;
if (WeatherX.x==12) cc = skyW;
if (WeatherX.x==15) cc = skyW;
if (WeatherX.x==16) cc = skyW;

if (WeatherX.y==4) cn = skyW;
if (WeatherX.y==7) cn = skyW;
if (WeatherX.y==8) cn = skyW;
if (WeatherX.y==9) cn = skyW;
if (WeatherX.y==12) cn = skyW;
if (WeatherX.y==15) cn = skyW;
if (WeatherX.y==16) cn = skyW;

float3 sss = lerp(cc, cn, WeatherX.z);

return sss;
}

float 	CalculateLuminance(in float3 color)
{
	return (color.r * 0.2126 + color.g * 0.7152 + color.b * 0.0722)*1.0;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

float4 PS_Draw(VS_OUTPUT IN, float2 vPos : VPOS) : COLOR
{
	

   float t = GameTime;
   float tf = CalculateGameTime(t);   
   float tf2 = CalculateGameTime0(t);    
   
 	float2	coord=vPos.xy*ScreenSize.y;
	        coord.y*=ScreenSize.z;  
   
   float4 r0 = float4(1.0, 1.0, 1.0, 1.0);
   float4 wpos = WorldPos(coord);
   float3 colorS = ColorTop();  
//Sky------------------------------------------------------------------------------	
   float4 w0 = wpos+CameraPosition;
   float4 cam = CameraPosition; 
   float p2 = length(wpos.xyz);  
	 	  
   float4 np0 = float4(normalize(wpos.xyz), 1.0);
	      np0.xyz = normalize(w0.xyz-float3(cam.xy, 0.0));	
		  
   float3 vec = normalize(float3(0.0, 0.0, 1.0));	  
   float3 hv = normalize(-vec+np0);
   float sgf = dot(hv, np0);    
   
 	if (dot(hv, np0) > 0.75*0.95)
		  sgf = (1.5*0.95) - sgf;
	float curve = 8.0*0.85;	
	      sgf = pow(sgf, curve);  

	float3 d1 = float3(1.0, 1.0, 1.0);	
 		   d1 = CalculateLuminance(d1) * colorS;
		   d1 *= (sgf * (20.0*1.4)) + 0.85;
			   
   float3 sun = CalculateSun(coord);		   
   float3 d2 = ColorHorizon(d1, sgf);
   
	d1 = d2;
	d1*= lerp(0.0, 1.0, tf);
	d1.xyz+=0.000001;
	float3 n0 = normalize(d1.xyz);
	float3 ct0=d1.xyz/n0.xyz;
	ct0=pow(ct0, 0.70);
    n0.xyz = pow(n0.xyz, 0.62);   
    d1.xyz = ct0*n0.xyz;  
	
    float np4 = (np0.z * (2.0)) + 0.40;
	      np4 *= lerp(0.0, 1.0, tf2);  		  
    float3 a = d1+sun;
    float3 a2 = np4;
    float3   sc = 1.0;
    float3   sn = 1.0;   
    float4 wx = WeatherAndTime;	
    if (wx.x==0,1) sc = a*0.6;
    if (wx.y==0,1) sn = a*0.6;   
    if (wx.x==4) sc = a2;
    if (wx.x==7) sc = a2;
    if (wx.x==8) sc = a2;
    if (wx.x==9) sc = a2;
    if (wx.x==12) sc = a2;
    if (wx.x==15) sc = a2;
    if (wx.x==16) sc = a2;
    if (wx.y==4) sn = a2;
    if (wx.y==7) sn = a2;
    if (wx.y==8) sn = a2;
    if (wx.y==9) sn = a2;
    if (wx.y==12) sn = a2;
    if (wx.y==15) sn = a2;
    if (wx.y==16) sn = a2;
    float3 sss = lerp(sc, sn, wx.z);			
    r0.xyz = sss;	

	r0.w=1.0;  
	return r0;

}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

technique Draw
{
    pass p0
    {
	VertexShader = compile vs_3_0 VS_Draw();
	PixelShader  = compile ps_3_0 PS_Draw();
	}
}


