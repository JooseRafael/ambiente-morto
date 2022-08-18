
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
//++++++++++++++++++++++++++++++++++    https://www.youtube.com/user/XMakarusX     ++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++          Visit http://enbdev.com              ++++++++++++++++++++++++++++++++++//
//++++++++++++++++++++++++++++++++++    Copyright (c) 2007-2018 Boris Vorontsov    ++++++++++++++++++++++++++++++++++//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

float4	tempF1; float4	tempF2; float4	tempF3;
float4 ScreenSize; float4 Timer;float ENightDayFactor;
float4 SunDirection; float EInteriorFactor; float FadeFactor; float FieldOfView;
float4	MatrixVP[4]; float4 MatrixInverseVP[4]; float4 MatrixVPRotation[4]; float4 MatrixInverseVPRotation[4];
float4 MatrixView[4];float4 MatrixInverseView[4];float4 CameraPosition;float GameTime;float4 CustomShaderConstants1[8];
float4	WeatherAndTime;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//textures
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
texture2D texColor;
texture2D texDepth;

sampler2D SamplerColor = sampler_state
{
	Texture   = <texColor>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerDepth = sampler_state
{
	Texture   = <texDepth>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = NONE;
	AddressU  = Wrap;
	AddressV  = Wrap;
	SRGBTexture=TRUE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};


struct VS_OUTPUT_POST
{
	float4 vpos  : POSITION;
	float2 txcoord : TEXCOORD0;
};

struct VS_INPUT_POST
{
	float3 pos  : POSITION;
	float2 txcoord : TEXCOORD0;
};

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

VS_OUTPUT_POST VS_PostProcess(VS_INPUT_POST IN)
{
	VS_OUTPUT_POST OUT;

	float4 pos=float4(IN.pos.x,IN.pos.y,IN.pos.z,1.0);

	OUT.vpos=pos;
	OUT.txcoord.xy=IN.txcoord.xy;

	return OUT;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

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
   
   float3 t0 = lerp(0.04, 0.1, x1);
          t0 = lerp(t0, 0.1, x2);
          t0 = lerp(t0, 0.8, x3); 
          t0 = lerp(t0, 0.9, x4);
          t0 = lerp(t0, 1.0, xE);
          t0 = lerp(t0, 1.0, x5);
          t0 = lerp(t0, 0.9, x6);		 
          t0 = lerp(t0, 0.5, x7);
		  t0 = lerp(t0, 0.4, xG);
		  t0 = lerp(t0, 0.3, xZ);
          t0 = lerp(t0, 0.2, x8); 
          t0 = lerp(t0, 0.04, x9);  		  
   return t0;	  
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

float4 PS_DX0(VS_OUTPUT_POST IN, float2 vPos : VPOS) : COLOR
{
   float4 r0=tex2D(SamplerColor, IN.txcoord.xy);
   float  d = tex2D(SamplerDepth, IN.txcoord.xy).x;
   float2 cd = IN.txcoord.xy;
   float4 tv;
   float4 wp;  
   float t = GameTime;
   float tf = CalculateGameTime0(t);       
   tv.xy = IN.txcoord.xy*2.0-1.0;
   tv.y = -tv.y;
   tv.z = d;
   tv.w = 1.0;
   wp.x = dot(tv, MatrixInverseVPRotation[0]);
   wp.y = dot(tv, MatrixInverseVPRotation[1]);
   wp.z = dot(tv, MatrixInverseVPRotation[2]);
   wp.w = dot(tv, MatrixInverseVPRotation[3]);
   wp.xyz/= wp.w;

   float4 np0 = float4(normalize(wp.xyz), 1.0);    
   float np1 = (np0.z * 1.2) + 1.0;
	     np1 *= lerp(0.0, 1.0, tf);	

float3 atm = r0;
float atm2 = saturate(np1);
float4 wx = WeatherAndTime;	
float3 scc = 1.0;
float3 ncc = 1.0; 
if (wx.x==0,1) scc = atm;
if (wx.y==0,1) ncc = atm;
if (wx.x==4) scc = atm2;
if (wx.x==7) scc = atm2;
if (wx.x==8) scc = atm2;
if (wx.x==9) scc = atm2;
if (wx.x==12) scc = atm2;
if (wx.x==15) scc = atm2;
if (wx.x==16) scc = atm2;
if (wx.y==4) ncc = atm2;
if (wx.y==7) ncc = atm2;
if (wx.y==8) ncc = atm2;
if (wx.y==9) ncc = atm2;
if (wx.y==12) ncc = atm2;
if (wx.y==15) ncc = atm2;
if (wx.y==16) ncc = atm2;
   float p0 = length(wp.xyz);  
   float3 dp = -1.0/exp(length(0.0-p0)/(25.0*1.30));
   float f1 = pow(tex2D(SamplerDepth, cd).x, 11000.0);
         f1 = lerp(1.0, f1, dp);   
   float4 vf = pow(f1.x, -100.0);
          vf.a = min(1.0, vf.a); 
   float3 wmix = lerp(scc, ncc, wx.z);
   float w0 = saturate(wmix);	  
          r0.xyz = lerp(r0, wmix, vf*0.65);
   return r0;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

technique PostProcess
{
   pass P0
   {
      VertexShader = compile vs_3_0 VS_PostProcess();
      PixelShader  = compile ps_3_0 PS_DX0();
   }
}

