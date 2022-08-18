
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
float4 ScreenSize; float4 Timer;float ENightDayFactor;
float4 SunDirection; float EInteriorFactor; float FadeFactor; float FieldOfView;
float4 MatrixVP[4]; float4 MatrixInverseVP[4]; float4 MatrixVPRotation[4]; float4 MatrixInverseVPRotation[4];
float4 MatrixView[4];float4 MatrixInverseView[4];float4 CameraPosition;float GameTime;float4 CustomShaderConstants1[8];
float4 WeatherAndTime; float4 FogFarColor;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//Textures
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

texture2D texEnv;
texture2D texOriginal;
texture2D texColor;
texture2D texDepth;
texture2D texNoise;
texture2D texShadow;
texture2D texPalette;
texture2D texFocus;
texture2D texCurr; 
texture2D texPrev; 
texture2D texNormal;
texture2D texStars < string ResourceName="stars.png"; >;
texture2D texRain < string ResourceName="Rain.png"; >;
texture2D texRain2 < string ResourceName="Rain2.png"; >;
texture2D texMoon < string ResourceName="Moon.png"; >;
texture2D texNs < string ResourceName="Noise.png"; >;
texture2D texNs2 < string ResourceName="Noise2.png";>;
texture2D texPd < string ResourceName="NoisePd.png"; >;
texture2D texRipples < string ResourceName="Ripples.png"; >;
texture2D texAlpha < string ResourceName="RipplesAlpha.png"; >;
texture2D texAlpha2 < string ResourceName="RipplesAlpha2.png"; >;
texture2D texwrl < string ResourceName="PuddlesRelief.png"; >;
texture2D texltn < string ResourceName="lightning.png"; >;
texture2D texltn2 < string ResourceName="lightning2.png"; >;
textureCUBE texCube < string ResourceName = "SkyLight.dds"; string ResourceType = "CUBE"; >;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//Sampler Inputs
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
sampler2D SamplerLightning2 = sampler_state
{
	Texture = <texltn2>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
    AddressU = Clamp;
    AddressV = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};


sampler2D SamplerLightning = sampler_state
{
	Texture = <texltn>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
    AddressU = Clamp;
    AddressV = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerEnv = sampler_state
{
	Texture   = <texEnv>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU  = Mirror;
	AddressV  = Mirror;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerOriginal = sampler_state
{
	Texture   = <texOriginal>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;//NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerWrl = sampler_state
{
	Texture = <texwrl>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Wrap;
	AddressV = Wrap;
	AddressW = Wrap;
	SRGBTexture=TRUE;
	MaxMipLevel=1;
	MipMapLodBias=0;
};

sampler2D SamplerWrp = sampler_state
{
	Texture = <texRipples>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Wrap;
	AddressV = Wrap;
	AddressW = Wrap;
	SRGBTexture=TRUE;
	MaxMipLevel=1;
	MipMapLodBias=0;
};

sampler2D SamplerAlpha = sampler_state
{
	Texture = <texAlpha>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Wrap;
	AddressV = Wrap;
	AddressW = Wrap;
	SRGBTexture=TRUE;
	MaxMipLevel=1;
	MipMapLodBias=0;
};

sampler2D SamplerAlpha2 = sampler_state
{
	Texture = <texAlpha2>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Wrap;
	AddressV = Wrap;
	AddressW = Wrap;
	SRGBTexture=TRUE;
	MaxMipLevel=1;
	MipMapLodBias=0;
};

samplerCUBE SamplerCM = sampler_state 
{
	Texture = <texCube>;
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

sampler2D SamplerNormal = sampler_state
{
	Texture   = <texNormal>;
	MinFilter = Anisotropic;
	MagFilter = Anisotropic;
	MipFilter = LINEAR;
	AddressU  = Wrap;
	AddressV  = Wrap;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerColor = sampler_state
{
	Texture   = <texColor>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=TRUE;
	MaxMipLevel=2;
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
	SRGBTexture=TRUE;
	MaxMipLevel=1;
	MipMapLodBias=0;
};

sampler2D SamplerNoise = sampler_state
{
	Texture   = <texNoise>;
	MinFilter = POINT;
	MagFilter = POINT;
	MipFilter = NONE;
	AddressU  = Wrap;
	AddressV  = Wrap;
	SRGBTexture=FALSE;
	MaxMipLevel=2;
	MipMapLodBias=0;
};

sampler2D SamplerPalette = sampler_state
{
	Texture   = <texPalette>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};


sampler2D SamplerCurr = sampler_state
{
	Texture   = <texCurr>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;//NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};


sampler2D SamplerPrev = sampler_state
{
	Texture   = <texPrev>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerFocus = sampler_state
{
	Texture   = <texFocus>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerStars = sampler_state
{
   Texture = <texStars>;
   MinFilter = LINEAR;
   MagFilter = LINEAR;
   MipFilter = NONE;
   AddressU = WRAP;
   AddressV = WRAP;
   SRGBTexture=FALSE;
   MaxMipLevel=0;
   MipMapLodBias=0;
}; 

sampler2D SamplerRain = sampler_state
{
   Texture = <texRain>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Wrap;
	AddressV = Wrap;
	AddressW = Wrap;
	SRGBTexture=TRUE;
	MaxMipLevel=0;
	MipMapLodBias=0;
}; 

sampler2D SamplerRain2 = sampler_state
{
   Texture = <texRain2>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Wrap;
	AddressV = Wrap;
	AddressW = Wrap;
	SRGBTexture=TRUE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerMoon = sampler_state
{
   Texture = <texMoon>;
   MinFilter = LINEAR;
   MagFilter = LINEAR;
   MipFilter = NONE;
   AddressU = Clamp;
   AddressV = Clamp;
   SRGBTexture=FALSE;
   MaxMipLevel=0;
   MipMapLodBias=0;
};

sampler2D SamplerNs = sampler_state
{
   Texture   = <texNs>;
   MinFilter = LINEAR;
   MagFilter = LINEAR;
   MipFilter = LINEAR;
   AddressU  = Wrap;
   AddressV  = Wrap;
   AddressW  = Wrap;   
   SRGBTexture=FALSE;
   MaxMipLevel=0;
   MipMapLodBias=0;
};

sampler2D SamplerNs2 = sampler_state
{
   Texture   = <texNs2>;
   MinFilter = LINEAR;
   MagFilter = LINEAR;
   MipFilter = LINEAR;
   AddressU  = Wrap;
   AddressV  = Wrap;
   AddressW  = Wrap;   
   SRGBTexture=FALSE;
   MaxMipLevel=0;
   MipMapLodBias=0;
};

sampler2D SamplerPuddle = sampler_state
{
   Texture   = <texPd>;
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

sampler2D SamplerShadow = sampler_state
{
	Texture   = <texShadow>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
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

float linearlizeDepth(float nonlinearDepth)
{
float2 dofProj = float2(0.0509804, 3098.0392);
float2 dofDist = float2(0.0, 0.0509804);
float4 depth = nonlinearDepth;
       depth.y = -dofProj.x + dofProj.y;
       depth.y = 1.0/depth.y;
       depth.z = depth.y * dofProj.y;
       depth.z = depth.z * -dofProj.x;
       depth.x = dofProj.y * -depth.y + depth.x;
       depth.x = 1.0/depth.x;
       depth.y = depth.z * depth.x;
       depth.x = depth.z * depth.x - dofDist.y;
       depth.x+= dofDist.x * -0.5;
       depth.x = max(depth.x, 0.0);
return depth.x;
}

float CL(in float3 c)
{
	return (c.r * 0.2126 + c.g * 0.7152 + c.b * 0.0722)*1.0;
}

float4 wpd0(in float2 cd, in float d)
{	
float4 tv;
       tv.xy = cd.xy*2.0-1.0;
       tv.y = -tv.y;
       tv.z = d;
       tv.w = 1.0;
float4 wp;
       wp.x = dot(tv, MatrixInverseVPRotation[0]);
       wp.y = dot(tv, MatrixInverseVPRotation[1]);
       wp.z = dot(tv, MatrixInverseVPRotation[2]);
       wp.w = dot(tv, MatrixInverseVPRotation[3]);
       wp.xyz/= wp.w;
       wp.xyz+= CameraPosition;
return wp;
}

float4 wpd1(float2 cd)
{
float4 tv; 
       tv.xy = cd.xy*2.0-1.0;   
       tv.y =-tv.y;   
       tv.z = tex2Dlod(SamplerDepth, float4(cd.xy,0,0)).x;
       tv.w = 1.0;
float4 wp;
	   wp.x = dot(tv,MatrixInverseVPRotation[0]);
	   wp.y = dot(tv,MatrixInverseVPRotation[1]);
	   wp.z = dot(tv,MatrixInverseVPRotation[2]);
	   wp.w = dot(tv,MatrixInverseVPRotation[3]);
	   wp.xyz/= wp.w;	   
return wp;
}

float2 wpd2(float3 cd)
{
float4 tv = float4(cd.xyz, 1.0);
float4 wp = 0.0;
       wp.x = dot(tv,MatrixVPRotation[0]);
       wp.y = dot(tv,MatrixVPRotation[1]);
       wp.z = dot(tv,MatrixVPRotation[2]);
       wp.w = dot(tv,MatrixVPRotation[3]);
       wp.xyz/= wp.w;
       wp.y = -wp.y;
       wp.xy = wp.xy*0.5+0.5;	  
return wp.xy;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
float3 fs0(float2 cd)
{
  float3 wpos = wpd1(cd.xy);
  float2 fscreen = float2(ScreenSize.y,ScreenSize.y*ScreenSize.z);
  float3 y = wpd1(cd.xy+float2(fscreen.x, 0.0)) - wpos;
  float3 f = wpd1(cd.xy+float2(0.0, fscreen.y)) - wpos;
  float3 n = cross(f,y);		 
  return normalize(n);
}

float4 reflection(float3 n, float2 cd)
{
   float3 wpos = wpd1(cd.xy);
   float3 v = float3(0.353, 0.353, 1.0);
          n = (normalize(n.xyz*v)*0.98);
   float3 n0 = reflect(wpos.xyz, n.xyz);	 
   float3 n1 = ((1000.0/0.01)*n0)/1000.0;
   float3 r0 = (wpos+n1);
   float2 r1 = wpd2(r0.xyz);
   float4 r2 = tex2Dlod(SamplerColor, float4(r1.xy, 0.0, 6.0)); 
	      r2.w = r1.y<0.0||r1.y>1.0 ? 0.0:1.0;
   return r2;
}

float Lighting_Brightness
<
	string UIName="Lighting - Brightness";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=6.0;
> = {1.20};

float Lighting_Saturate
<
        string UIName="Lighting - Saturate";
        string UIWidget="Spinner";
        float UIMin=0.0;
        float UIMax=5.0;
> = {1.20};

float Lighting_Contrast
<
	string UIName="Lighting - Contrast";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=5.0;
> = {1.30};

float Shadow_Brightness
<
	string UIName="Shadow - Brightness";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=2.0;
> = {0.65};

float Shadow_Saturate
<
        string UIName="Shadow - Saturate";
        string UIWidget="Spinner";
        float UIMin=0.0;
        float UIMax=5.0;
> = {0.80};

float Shadow_Contrast
<
	string UIName="Shadow - Contrast";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=5.0;
> = {1.80};

float Peds_SideLighting
<
	string UIName="Peds - Fresnel Lighting";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10.0;
> = {2.00};


float Lighting_Level
<
	string UIName="Cars - Shadow";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10.0;
> = {0.60};

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
float4 PS_DX0(VS_OUTPUT_POST IN, float2 vPos : VPOS) : COLOR
{
//=============================================================
   float4 tv0;
   float4 tv1;
   float4 wpos;
   float4 npos; 
   float2 cd = IN.txcoord.xy;
   float4 r0 = tex2D(SamplerColor, cd);	
		  r0.a = min(1.0, r0.a);
   float4 r1 = tex2D(SamplerColor, cd);
   float4 r2 = r1;  
   float4 r3 = r1;  
   float4 r4 = r1;
	      r4.xyz+=0.000001; 
   float4 r5 = tex2D(SamplerOriginal, cd);	  		  
   float4 shadow = tex2D(SamplerShadow, cd);		  
   float3 normal = tex2D(SamplerNormal, IN.txcoord.xy)*2.0-1.0;
   float  d0 = tex2D(SamplerDepth, cd).x;
   float  d1 = linearlizeDepth(tex2D(SamplerDepth, cd).x)*1.0;
//=============================================================	
   float t0 = GameTime;
   float3 sv0 = SunDirection.xyz;
   float3 sv1 = normalize(float3(0.0, 0.0, 0.0));
   float s1 = smoothstep(0.0, 1.0, t0);	
   float s2 = smoothstep(1.0, 23.0, t0);
   float s3 = smoothstep(23.0, 24.0, t0);	
   float3 sunvec = lerp(sv1, sv0, s1);
          sunvec = lerp(sunvec, sv0, s2); 
          sunvec = lerp(sunvec, sv1, s3); 	  
//=============================================================	
   float x1 = smoothstep(0.0, 4.0, t0);
   float x2 = smoothstep(4.0, 5.0, t0);
   float x3 = smoothstep(5.0, 6.0, t0);
   float x4 = smoothstep(6.0, 7.0, t0);
   float xE = smoothstep(8.0, 11.0, t0);
   float x5 = smoothstep(16.0, 17.0, t0);
   float x6 = smoothstep(18.0, 19.0, t0);
   float x7 = smoothstep(19.0, 20.0, t0);
   float xG = smoothstep(20.0, 21.0, t0);  
   float xZ = smoothstep(21.0, 22.0, t0);
   float x8 = smoothstep(22.0, 23.0, t0);
   float x9 = smoothstep(23.0, 24.0, t0); 	   	  
//=============================================================	 		 
   tv1.xy = IN.txcoord.xy*2.0-1.0;
   tv1.y = -tv1.y;
   tv1.z = d0;
   tv1.w = 1.0;
   wpos.x = dot(tv1, MatrixInverseVPRotation[0]);
   wpos.y = dot(tv1, MatrixInverseVPRotation[1]);
   wpos.z = dot(tv1, MatrixInverseVPRotation[2]);
   wpos.w = dot(tv1, MatrixInverseVPRotation[3]);
   wpos.xyz/= wpos.w;
//=============================================================	
   tv0.xy = normal.xy;
   tv0.x = -tv0.x;
   tv0.z = normal.z;
   npos.x = dot(tv0.xyz, MatrixInverseView[0])*1.0;
   npos.y = dot(tv0.xyz, MatrixInverseView[1])*1.0;
   npos.z = dot(tv0.xyz, MatrixInverseView[2])*1.0;
   float3 n2 = normalize(npos.xyz);
   float3 n3 = normalize(wpos.xyz);
//=============================================================	    
   float4 wpos0 = float4(normalize(wpos.xyz), 1.0);	 
   float3 c1 = normalize(n3.xyz-sunvec);
          c1*= smoothstep(16.0, 0.0, normalize(-wpos0.z-sunvec));
		  c1*= saturate(1.0 - dot(r0.xyz, 0.0));
//=============================================================  
   float  d11 = dot(n2, sunvec);
          d11 = d11 * 1.0 + 0.0;		
   float  d12 = pow(d11 * 0.5 + 0.5, Lighting_Level);
		  
   float  d5 = dot(n2.xzy, normalize(float3(0.0, -1.0, 0.0))); 
          d5 = d5 * 1.0 + 0.0;		  
   float  d5s = dot(n2.xzy, normalize(float3(0.0, 1.0, 0.0))); 
          d5s = d5s * 1.0 + 0.0;		  
   float mix3 = pow(d5 * 0.5 + 0.5, 1.5);
   float mix3s = pow(d5s * 0.5 + 0.5, 1.5);  	  
   float4 sl = mix3*4.0; sl.a = min(1.0, sl.a);
   float4 sl2 = mix3s*4.0; sl2.a = min(1.0, sl2.a); 		  
//=============================================================
   float3 timelt = lerp(0.0, 0.00, x1);
          timelt = lerp(timelt, float3(0.686, 0.392, 0.235)*float3(0.745, 0.588, 0.157)*7.0*0.3, x2);
          timelt = lerp(timelt, float3(0.686, 0.392, 0.235)*float3(0.745, 0.588, 0.157)*7.0, x3);
          timelt = lerp(timelt, 4.1*float3(1.0, 0.745, 0.51), x4);
          timelt = lerp(timelt, 4.1*float3(1.0, 0.98, 0.745)*0.90, xE);
          timelt = lerp(timelt, 4.1*float3(1.0, 0.98, 0.745)*0.90, x5);
          timelt = lerp(timelt, float3(0.686, 0.392, 0.235)*float3(0.745, 0.588, 0.157)*7.0, x6);		 
          timelt = lerp(timelt, 0.0, x7);
		  timelt = lerp(timelt, 0.00, xG);
		  timelt = lerp(timelt, 0.00, xZ);
          timelt = lerp(timelt, 0.00, x8);
          timelt*= 1.0;		  
   float3 timesh = lerp(0.0, 0.00, x1);
          timesh = lerp(timesh, 0.10, x2); 
          timesh = lerp(timesh, 0.10, x3);
          timesh = lerp(timesh, 0.20, x4);
          timesh = lerp(timesh, 0.2, xE);
          timesh = lerp(timesh, 0.2, x5);
          timesh = lerp(timesh, 0.40, x6);		 
          timesh = lerp(timesh, 0.10, x7);
		  timesh = lerp(timesh, 0.10, xG);
		  timesh = lerp(timesh, 0.10, xZ);
          timesh = lerp(timesh, 0.10, x8);  
          timesh = lerp(timesh, 0.0, x9); 		  	  
//=============================================================	
	float4 sh0 = shadow; 
	float4 sh1 = shadow;	
	float4 sh2 = shadow;	
    sh0*= saturate(0.0-(d11) * 2.0);
	sh0.xyz = lerp(1.0, sh0, timesh);
    sh1*= saturate(0.0-(d11) * 2.0);
	sh1.xyz = lerp(1.0, sh1, timesh*Lighting_Level);
	sh2*= saturate(-d11);
	float3 albedo = r0;
	       albedo = clamp(albedo, 0.001, 1.0);
	       albedo.xyz = albedo.xyz/(0.5+dot(albedo.xyz, 0.35));
		   
	float3 st0 = normalize(albedo.xyz);
	float3 cs0=albedo.xyz/st0.xyz;
	       cs0=pow(cs0, Shadow_Contrast);
		   st0.xyz = pow(st0.xyz, Shadow_Saturate);
	       albedo.xyz = cs0*st0.xyz;			   
//=============================================================	   
	   float4 mask0 = tex2D(SamplerNormal, cd).w;	
	if (mask0.w<0.99) mask0 = 0.0;	
	   float4 mask1 = tex2D(SamplerNormal, cd).w;	
	if (mask1.w<0.99) mask1 = 50.0*1.9;
		              mask1*= 0.01;
	   float4 mask2 = tex2D(SamplerNormal, cd).w;	
	if (mask2.w<0.0) mask2 = 0.0;			
	   float4 mask3 = tex2D(SamplerNormal, cd).w;
	if (mask3.w==254/255.0) mask3 = 0.0;			
	   float4 mask4 = tex2D(SamplerNormal, cd).w;
	if (mask4.w==254/255.0) mask4 = 50.0*1.9;	
		                    mask4*= 0.01;			
	   float4 mask5 = tex2D(SamplerNormal, cd).w;
	if (mask5.w==253/255.0) mask5 = 5.0*100.0; 
		                    mask5*= 0.001;	
  float4 mask6 = r5; 
	     mask6.xyz+= 0.000001;
  float3 st1 = normalize(mask6.xyz);
  float3 cs1 = mask6.xyz/st1.xyz;
	     cs1 = pow(cs1, 1.0);
		 st1.xyz = pow(st1.xyz, 1.0);
	     mask6.xyz = cs1*st1.xyz;		 
  float4 mask7 = r5; 
	     mask7.xyz+= 0.000001;
  float3 st2 = normalize(mask7.xyz);
  float3 cs2 = mask7.xyz/st2.xyz;
	     cs2 = pow(cs2, 2.0);
		 st2.xyz = pow(st2.xyz, 0.1);
	     mask7.xyz = cs2*st2.xyz;
	float graymask0 = saturate(mask6);	
	float graymask1 = saturate(mask7);								
//=============================================================
   float3 timesk = lerp(float3(0.843, 0.784, 0.667)*0.03, float3(0.843, 0.784, 0.667)*0.03, x1);
          timesk = lerp(timesk, float3(0.843, 0.784, 0.667)*0.6, x2);
          timesk = lerp(timesk, float3(0.843, 0.784, 0.667), x3);
          timesk = lerp(timesk, float3(0.843, 0.784, 0.667), x4);
          timesk = lerp(timesk, float3(0.843, 0.784, 0.667)*1.2, xE);
          timesk = lerp(timesk, float3(0.843, 0.784, 0.667)*1.2, x5);
          timesk = lerp(timesk, float3(0.843, 0.784, 0.667), x6); 
          timesk = lerp(timesk, float3(0.843, 0.784, 0.667)*0.1, x7);
		  timesk = lerp(timesk, float3(0.843, 0.784, 0.667)*0.05, xG);
		  timesk = lerp(timesk, float3(0.843, 0.784, 0.667)*0.05, xZ);
          timesk = lerp(timesk, float3(0.843, 0.784, 0.667)*0.03, x8);		  	  
   float3 timesk2 = lerp(float3(0.314, 0.431, 0.588)*1.0, float3(0.314, 0.431, 0.588)*1.0, x1);
          timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588)*0.6, x2);
          timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588), x3);
          timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588), x4);
          timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588), xE);
          timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588), x5);
          timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588), x6);	 
          timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588)*0.6, x7);
		  timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588)*0.6, xG);
		  timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588)*0.7, xZ);
          timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588)*0.7, x8);	  
          timesk2 = lerp(timesk2, float3(0.314, 0.431, 0.588)*1.0, x8);		  
//=============================================================	 
float4 wx = WeatherAndTime;
float4 lc;
float4 ln;
float3 slc;
float3 sln;

	   float4 mask0z = tex2D(SamplerNormal, cd).w;	
	if (mask0z.w<0.50) mask0z = 0.0;	
	   float4 mask1z = tex2D(SamplerNormal, cd).w;	
	if (mask1z.w<0.99) mask1z = 50.0*40.0;
		              mask1z*= 0.001;
	float mask2z = lerp(1.0, 0.0, mask1z);				  
          mask2z = saturate(mask2z);					  
	float mask1x = saturate(mask1z*mask0z)*1.0;	

if (wx.x==0,1) slc = lerp(((timesk2*sl)+(sl2*timesk))*1.0, 1.0, mask1x);
if (wx.y==0,1) sln = lerp(((timesk2*sl)+(sl2*timesk))*1.0, 1.0, mask1x);
if (wx.x==4) slc = 1.0;
if (wx.x==7) slc = 1.0;
if (wx.x==8) slc = 1.0;
if (wx.x==9) slc = 1.0;
if (wx.x==12) slc = 1.0;
if (wx.x==15) slc = 1.0;
if (wx.x==16) slc = 1.0;
if (wx.y==4) sln = 1.0;
if (wx.y==7) sln = 1.0;
if (wx.y==8) sln = 1.0;
if (wx.y==9) sln = 1.0;
if (wx.y==12) sln = 1.0;
if (wx.y==15) sln = 1.0;
if (wx.y==16) sln = 1.0;
    float3 wmix = lerp(slc, sln, pow(wx.z,5.0));
	
	
	
	
	
	if(d1 < 3000)r0.xyz*= wmix*sh0*albedo*Shadow_Brightness;
	if(d1 < 3000)r0.xyz = lerp(r1, r0, mask0z); 
	if(d1 < 3000)r0.xyz = lerp(r1, r0, mask2z);	
	if(d1 < 3000)r1.xyz*= wmix*sh1*Shadow_Brightness;
    if(d1 < 3000)r0.xyz = lerp(r0.xyz, r1.xyz, mask1x);
	if(d1 < 3000)r0.xyz = lerp(r2.xyz, r0.xyz, mask2);	
	
	float3 st3 = normalize(r2.xyz);
	float3 cs3=r2.xyz/st3.xyz;
	       cs3=pow(cs3, Lighting_Contrast);
		   st3.xyz = pow(st3.xyz, Lighting_Saturate);
	       r2.xyz = cs3*st3.xyz;	
 
   float cz = 100.0*0.05;
   float f0 = pow(dot(c1.z, n2.z), 0.15);
         f0 = pow(f0, cz);
   float f1 = (f0*f0);
         f1/= 2.5;		 
		
   float3 timebl = lerp(0.0, 0.00, x1);
          timebl = lerp(timebl, float3(0.667, 0.51, 0.392)*float3(0.706, 0.51, 0.157)*7.0, x2);
          timebl = lerp(timebl, float3(0.667, 0.51, 0.392)*float3(0.706, 0.51, 0.157)*7.0*8.0, x3);
          timebl = lerp(timebl, 4.1*float3(1.0, 0.745, 0.51)*2.0, x4);
          timebl = lerp(timebl, 4.1*float3(1.0, 0.98, 0.745)*0.30, xE);
          timebl = lerp(timebl, 4.1*float3(1.0, 0.98, 0.745)*0.30, x5);
          timebl = lerp(timebl, float3(0.667, 0.51, 0.392)*float3(0.706, 0.51, 0.157)*7.0*8.0, x6);	 
          timebl = lerp(timebl, 0.0, x7);
		  timebl = lerp(timebl, 0.00, xG);
		  timebl = lerp(timebl, 0.00, xZ); 
          timebl = lerp(timebl, 0.00, x8);
          timebl*= 1.0;	

	if(d1 < 3000)r0.xyz+= lerp(r2*sh2*graymask0*timelt*(mask3*mask0)*Lighting_Brightness, graymask1*10.0*f1*timebl*Lighting_Brightness*0.80, f1*graymask0*sh2*mask2z);
	 if(d1 < 3000)r0.xyz+= r3*sh2*(timelt*0.0)*(mask1+mask4)*Lighting_Brightness;
	 if(d1 < 3000)r0.xyz = lerp(r3, r0, mask0z);
//=============================================================
	float3 wp = reflect(n3.xyz, n2.xyz);
    float4 pr = texCUBE(SamplerCM, wp.xzy);  
	float4 n0 = tex2D(SamplerNormal, cd);    
	float smix;		 
	      smix=saturate(npos.z)*saturate(-n0.z);
	      smix*=smix;	
	float nf1 = saturate(1.65+smix*0.3 -abs(n0.z));
	      nf1 = pow(nf1, 10.0);			 
   r0.xyz = lerp(r0, r0+pr*Peds_SideLighting, r2*nf1*mask5);
//=============================================================   
  float3 st4 = normalize(r4.xyz);
  float3 cs4=r4.xyz/st4.xyz;
	     cs4=pow(cs4, 2.1);
		 st4.xyz = pow(st4.xyz, 1.6);
	     r4.xyz = lerp(r4*0.6, cs4*st4.xyz, mask0*mask3);
   float4 t2;		 
          t2 = lerp(r4*1.5, r4*1.5, x1);
          t2 = lerp(t2, r0, x2);
          t2 = lerp(t2, r0, x3);
          t2 = lerp(t2, r0, x4);
          t2 = lerp(t2, r0, xE);
          t2 = lerp(t2, r0, x5);
          t2 = lerp(t2, r0, x6);		 
          t2 = lerp(t2, r0, x7);
          t2 = lerp(t2, r4*1.5, x8);
if (wx.x==0,1) lc = t2;
if (wx.y==0,1) ln = t2;
if (wx.x==4) lc = r4;
if (wx.x==7) lc = r4;
if (wx.x==8) lc = r4;
if (wx.x==9) lc = r4;
if (wx.x==12) lc = r4;
if (wx.x==15) lc = r4;
if (wx.x==16) lc = r4;
if (wx.y==4) ln = r4;
if (wx.y==7) ln = r4;
if (wx.y==8) ln = r4;
if (wx.y==9) ln = r4;
if (wx.y==12) ln = r4;
if (wx.y==15) ln = r4;
if (wx.y==16) ln = r4;
   float4 wmix0 = lerp(lc, ln, pow(wx.z, 5.0)); 
//=============================================================     
  //r0.xyz = mask2;
   //return r0;
   return wmix0;
   
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

float4 	WorldPos(in float2 coord)
{
   float  d0 = tex2D(SamplerDepth, coord.xy).x;
   float4 tvec; 
          tvec.xy = coord.xy*2.0-1.0;   
          tvec.y = -tvec.y;   
          tvec.z = d0;
          tvec.w = 1.0;
   float4 wpos;
          wpos.x = dot(tvec, MatrixInverseVPRotation[0]);
          wpos.y = dot(tvec, MatrixInverseVPRotation[1]);
          wpos.z = dot(tvec, MatrixInverseVPRotation[2]);
          wpos.w = dot(tvec, MatrixInverseVPRotation[3]);
          wpos.xyz/= wpos.w;
   return wpos;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// etc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

float	Sky_Saturate
<
	string UIName="Sky - Saturate";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=6.0;
> = {1.30};

float	Sky_Contrast
<
	string UIName="Sky - Contrast";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=6.0;
> = {1.35};

float	Sky_Brightness
<
	string UIName="Sky - Brightness";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=6.0;
> = {1.00};


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
   float3 f3 = lerp(0.0, (f0*t4), smoothstep(0.0, 0.005, np1.z)) +(f1*t0)+(f2*t3)+(fnight*t5);
 
float3 SunCurrent;
float3 SunNext;
float4 wx = WeatherAndTime;

if (wx.x==0,1) SunCurrent = f3;
if (wx.y==0,1) SunNext = f3;

if (wx.x==4) SunCurrent = 0.0;
if (wx.x==7) SunCurrent = 0.0;
if (wx.x==8) SunCurrent = 0.0;
if (wx.x==9) SunCurrent = 0.0;
if (wx.x==12) SunCurrent = 0.0;
if (wx.x==15) SunCurrent = 0.0;
if (wx.x==16) SunCurrent = 0.0;

if (wx.y==4) SunNext = 0.0;
if (wx.y==7) SunNext = 0.0;
if (wx.y==8) SunNext = 0.0;
if (wx.y==9) SunNext = 0.0;
if (wx.y==12) SunNext = 0.0;
if (wx.y==15) SunNext = 0.0;
if (wx.y==16) SunNext = 0.0;	  

float3 wmix = lerp(SunCurrent, SunNext, wx.z);	   
 
   return wmix;
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
   
float4 wx = WeatherAndTime;

if (wx.x==0,1) cc = float3(0.251, 0.51, 1.0);
if (wx.y==0,1) cn = float3(0.251, 0.51, 1.0);

float3 skyW = float3(0.259, 0.259, 0.259);

if (wx.x==4) cc = skyW;
if (wx.x==7) cc = skyW;
if (wx.x==8) cc = skyW;
if (wx.x==9) cc = skyW;
if (wx.x==12) cc = skyW;
if (wx.x==15) cc = skyW;
if (wx.x==16) cc = skyW;

if (wx.y==4) cn = skyW;
if (wx.y==7) cn = skyW;
if (wx.y==8) cn = skyW;
if (wx.y==9) cn = skyW;
if (wx.y==12) cn = skyW;
if (wx.y==15) cn = skyW;
if (wx.y==16) cn = skyW;

float3 wmix0 = lerp(cc, cn, wx.z);

return wmix0;
}

float Noise1(in float3 p)
{	
    float t = (Timer.x * 5.0) * 30.0;
    p.z += 0.5 * t;
	p.xyz += 1.0 * t;
    float2 uv = (p * float2(17.2, 17.2));	
	float2 c = (uv  + 0.5) / 20.0;
    float r0 = tex2Dlod(SamplerNs, float4(c*0.04, 0.0, 0.0)).x;
	return r0;
}

float Noise2(in float3 p)
{	
    float t = (Timer.x * 30.0) * 30.0;
    p.z += 0.5 * t;
	p.xyz += 1.0 * t;
    float2 uv = (p * float2(17.2, 17.2));
	float2 c = (uv + 0.5) / 20.0;
    float r0 = tex2Dlod(SamplerNs2, float4(c*0.04, 0.0, 0.0)).x;
	return r0;
}

float4 	GWP(in float2 c, in float d)
{
   float4 t; 
   float4 w;
   t.xy=c.xy*2.0-1.0;
   t.y=-t.y;
   t.z=d;
   t.w=1.0;
   w.x=dot(t, MatrixInverseVPRotation[0]);
   w.y=dot(t, MatrixInverseVPRotation[1]);
   w.z=dot(t, MatrixInverseVPRotation[2]);
   w.w=dot(t, MatrixInverseVPRotation[3]);
   w.xyz/=w.w;
   //w.xyz+=CameraPosition;
   return w;
}

float SunLight(in float2 c, in float d)
{
    float3 v0=SunDirection.xyz;
    float3 v2=normalize(float3(-0.09, -0.94, 0.31));	
   float t0 = GameTime;
   float x1 = smoothstep(0.0, 4.0, t0);	
   float x2 = smoothstep(4.0, 23.0, t0);
   float x3 = smoothstep(23.0, 24.0, t0);	
   float3 sv = lerp(v2, v0, x1);
          sv = lerp(sv, v0, x2);
          sv = lerp(sv, v2, x3);	
   float4 WorldP = GWP(c, d);
   float c0 = 6.07;
   float3 np = normalize(WorldP.xyz);
   float f0 = 0.95 - dot(-sv, np);
         f0 = pow(f0, c0);
   float r0 = (f0*f0);
   return r0/40000.0;
}

float Coverage(in float v, in float d, in float c)
{
	c = clamp(c - (1.0 - v), 0.0, 1.0 - d)/(1.0 - d);
	c = max(0.0, c * 1.1 - 0.1);
	c = c = c * c * (3.0 - 2.0 * c);
	return c;
}

float fcc
<
        string UIName="Clouds - Coverage";
        string UIWidget="Spinner";
        float UIMin=0.0;
        float UIMax=0.50;
> = {0.40};

float4 GenerateClouds(float4 worldpos, in float3 sunlight)
{
	float t = (Timer.x * 1000.0)*30.0;	
	float3 vecm = normalize(float3(-0.09, -0.94, 0.31));
    float3 sv0=-SunDirection.xyz;
    float3 sv2=-vecm;
	
   float t0 = GameTime;
   float x1 = smoothstep(0.0, 4.0, t0);
   float x2 = smoothstep(4.0, 5.0, t0);
   float x3 = smoothstep(5.0, 6.0, t0);
   float x4 = smoothstep(6.0, 7.0, t0);
   float xE = smoothstep(8.0, 11.0, t0);
   float x5 = smoothstep(16.0, 17.0, t0);
   float x6 = smoothstep(18.0, 19.0, t0);
   float x7 = smoothstep(19.0, 20.0, t0);
   float xG = smoothstep(20.0, 21.0, t0);  
   float xZ = smoothstep(21.0, 22.0, t0);
   float x8 = smoothstep(22.0, 23.0, t0);
   float x9 = smoothstep(23.0, 24.0, t0);
   float x10 = smoothstep(4.0, 23.0, t0);
   float x11 = smoothstep(23.0, 24.0, t0);	
   
   float3 sv = lerp(sv2, sv0, x1);
          sv = lerp(sv, sv0, x10);
          sv = lerp(sv, sv2, x11);	  	
	
	float3 wp = worldpos/1.5;
		   wp.x *= 0.04;
		   wp.y *= 0.10;	
		   wp.y -= t * 0.01;	
	float3 wp1 = wp * float3(1.0, 0.5, 1.0) + float3(0.0, t * 0.01, 0.0);
	
	float noise  = 	Noise2(wp * float3(1.0, 0.5, 1.0) + float3(0.0, t * 0.01, 0.0));
		   wp *= 3.0;
		   wp.xy -= t * 0.04;
		   wp.x += 2.0;
	float3 wp2 = wp;
	
		  noise += (2.0 - abs(Noise2(wp) * 0.8)) * 0.25;
		   wp *= 10.0;
		   wp.xy -= t * 0.035;
	float3 wp3 = wp;
		  
          noise += (2.0 - abs(Noise1(wp) * 1.0)) * 0.03;
		  noise /= 0.80;
		  
float4 wx = WeatherAndTime;
float CovCurrent;
float CovNext;
if (wx.x==0,1) CovCurrent = fcc;
if (wx.y==0,1) CovNext = fcc;

if (wx.x==0) CovCurrent = 0.00;
if (wx.x==2) CovCurrent = 0.00;
if (wx.x==6) CovCurrent = 0.00;
if (wx.x==11) CovCurrent = 0.00;
if (wx.x==13) CovCurrent = 0.00;  
if (wx.x==17) CovCurrent = 0.00; 

if (wx.x==4) CovCurrent = 1.70;
if (wx.x==7) CovCurrent = 1.70;
if (wx.x==8) CovCurrent = 1.70;
if (wx.x==9) CovCurrent = 1.70;
if (wx.x==12) CovCurrent = 1.70;
if (wx.x==15) CovCurrent = 1.70;
if (wx.x==16) CovCurrent = 1.70;

if (wx.y==0) CovNext = 0.00;
if (wx.y==2) CovNext = 0.00;
if (wx.y==6) CovNext = 0.00;
if (wx.y==11) CovNext = 0.00;
if (wx.y==13) CovNext = 0.00;
if (wx.y==17) CovNext = 0.00;

if (wx.y==4) CovNext = 1.70;
if (wx.y==7) CovNext = 1.70;
if (wx.y==8) CovNext = 1.70;
if (wx.y==9) CovNext = 1.70;
if (wx.y==12) CovNext = 1.70;
if (wx.y==15) CovNext = 1.70;
if (wx.y==16) CovNext = 1.70;	

float3 wmix0 = lerp(CovCurrent, CovNext, pow(wx.z, 0.9));
		  	  
	float coverage = wmix0;	  	
	float dn = 0.1 - 0.3 * 0.3;  
	noise = Coverage(coverage, dn, noise);
	
	float d0 = Noise2(wp1 + sv.xyz * 0.70 / 2.3);
		  d0 += (2.0 - abs(Noise2(wp2 + sv.xyz * 0.70 / 2.3) * 0.8)) * 0.25;	
		  d0 += (2.0 - abs(Noise1(wp3 + sv.xyz * 0.70 / 2.3) * 1.0)) * 0.03; 
	      d0 = Coverage(0.84, dn, d0);
	
	float bf = lerp(clamp(pow(noise, 0.5) * 1.0, 0.0, 1.0), 0.5, pow(sunlight, 1.0));
		  d0 *= bf;	
			  
   float3 lt = lerp(float3(0.0275, 0.0353, 0.0471), 0.01, x1);
          lt = lerp(lt, float3(1.0, 0.627, 0.235)*0.6, x2);
          lt = lerp(lt, float3(1.0, 0.627, 0.235), x3);
          lt = lerp(lt, float3(1.0, 0.627, 0.235), x4);
          lt = lerp(lt, float3(1.0, 1.0, 1.0), xE);
          lt = lerp(lt, float3(1.0, 1.0, 1.0), x5);
          lt = lerp(lt, float3(1.0, 0.627, 0.235), x6);		 
          lt = lerp(lt, float3(1.0, 0.627, 0.235), x7);
		  lt = lerp(lt, 0.5, xG);
		  lt = lerp(lt, 0.3, xZ);
          lt = lerp(lt, 0.1, x8);
          lt = lerp(lt, float3(0.0275, 0.0353, 0.0471), x9);

   float3 sh = lerp(float3(0.051, 0.0784, 0.118)*0.05, float3(0.03, 0.04, 0.05)*0.2, x1);
          sh = lerp(sh, float3(0.157, 0.165, 0.216)*0.6, x2);
          sh = lerp(sh, float3(0.157, 0.165, 0.216), x3);
          sh = lerp(sh, float3(0.157, 0.165, 0.216), x4);
          sh = lerp(sh, float3(0.392, 0.392, 0.392), xE);
          sh = lerp(sh, float3(0.392, 0.392, 0.392), x5);
          sh = lerp(sh, float3(0.157, 0.165, 0.216), x6);	 
          sh = lerp(sh, float3(0.157, 0.165, 0.216), x7);
		  sh = lerp(sh, float3(0.157, 0.165, 0.216), xG);
		  sh = lerp(sh, float3(0.157, 0.165, 0.216)*0.6, xZ);
          sh = lerp(sh, float3(0.157, 0.165, 0.216)*0.4, x8);
          sh = lerp(sh, float3(0.051, 0.0784, 0.118)*0.05, x9);	

   float3 fmix = lerp(0.02, 0.1, x1);
          fmix = lerp(fmix, 0.1, x2);
          fmix = lerp(fmix, 0.2, x3);			  
          fmix = lerp(fmix, 0.5, x4);
          fmix = lerp(fmix, 0.5, xE);
          fmix = lerp(fmix, 0.5, x5);
          fmix = lerp(fmix, 0.3, x6);
          fmix = lerp(fmix, 0.2, x7); 
          fmix = lerp(fmix, 0.1, xG);	
          fmix = lerp(fmix, 0.05, xZ);
          fmix = lerp(fmix, 0.03, x8);	  
          fmix = lerp(fmix, 0.02, x9);		  
		  
   float3 fmix2 = lerp(0.01, 0.1, x1);
          fmix2 = lerp(fmix2, 0.05, x2);
          fmix2 = lerp(fmix2, 0.1, x3);			  
          fmix2 = lerp(fmix2, 0.2, x4);
          fmix2 = lerp(fmix2, 0.2, xE);
          fmix2 = lerp(fmix2, 0.2, x5);
          fmix2 = lerp(fmix2, 0.1, x6);
          fmix2 = lerp(fmix2, 0.08, x7); 
          fmix2 = lerp(fmix2, 0.05, xG);	
          fmix2 = lerp(fmix2, 0.03, xZ);
          fmix2 = lerp(fmix2, 0.02, x8);	  
          fmix2 = lerp(fmix2, 0.01, x9);
		  
   float3 ColorSun = lerp(float3(0.0392, 0.0588, 0.0784)*0.2, float3(0.05, 0.04, 0.02)*0.0, x1);
          ColorSun = lerp(ColorSun, float3(0.40, 0.08, 0.0)*0.8, x2);
          ColorSun = lerp(ColorSun, float3(0.45, 0.07, 0.0), x3);
          ColorSun = lerp(ColorSun, float3(0.392, 0.118, 0.0588), x4);
          ColorSun = lerp(ColorSun, float3(0.392, 0.275, 0.118), xE);
          ColorSun = lerp(ColorSun, float3(0.392, 0.275, 0.118), x5);
          ColorSun = lerp(ColorSun, float3(0.392, 0.08, 0.0), x6); 
          ColorSun = lerp(ColorSun, float3(0.392, 0.08, 0.0), x7);
		  ColorSun = lerp(ColorSun, float3(0.392, 0.118, 0.0588), xG);
		  ColorSun = lerp(ColorSun, float3(0.392, 0.118, 0.0588)*0.4, xZ);
          ColorSun = lerp(ColorSun, float3(0.392, 0.118, 0.0588)*0.2, x8);		  
          ColorSun = lerp(ColorSun, float3(0.0392, 0.0588, 0.0784)*0.2, x9);
		  
float3 SunCurrent;
float3 SunNext;

 float3 m0 = float3(min(1.0, d0), min(1.0, d0), min(1.0, d0));
 float3 color = lerp((pow(sunlight, 0.3)*ColorSun*50.0)+lt*2.0, sh, m0);
 float3 colorPasm = lerp(fmix, fmix2, m0);
 
if (wx.x==0,1) SunCurrent = color;
if (wx.y==0,1) SunNext = color;

if (wx.x==4) SunCurrent = colorPasm;
if (wx.x==7) SunCurrent = colorPasm;
if (wx.x==8) SunCurrent = colorPasm;
if (wx.x==9) SunCurrent = colorPasm;
if (wx.x==12) SunCurrent = colorPasm;
if (wx.x==15) SunCurrent = colorPasm;
if (wx.x==16) SunCurrent = colorPasm;

if (wx.y==4) SunNext = colorPasm;
if (wx.y==7) SunNext = colorPasm;
if (wx.y==8) SunNext = colorPasm;
if (wx.y==9) SunNext = colorPasm;
if (wx.y==12) SunNext = colorPasm;
if (wx.y==15) SunNext = colorPasm;
if (wx.y==16) SunNext = colorPasm;	  

float3 wmix = lerp(SunCurrent, SunNext, wx.z);	
		  
	float4 r0 = float4(wmix.rgb, (noise*noise*noise));
	return r0;
}

float IntersectPlane(float3 pos, float3 dir)
{
	return -pos.z/dir.z;
}

float Stars_Brightness
<
        string UIName="Stars - Brightness";
        string UIWidget="Spinner";
        float UIMin=0.0;
        float UIMax=5.0;
> = {1.00};

float starscover
<
        string UIName="Stars - Size";
        string UIWidget="Spinner";
        float UIMin=0.0;
        float UIMax=0.20;
> = {0.05};

float Stars(in float2 coord)
{	
  //float t = (Timer.x * 10.0) * 0.25;
    float2 uv = coord.xy;
	float2 p = (uv) / (starscover*2.0);
	float fstars = tex2Dlod(SamplerStars, float4(p, 0.0, 0.0));
	return fstars;
}

float Moon(in float2 coord)
{	
    float2 uv = coord.xy + 0.1 * 1.5;
	float2 p = (uv + 0.5) / (-3.25*2.0);
	float fmoon = tex2Dlod(SamplerMoon, float4(p, 0.0, 0.0));
	return fmoon*1.5;
}


////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
float4 PS_DX1(VS_OUTPUT_POST IN, float2 vPos : VPOS) : COLOR
{
   float t = GameTime;
   float tf = CalculateGameTime(t);   
   float tf2 = CalculateGameTime0(t);    
   
   float2 coord = IN.txcoord.xy;
   float4 r0 = tex2D(SamplerColor, coord.xy);
   float  d0 = tex2D(SamplerDepth, coord.xy).x;
   float4 wpos = WorldPos(coord);
   float3 colorS = ColorTop();  
//Sky------------------------------------------------------------------------------	
   float4 w0 = wpos+CameraPosition;
   float4 cam = CameraPosition; 
   float p2 = length(wpos.xyz);  
   float3 dp = -1.0/exp(length(0.0-p2)/250.0);  
   float fog = pow(tex2D(SamplerDepth, coord).x, 117000.0);
         fog = lerp(1.0, fog, dp); 		 
   float4 f0 = pow(fog.x, -1000.0);
          f0.a = min(1.0, f0.a); 	  
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
 		   d1 = CL(d1) * colorS;
		   d1 *= (sgf * (20.0*1.4)) + 0.85;
		   
   float3 sun = CalculateSun(coord);		   
   float3 d2 = ColorHorizon(d1, sgf);
   
	d1 = d2;
	d1*= lerp(0.0, 1.0, tf);
	d1.xyz+=0.000001;
	float3 n0 = normalize(d1.xyz);
	float3 ct0=d1.xyz/n0.xyz;
	ct0=pow(ct0, Sky_Contrast);
    n0.xyz = pow(n0.xyz, Sky_Saturate);   
    d1.xyz = ct0*n0.xyz;  
    d1.xyz*= Sky_Brightness;
	
    float np4 = (np0.z * (2.0)) + 0.40;
	      np4 *= lerp(0.0, 1.0, tf2);  		  
    float3 a = d1+sun;
    float3 a2 = np4;
    float3   sc = 1.0;
    float3   sn = 1.0;   
    float4 wx = WeatherAndTime;	
    if (wx.x==0,1) sc = a;
    if (wx.y==0,1) sn = a;   
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
    float3 wmix0 = lerp(sc, sn, wx.z);			
    r0.xyz = lerp(r0, wmix0, f0);	
   float4  worldpos = GWP(coord, d0);
   float L = worldpos.xyz;
//NightSky-------------------------------------------------------------------------
   float sm1 = smoothstep(0.0, 4.0, t),	
         sm2 = smoothstep(4.0, 5.0, t),	 
         sm3 = smoothstep(5.0, 6.0, t),			 
         sm4 = smoothstep(6.0, 23.0, t),
         sm5 = smoothstep(23.0, 24.0, t);
   float3 ti = lerp(1.0, 1.0, sm1);
          ti = lerp(ti, 0.3, sm2);   
          ti = lerp(ti, 0.0, sm3); 		  
          ti = lerp(ti, 0.0, sm4);
          ti = lerp(ti, 1.0, sm5);    
   float3 p0 = (MatrixInverseVP[3].xyz/MatrixInverseVP[3].w);		  
   float3 ns0 = normalize(worldpos.xzy);
   float3 ns1 = normalize(worldpos.xyz);
          ns1.z = abs(ns1.z)*1.40;	
   float3x3 vec0 = float3x3(1, 0, 0, 0, cos(2.91), sin(2.91), 0, -sin(2.91), cos(2.91));   
		   ns0 = mul(ns0, vec0);		   
	float3 np1 = normalize(ns0);
    float3 np2 = normalize(-p0+ns1);
    float4 np3 = float4(normalize(worldpos.xyz), 1.0); 		
    if( 10000 < L || d0 == 1)
      {
          float4 s0 = Stars(np2 * 0.5 + 0.5)*3.2*Stars_Brightness;
		         s0*= step(0.3, s0);
				 s0*= smoothstep(0.0, 0.3, np3.z);
			     s0*= saturate(1.0 - dot(r0.xyz, 4.0));
				 s0.a = min(1.0, ti*s0.a);
                 r0.xyz = lerp(r0, s0.xyz, s0.a);
      }
    if( 6000 < L || d0 == 1)
      {
          float4 m0 = Moon(np1*47.5)*0.78;
		  		 m0*= smoothstep(0.0, 0.3, -np3.y);
				 m0.a = min(1.0, ti*m0.a);
                 r0.xyz = lerp(r0, 1.5*m0.xyz, m0.a);
      }
//Clouds---------------------------------------------------------------------------
   float d3 = pow(tex2D(SamplerDepth, coord).x, 17000.0);
   float3 sc0 = SunLight(coord, d3);
   float3 w1 = float3(1.0, 1.0, 10.0);
   float4 ns = float4(normalize(worldpos.xyz), 1.0);    
   float ip = IntersectPlane(w1, ns.xyz);
   float4 r1 = r0;
   if (ip <= 1.0)
   if( 10000.0 < L || d0 == 1)
      {	
        float4 c0 = GenerateClouds(float4((ns.xyz*ip*5.0), 1.0), sc0);	
		   	   c0.a = min(1.0, c0.a);
        r1.xyz = lerp(r0, c0.rgb*3.0, c0.a);
       }
        r0.xyz = lerp(r0, r1, smoothstep(0.0, 0.10, pow(ns.z, 1.85)));	
//---------------------------------------------------------------------------------		  
	return r0;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

float Rain0(in float2 cd)
{	
    float2 uv = cd.yx;
	float2 p = uv/40.0;
	float r = tex2Dlod(SamplerRain2, float4(p*float2(10.0, 2.0), 0.0, 0.0)).y;
	return r;
}

float Rain1(in float2 cd)
{	
    float2 uv = cd.yx;
	float2 p = uv/1.2;
	float r = tex2Dlod(SamplerRain, float4(p*float2(3.5, 0.3), 0.0, 0.0)).y;	
	return r;
}

float Rain2(in float2 cd)
{	
    float2 uv = cd.yx;
	float2 p = uv/0.5;
	float r = tex2Dlod(SamplerRain, float4(p*float2(3.0, 0.5), 0.0, 0.0)).y;
	float r2 = tex2Dlod(SamplerRain2, float4(p*float2(2.0, 0.5), 0.0, 0.0)).y;
	float j = (r+r2);
	return j;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

float3 ComputeRipple(float2 uv, float CurrentTime, float Weight)
{
//https://seblagarde.wordpress.com/2013/01/03/water-drop-2b-dynamic-rain-and-its-effects/	
    float4 rp = tex2Dlod(SamplerWrp, float4(uv, 0.0, 0.0))*1.4;
           rp.yz = rp.yz * 2.0-1.0;	
    float a0 = tex2Dlod(SamplerAlpha, float4(uv, 0.0, 0.0))*1.4;
    float a1 = tex2Dlod(SamplerAlpha2, float4(uv, 0.0, 0.0))*1.4;	

	float4 rp3 = lerp(0.0, rp, a0);
		   
float4  wx = WeatherAndTime;	   
float  wc;
float  wn;
if (wx.x==0,1) wc = 0.0;
if (wx.y==0,1) wn = 0.0;
if (wx.x==8)   wc = 2.0;
if (wx.x==16)  wc = 2.0;
if (wx.y==8)   wn = 2.0;
if (wx.y==16)  wn = 2.0;	  
      float rp0 = lerp(wc, wn, pow(wx.z, 0.5));
	   float DropFrac = frac(rp3.w+CurrentTime);
	   float TimeFrac = DropFrac - 1.0 + rp3.x;
	   float DropFactor = saturate(0.20 + Weight * 0.80 - DropFrac);
	   float PI = 3.14159265359;
	   float FF = DropFactor*rp3.x * sin(clamp(TimeFrac * 12.0, 0.0, 5.0) * PI);
	   
	   float3 cpl = lerp(0.0, float3(rp3.yz*FF*rp0, 1.0), a1);
	   return cpl;
}

float3 CR(float3 n, float3 uv)
{
float Tr = Timer.x * 14000.0;
float4 TimeMul = float4(0.6, 0.7, 0.8, 1.0); 
float4 TimeAdd = float4(0.22, 0.44, 0.66, 0.9);
float4 Times = (Tr * TimeMul + TimeAdd) * 2.4;
       Times = frac(Times);
float4 Weights = 1.0 - float4(0, 0.25, 0.5, 0.75);
       Weights = saturate(Weights * 4.0);
float3 Ripple1 = ComputeRipple(uv * 0.18 + float2(0.25, 0.0),  Times.x, Weights.x);
float3 Ripple2 = ComputeRipple(uv * 0.25 + float2(-0.55, 0.3), Times.y, Weights.y);
float3 Ripple3 = ComputeRipple(uv * 0.28 + float2(0.6, 0.85),  Times.z, Weights.z);
float3 Ripple4 = ComputeRipple(uv * 0.30 + float2(0.5, -0.75),  Times.w, Weights.w);
//float3 Ripple5 = ComputeRipple(uv * 0.32 + float2(0.3, -0.42),  Times.w, Weights.w);
//float3 Ripple6 = ComputeRipple(uv * 0.20 + float2(0.35, 0.1),  Times.x, Weights.x);
	n=normalize(float3(n.xy+Ripple1.xy,n.z));
	n=normalize(float3(n.xy+Ripple2.xy,n.z));
	n=normalize(float3(n.xy+Ripple3.xy,n.z));
	n=normalize(float3(n.xy+Ripple4.xy,n.z));	
//	n=normalize(float3(n.xy+Ripple5.xy,n.z));
//	n=normalize(float3(n.xy+Ripple6.xy,n.z));	
	return n;	
}

float3 ComputeRelief(float2 uv, float CurrentTime, float3 Weight)
{
    float tx = (Timer.x * 1000.0) * 1.05;
    float4 rf = tex2Dlod(SamplerWrl,  float4(uv.xy + float2(0.5, 0.5*tx), 0.0, 0.0));
		   rf.yz = rf.yz * 2.0-1.0;		   
       float RelfSt = step(0.2, rf.a);		   
	   float DropFrac = frac(rf.w+CurrentTime);
	   float3 DropFactor = saturate(0.2 + Weight * 0.8 - DropFrac);	
	   float3 FF = DropFactor * rf.x*RelfSt;
    return float3(rf.yz*FF, 1.0);
}

float3 PR(float3 n, float3 uv)
{
float4 t = 1.99;
float3 w = float3(1.2, 1.2, 1.2);
float tx = (Timer.x * 1000.0) * 1.05;	   
float3 Relief1 = ComputeRelief(uv * 0.05 + float2(0.55, 0.89),  t.x, w);
float3 Relief2 = ComputeRelief(uv * 0.25 + float2(0.22, 0.33)*tx,  t.x, w);	
	n=normalize(float3(n.xy+Relief1.xy,n.z)); 
	n=normalize(float3(n.xy+Relief2.xy,n.z));
    return n;
}

bool Puddles_alltime
<
         string UIName="Puddles_alltime";	
> = {true};

float Puddles_alltime_coverage
<
        string UIName="Puddles_alltime - coverage";
        string UIWidget="Spinner";
        float UIMin=0.0;
        float UIMax=2.0;
> = {0.65};

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

float Lightning(in float2 coord)
{	
    float2 uv = coord.xy + 0.1 * 1.5;
	float2 p = (uv + 0.5) / (-3.25*2.0);
	float fl0 = tex2Dlod(SamplerLightning, float4(p, 0.0, 0.0));
	float fl1 = tex2Dlod(SamplerLightning2, float4(p, 0.0, 0.0));	
	float ti = max(sin(Timer.x*100000.0*2.0),0);			
	float mix = lerp(fl0, fl1, ti);	
	return mix*1.5;
}

float Lightning2(in float2 coord, float v1, float v2)
{	
    float2 uv = coord.xy + 0.1 * 1.5;
	float2 p = (float2(v1, v2)+ uv + 0.5) / (-3.25*2.0);
	float fl0 = tex2Dlod(SamplerLightning, float4(p, 0.0, 0.0));
	float fl1 = tex2Dlod(SamplerLightning2, float4(p, 0.0, 0.0));	
	float ti = max(sin(Timer.x*100000.0*2.0),0);			
	float mix = lerp(fl0, fl1, ti);	
	return mix*1.5;
}

float4 PS_DX2(VS_OUTPUT_POST IN, float2 vPos : VPOS) : COLOR
{
   float4 r0 = tex2D(SamplerColor, IN.txcoord.xy);	
   float4 r1 = r0;	   
   float2 cd = IN.txcoord.xy;
   float d0 = tex2D(SamplerDepth, IN.txcoord.xy).x;
   float d1 = linearlizeDepth(tex2D(SamplerDepth, IN.txcoord.xy).x)*1.0;	   
   float4 tv0;
   float4 wpos;
   tv0.xy = IN.txcoord.xy*2.0-1.0;
   tv0.y =-tv0.y;
   tv0.z = d0;
   tv0.w = 1.0;
   wpos.x = dot(tv0, MatrixInverseVPRotation[0]);
   wpos.y = dot(tv0, MatrixInverseVPRotation[1]);
   wpos.z = dot(tv0, MatrixInverseVPRotation[2]);
   wpos.w = dot(tv0, MatrixInverseVPRotation[3]);
   wpos.xyz/= wpos.w;
   float4 npos;
   float4 n0 = tex2D(SamplerNormal, IN.txcoord.xy)*2.0-1.0;
   float4 tv1;
   tv1.xy = n0.xy;
   tv1.x = -tv1.x;
   tv1.z = n0.z;
   npos.x = dot(tv1.xyz, MatrixInverseView[0]);
   npos.y = dot(tv1.xyz, MatrixInverseView[1]);
   npos.z = dot(tv1.xyz, MatrixInverseView[2]);
   float3 n2 = normalize(npos.xyz);	  
          n2 = CR(n2.xyz, wpos+CameraPosition);  
          n2 = PR(n2.xyz, wpos+CameraPosition);  		  
    float3 vp = float3(0.0, 0.0, 1.03);
    float4 wc = wpd0(cd, d0);
	float rx = tex2Dlod(SamplerPuddle, float4(wc.xy*0.012, 0.0, 0.0));	
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
float rc;
float rn;
float4 Wx = WeatherAndTime;	
if (Wx.x==0,1) rc = 2.0;
if (Wx.y==0,1) rn = 2.0;
if (Wx.x==8)   rc = 0.2;
if (Wx.y==8)   rn = 0.2;
if (Wx.x==16)  rc = 0.2;
if (Wx.y==16)  rn = 0.2;	  
float wmix1 = lerp(rc, rn, pow(Wx.z, 8.5));
	float coverage = lerp(wmix1, Puddles_alltime_coverage, Puddles_alltime);
		  coverage = lerp(coverage, 0.87, 2.0);		  		  
	float density = 0.1 - 2.0 * 0.3;		  
	      rx = Coverage(coverage, density, rx);	
	float rh = smoothstep(0.999, 1.0, dot(-vp, -fs0(cd)));
	float4 rx2 = rx; 
	       rx2.a = min(1.0, rx2.a);  
	  float4 skyenv=tex2D(SamplerEnv, cd);		   
    float3 ssr = skyenv;	
	float4 rfl = reflection(n2.xyz, IN.txcoord.xy);		
           rfl = (rfl * 1.05); 
    float3 wp0 = wpd1(cd.xy);
    float3 p3 = float3(0.353, 0.353, 1.0);
           n2 = (normalize(n2.xyz*p3)*0.98);
    float3 nz = reflect(wp0.xyz, n2.xyz); 
    float3 rw = ((1000.0/0.01)*nz)/1000.0;	
    float3 ref11 = (wp0+rw);
    float2 rd = wpd2(ref11.xyz);		
	float  nf0 = saturate(8.0*(rd.y));
	       nf0 = pow(nf0, 1.5);
           ssr=lerp(ssr, rfl.xyz, rfl.w*nf0);
	float4 n1 = tex2D(SamplerNormal, IN.txcoord.xy);     
	float  rmix;		 
	       rmix=saturate(n1.z)*saturate(n1.z);
	       rmix*=rmix;					
		ssr.xyz=lerp(ssr, skyenv, rmix);
		ssr.xyz=lerp(ssr, r0*0.15, pow(rmix, 1.2));	
float wcr;
float wnr;
if (Wx.x==0,1) wcr = 0.0;
if (Wx.y==0,1) wnr = 0.0;
if (Wx.x==8)   wcr = 1.0;
if (Wx.x==16)  wcr = 1.0;
if (Wx.y==8)   wnr = 1.0;
if (Wx.y==16)  wnr = 1.0;	  
	float wmix = lerp(wcr, wnr, pow(Wx.z, 0.07));
	float4 nw = tex2D(SamplerNormal, IN.txcoord.xy).w;	  
	if (nw.w<1.0) nw = 0.0;
		  
	float3 n3 = CR(float3(0.0, 1.0, 1.0), wpos+CameraPosition);  		 
	float3 sv = normalize(float3(0.0, 0.0, 1.0));	 	
    float t1 = lerp(0.04, 0.04, x1);
          t1 = lerp(t1, 0.04, x2);
          t1 = lerp(t1, 0.09, x3);
          t1 = lerp(t1, 0.09, x4);
          t1 = lerp(t1, 0.09, xE);
          t1 = lerp(t1, 0.09, x5);
          t1 = lerp(t1, 0.09, x6);		 
          t1 = lerp(t1, 0.08, x7);
		  t1 = lerp(t1, 0.07, xG);
		  t1 = lerp(t1, 0.06, xZ);
          t1 = lerp(t1, 0.05, x8);  
          t1 = lerp(t1, 0.04, x9);	
          t1*= 1.1;		  
   float factor = t1 - dot(sv, n2);
         factor = pow(factor, 12.0);
   float fr = (factor*factor); 
         fr/= 4.5;	
	float n6 = saturate(n3);
	float lp = lerp(0.0, fr, n6);
	
	if(d1 < 3000)r0.xyz = lerp(r0, 0.75*ssr+lp, nw*rh*rx2*1.0);  
	float4 rp = float4(normalize(wpos.xyz), 1.0);
	float3 rpos = normalize(rp.xyz)*1.0;
	float3x3 vec = float3x3(1,0,0,0,cos(0.25),sin(0.25),0,-sin(0.25),cos(0.25));				   
		   rpos = mul(rpos, vec);	
	float2 f0=float2(acos(rpos.z),atan(rpos.x/rpos.y));		  		  
	float2 RainTimer = float2(Timer.x*31000.0, 0.0); 
	float4 raindx2 = Rain2(-f0.xy+RainTimer*0.7);
	float4 raindx1 = Rain1(-f0.xy+RainTimer*1.0);	  
    float4 raindx0 = Rain0(-f0.xy+RainTimer*2.9);				 
	float a0 = smoothstep(1.0, 0.90, abs(rpos.z));
    float pow2 = length(wpos.xyz);
    float distance = 0.0;	   
	   float3 rd0 = 1.0-1.0/exp(length(distance-pow2)/0.1);
	   float3 rd1 = 1.0-1.0/exp(length(distance-pow2)/8.0);
	   float3 rd2 = 1.0-1.0/exp(length(distance-pow2)/95.0);		   
	   float3 dr0 = step(0.15, rd0*raindx0);
	   float3 dr1 = step(0.1, rd1*raindx1);
	   float3 dr2 = step(0.1, rd2*raindx2);
	   //float4   rn0 = tex2Dlod(SamplerColor, float4(IN.txcoord.xy+raindx2*float2(0.0, -0.10), 0.0, 4.0));	
	   //		  	rn0*= 0.63; 	   
	   float4   rn1 = tex2Dlod(SamplerColor, float4(IN.txcoord.xy+raindx1*float2(0.0, -0.13), 0.0, 4.0));  
	   		  	rn1*= 0.63; 	   
	   float4   rn2 = tex2Dlod(SamplerColor, float4(IN.txcoord.xy+raindx0*float2(0.0, -0.15), 0.0, 4.0));
	   	   		rn2*= 0.63; 
	   		  	//rn0.a = min(0.9, rn0.a); 
		  	    rn1.a = min(0.9, rn1.a);	   
		  	    rn2.a = min(0.9, rn2.a);
    float t2 = lerp(0.0, 0.0, x1);
          t2 = lerp(t2, 0.10, x2);
          t2 = lerp(t2, 0.20, x3);
          t2 = lerp(t2, 0.25, x4);
          t2 = lerp(t2, 0.25, xE);
          t2 = lerp(t2, 0.25, x5);
          t2 = lerp(t2, 0.25, x6);		 
          t2 = lerp(t2, 0.15, x7);
		  t2 = lerp(t2, 0.12, xG);
		  t2 = lerp(t2, 0.10, xZ);
          t2 = lerp(t2, 0.10, x8);  
          t2 = lerp(t2, 0.07, x9);						
	     //r0.xyz = lerp(r0, 1.05*rn0+t2, wmix*dr2*a0*rn0.a);
		 r0.xyz = lerp(r0, 1.07*rn1+t2, wmix*dr1*a0*rn1.a);
		 r0.xyz = lerp(r0, 1.09*rn2+t2, wmix*dr0*a0*rn2.a);
float wc1;
float wn1;
if (Wx.x==0,1) wc1 = 0.0;
if (Wx.y==0,1) wn1 = 0.0;
if (Wx.x==8)   wc1 = 1.0;
if (Wx.y==8)   wn1 = 1.0;
if (Wx.x==16)  wc1 = 1.0;
if (Wx.y==16)  wn1 = 1.0;
float wmix2 = lerp(wc1, wn1, Wx.z);
float3 fx0 = FogFarColor;
	float3 st3 = normalize(fx0.xyz);
	float3 cs3=fx0.xyz/st3.xyz;
	       cs3=pow(cs3, 6.0);
	   fx0.xyz = cs3;	
float fx = saturate(fx0);
float3 npos0 = normalize(wpos.xyz);

   float vecx2 = 5.45;
   float vecx3 = 3.85; 
   float vecx4 = 1.65;	
   float3 ns0 = normalize(wpos.xzy);
   float3x3 vec0 = float3x3(1, 0, 0, 
                            0, cos(vecx2), sin(vecx2),
							0, -sin(vecx2), cos(vecx2));   
		   ns0 = mul(ns0, vec0);		   
	float3 np1 = normalize(ns0);
	
   float3 ns1 = normalize(wpos.xyz);
   float3 ns2 = normalize(wpos.xzy);   
   float3x3 vec1 = float3x3(1, 0, 0, 
                            0, cos(vecx3), sin(vecx3),
							0, -sin(vecx3), cos(vecx3));   
							
   float3x3 vec2 = float3x3(0, 0, 1, 
                            cos(vecx4), sin(vecx4), 0,
							-sin(vecx4), cos(vecx4), 0);   							
		   ns1 = mul(ns1, vec1);
		   ns2 = mul(ns2, vec2);
		   	   
	float3 np2 = normalize(ns1.xyz);	
	float3 np3 = normalize(ns2.xyz);
float4 fx2 = Lightning(np1*10.0)*3.0;
	   fx2*= smoothstep(0.0, 0.3, npos0.z);
float4 fx3 = Lightning(np2*10.0)*3.0;   
	   fx3*= smoothstep(0.0, 0.3, npos0.z);
float4 fx4 = Lightning2(np3*10.0, -7.9, -7.3)*3.0;   	
float4 fx5 = Lightning2(np3*10.0, 1.0, -6.4)*3.0;   	   
 
	  	r0.xyz*= lerp(1.0, 1.8, wmix2*fx);
float L = wpos.xyz;	  

	  float ti = max(sin(Timer.x*100000.0*0.5),0);
	  float ti2 = max(sin(Timer.x*100000.0*1.0),0);	  
	float3 mix0 = lerp(fx2, fx3, ti);  
	float3 mix1 = lerp(fx4, fx5, ti);  	
	
	float3 mix2 = lerp(mix0, mix1, ti2); 

	   if( 6000 < L || d0 == 1)
     {
		 r0.xyz =lerp(r0, r0+mix2, wmix2*fx);
      }	
	  
  return r0;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

technique PostProcess
{
   pass P0
   {
      VertexShader = compile vs_3_0 VS_PostProcess();
      PixelShader  = compile ps_3_0 PS_DX0();
   }
}

technique PostProcess2
{
	pass P0
	{
      VertexShader = compile vs_3_0 VS_PostProcess();
      PixelShader  = compile ps_3_0 PS_DX1();
	}
}

technique PostProcess3
{
   pass P0
   {
      VertexShader = compile vs_3_0 VS_PostProcess();
      PixelShader  = compile ps_3_0 PS_DX2();
   }
}
