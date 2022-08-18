
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

float4 tempF1; float4 tempF2; float4 tempF3; float4 ScreenSize; float ENightDayFactor; float EInteriorFactor;
float4 WeatherAndTime; float4 Timer; float FieldOfView; float GameTime; float4 SunDirection; 
float4 CustomShaderConstants1[8]; float4 MatrixVP[4]; float4 MatrixInverseVP[4]; float4 MatrixVPRotation[4];
float4 MatrixInverseVPRotation[4]; float4 MatrixView[4]; float4 MatrixInverseView[4]; float4 CameraPosition;
float4x4 MatrixWVP; float4x4 MatrixWVPInverse; float4x4 MatrixWorld; float4x4 MatrixProj; float4 diffColor;
float4 specColor; float4 ambColor; float4 FogParam; float4 FogFarColor; float4 lightDiffuse[8]; float4 lightSpecular[8];
float4 lightDirection[8]; float4 VehicleParameters1;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//Textures
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

texture2D texOriginal;
texture2D texRefl;
texture2D texEnv;
texture2D texNoise < string ResourceName = "NoiseCar.png"; >;
sampler2D SamplerOriginal = sampler_state { Texture   = <texOriginal>; };
textureCUBE texColeso < string ResourceName = "SkyCar.dds"; string ResourceType = "CUBE"; >;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//Sampler Inputs
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

sampler2D SamplerRefl = sampler_state
{
	Texture   = <texRefl>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU  = Clamp;
	AddressV  = Clamp;
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

samplerCUBE SamplerCube = sampler_state
{
	Texture   = <texColeso>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU  = Wrap;
	AddressV  = Wrap;
	SRGBTexture=TRUE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D NoiseSampler = sampler_state 
{
	Texture = <texNoise>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Wrap;
	AddressV = Wrap;
	AddressW = Wrap;
	SRGBTexture=TRUE;
	MaxMipLevel=2;
	MipMapLodBias=2;
};

struct PS_OUTPUT3
{
	float4 Color[3] : COLOR0;
};

struct VS_INPUT_N
{
	float3	pos : POSITION;
	float3	normal : NORMAL;
	float2	txcoord0 : TEXCOORD0;
};

struct VS_OUTPUT
{
	float4	pos : POSITION;
	float2	txcoord0 : TEXCOORD0;
	float3	viewnormal : TEXCOORD3;
	float3	eyedir : TEXCOORD4;
	float3	wnormal : TEXCOORD5;
	float4	vposition : TEXCOORD6;
	float3	normal : TEXCOORD7;
};

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
VS_OUTPUT VS_Draw(VS_INPUT_N IN)
{
    VS_OUTPUT OUT;
	float4	pos;
	pos.xyz=IN.pos.xyz;
	pos.w=1.0;
	float4	tpos;
	tpos=mul(pos, MatrixWVP);
	OUT.pos=tpos;
	OUT.vposition=tpos;
	OUT.txcoord0=IN.txcoord0;	
	float3	wnormal=normalize(mul(IN.normal.xyz, MatrixWorld));
	float3	normal;
	normal.x=dot(wnormal.xyz, MatrixView[0]);
	normal.y=dot(wnormal.xyz, MatrixView[1]);
	normal.z=dot(wnormal.xyz, MatrixView[2]);
	OUT.viewnormal=normalize(normal.xyz);
	OUT.normal=normalize(mul(IN.normal.xyz, MatrixWVP));
	OUT.wnormal=wnormal;
	float3	campos=CameraPosition;
	OUT.eyedir=(mul(pos, MatrixWorld) - campos);
    return OUT;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

float ReflContrast
<
	string UIName="Reflection - Contrast";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=4.0;
> = {2.0};

float ReflContrast2
<
	string UIName="Reflection - Contrast(Night)";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=4.0;
> = {1.5};

float ReflSaturate
<
	string UIName="Reflection - Saturate";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=4.0;
> = {0.65};

float refd
<
	string UIName="Reflection - Day";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=5.0;
> = {1.0};

float refn
<
	string UIName="Reflection - Night";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=5.0;
> = {1.0};

float refw
<
	string UIName="Reflection - Weather";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=5.0;
> = {1.0};

float skyref
<
        string UIName="Sky - Brightness";
        string UIWidget="Spinner";
        float UIMin=0.0;
        float UIMax=3.0;
> = {1.4};

float refglass
<
	string UIName="Glass - Reflection";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=4.0;
> = {1.0};

float carglass
<
	string UIName="Glass - Tinting";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=1.0;
> = {1.0};

bool EReflEnable
<
	string UIName="SreenSpaceRefl: Enable";
> = {true};

bool cre
<
	string UIName="Chrome: Enable";
> = {true};

float carchrome
<
	string UIName="Chrome - Brightness";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=5.0;
> = {1.0};

float EFresnelFactor
<
	string UIName="Vehicle: FresnelFactor";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=1.0;
> = {0.84};

float EF1
<
	string UIName="Fresnel: coverage1";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=3.0;
> = {0.4};

float EF2
<
	string UIName="Fresnel: coverage2";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=3.0;
> = {0.55};

float CarContrast
<
	string UIName="Vehicle: Contrast";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10000.0;
> = {1.0};

float CarSaturate
<
	string UIName="Vehicle: Saturate";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10000.0;
> = {1.0};

float wLhights
<
	string UIName="wheels and interior";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=5.0;
> = {1.5};

float ftColor3
<
	string UIName="Car - Brightness";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=2.0;
> = {0.14};

float flighting
<
	string UIName="Lighting - Brightness";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10.0;
> = {1.0};

float LightingContrast
<
	string UIName="Lighting - Contrast";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10000.0;
> = {1.0};

float LightingSaturate
<
	string UIName="Lighting - Saturate";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10000.0;
> = {1.0};


float L1
<
	string UIName="Lighting - L1";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10.0;
> = {1.0};


float L2
<
	string UIName="Lighting - L2";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=20.0;
> = {1.0};

float L3
<
	string UIName="Lighting - L2(weather)";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=20.0;
> = {1.0};

float RefSize
<
	string UIName="RefSize";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10.0;
> = {1.0};


////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

float4 wpd1(float2 cd, float d)
{
float4 tv; 
       tv.xy = cd.xy*2.0-1.0;   
       tv.y = -tv.y;   
       tv.z = d;
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

float4 reflection(float3 n, float2 cd, float d)
{
  float3 wpos = wpd1(cd.xy, d);
  float3 v = {0.35, 0.35, 1.15};
         n = normalize(n.xyz*v);		 
  float3 n0 = reflect(wpos.xyz, n.xyz);		 
  float3 n1 = ((1000.0/0.01)*n0)/1000.0;		
  float3 r0 = (wpos+n1);
  float2 r1 = wpd2(r0.xyz);	
  float4 r2 = tex2Dlod(SamplerRefl, float4(r1.xy, 0.0, RefSize));
	     r2.xyz+=0.000001;
  float3 st0 = normalize(r2.xyz);
  float3 ct0 = r2.xyz/st0.xyz;
	     ct0=pow(ct0, 0.70);	
	     st0.xyz = pow(st0.xyz, 0.95);
	     r2.xyz = ct0*st0.xyz;	
	     r2.w = r1.y<0.0||r1.y>1.0 ? 0.0:1.0;
	     r2.w*= r1.x<0.0||r1.x>1.0 ? 0.0:1.0;			 
  return r2;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

PS_OUTPUT3 PS_Draw(VS_OUTPUT IN, in float2 vpos : VPOS)
{
	float4	cd;
    float2  cd0 = IN.txcoord0.xy;
	float4  wx = WeatherAndTime;	
	float3	n0 = normalize(IN.normal.xyz);
	float3	wn = normalize(IN.wnormal.xyz);
	float3	ed = normalize(-IN.vposition.xyz);
	float3	wed = normalize(-IN.eyedir.xyz);
	float4	tex = tex2D(SamplerOriginal, cd0);
	float4	r0;
	float4	rfl;
	float2	rfl0;
	
	float4 fp1;
           fp1.xyz = diffColor*specColor*specColor.w*0.085;
           fp1.w = min(fp1.x, min(fp1.y, fp1.z));
           fp1.w = saturate(1.0-fp1.w);
           fp1.w*= fp1.w;
	cd.zw = 0.0;		
	cd.w = fp1.w*4.0;
		   
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

   float3 rv1 = {-0.2, -1.1, 1.0};

   float opacity = saturate(tex.a*diffColor.a);
   float3 refl = reflect(ed, n0);	
	      rfl0.xy = (IN.vposition.xy /IN.vposition.w)*float2(0.5, -0.5) + 0.5;
	      rfl0.xy+= (refl.xy*float2(-1.0, 1.0)*0.5);	
	cd.xy = rfl0;
	float smix;
		  smix = saturate(wn.z)*saturate(-n0.z);
		  smix*= smix;	
//-------------------------------------
//-------------------------------------						
	cd.y = 1.0-((wn.z*0.5)+0.5);
	cd.x = 0.5;	
	
	float3 wp = reflect(wed, wn);
	float4 tcube = texCUBE(SamplerCube, -wp.xzy);   	  		   
	float2 txcoord1;
	       txcoord1.xy = (IN.vposition.xy /IN.vposition.w)*float2(0.5, -0.5) + 0.5;
	       txcoord1.xy+= 0.5*float2(ScreenSize.y, ScreenSize.y*ScreenSize.z);
				
    float4 b0 = 1.0;	
	float3 nont01 = b0.xyz*0.0;
	nont01.xyz = min(nont01, b0);
	float3 tr2 = b0.xyz*10.0/saturate(opacity + 0.02 + VehicleParameters1.y);
	       b0.xyz = lerp(tr2, nont01, saturate(opacity+VehicleParameters1.y));	
    float mask8 = max(b0.x, max(b0.y, b0.z));	
		  mask8 = saturate(mask8*0.3)*1.0;		

    float4 noise = tex2Dlod(NoiseSampler, float4(txcoord1.xy*4.0, 0.0, 0.0));	
	       noise = lerp(noise, 0.0, mask8);  		   
	float d = (IN.vposition.z/IN.vposition.w);			   
	rfl = reflection(normalize(wn), txcoord1+(0.010*noise), d); 		   
    float4 sky = tex2D(SamplerEnv, rv1+cd.xy+(0.005*noise));	  
    float3 tc0 = lerp(0.01, 0.1, x1); 
           tc0 = lerp(tc0, 0.3, x2); 
           tc0 = lerp(tc0, 1.0, x3); 
           tc0 = lerp(tc0, 1.0, x4); 
           tc0 = lerp(tc0, 1.0, xE); 
           tc0 = lerp(tc0, 1.0, x5); 
           tc0 = lerp(tc0, 1.0, x6); 		 
           tc0 = lerp(tc0, 0.8, x7); 
		   tc0 = lerp(tc0, 0.6, xG); 
		   tc0 = lerp(tc0, 0.3, xZ); 
           tc0 = lerp(tc0, 0.2, x8);	
           tc0 = lerp(tc0, 0.01, x9);	  
	  
	float4 sc = tcube*sky;
	float sc0 = saturate(tcube)*tc0;	  
	  
float4 cubeC;
float4 cubeN;	  
if (wx.x==0,1) cubeC = sc;
if (wx.y==0,1) cubeN = sc;
if (wx.x==4) cubeC = sc0;
if (wx.x==7) cubeC = sc0;
if (wx.x==8) cubeC = sc0;
if (wx.x==9) cubeC = sc0;
if (wx.x==12) cubeC = sc0;
if (wx.x==15) cubeC = sc0;
if (wx.x==16) cubeC = sc0;
if (wx.y==4) cubeN = sc0;
if (wx.y==7) cubeN = sc0;
if (wx.y==8) cubeN = sc0;
if (wx.y==9) cubeN = sc0;
if (wx.y==12) cubeN = sc0;
if (wx.y==15) cubeN = sc0;
if (wx.y==16) cubeN = sc0;	  
    float4 wmix0 = lerp(cubeC, cubeN, wx.z);		  
    float4 skymix = wmix0;

	float3 st0 = normalize(rfl.xyz);
	float3 ct0=rfl.xyz/st0.xyz;
	       ct0=pow(ct0, 1.10);
	       st0.xyz = pow(st0.xyz, 1.15);	  
	rfl.xyz = ct0*st0.xyz;	

	float3 st1 = normalize(skymix.xyz);
	float3 ct1=skymix.xyz/st1.xyz;
	       ct1=pow(ct1, 0.5);
	       st1.xyz = pow(st1.xyz, 0.5);	  
	skymix.xyz = ct1*st1.xyz;
	skymix.xyz*= skyref;
	
	float3 wpos = wpd1(txcoord1.xy, d);
    float3 v = {0.35, 0.35, 1.15};
    float3 n1 = normalize(wn*v);
    float3 n2 = reflect(wpos.xyz, n1.xyz);
    float3 rw = ((1000.0/0.01)*n2)/1000.0;
    float3 ref = (wpos+rw);
    float2 rd = wpd2(ref.xyz);	

	float nf02 = saturate(8.0*(rd.y));
          nf02*= saturate(8.0*(rd.x));					
	      nf02 = pow(nf02, 1.0);			
	float nf03 = saturate(8.0+rd.y*(-8.0));
          nf03*= saturate(8.0+rd.x*(-8.0));					
	      nf03 = pow(nf03, 1.0);

    float4 rfl1 = skymix;
		   rfl.xyz = lerp(rfl1, rfl,  pow(0.01*0.5, smix));
		   
	float nf = saturate(1.0+smix*0.3 -abs(n0.z));
	float ff = pow(nf, 2.5);

		  rfl.xyz = lerp(rfl, rfl1,  pow(0.01*0.01, ff));   
		  if (EReflEnable==true) rfl1 = lerp(rfl1, rfl, rfl.w*nf02*nf03);
		  	  
	float3 st3 = normalize(rfl1.xyz);
	float3 ct3 = rfl1.xyz/st3.xyz;	
	 float tc1 = lerp(ReflContrast2, ReflContrast2, x1); 
           tc1 = lerp(tc1, ReflContrast, x2); 
           tc1 = lerp(tc1, ReflContrast, x3); 
           tc1 = lerp(tc1, ReflContrast, x4); 
           tc1 = lerp(tc1, ReflContrast, xE); 
           tc1 = lerp(tc1, ReflContrast, x5); 
           tc1 = lerp(tc1, ReflContrast, x6); 		 
           tc1 = lerp(tc1, ReflContrast, x7); 
		   tc1 = lerp(tc1, ReflContrast, xG); 
		   tc1 = lerp(tc1, ReflContrast, xZ); 
           tc1 = lerp(tc1, ReflContrast, x8);	
           tc1 = lerp(tc1, ReflContrast2, x9);	  	
	       ct3 = pow(ct3, tc1);
	       st3.xyz = pow(st3.xyz, ReflSaturate);		   
	       rfl1.xyz = ct3*st3.xyz;	
//-------------------------------------
//-------------------------------------	
    float3 sv0 = SunDirection.xyz;
    float3 sv1 = normalize(float3(0.0, 0.0, 0.0));

    float yy1 = smoothstep(0.0, 1.0, t0);	
    float yy2 = smoothstep(1.0, 23.0, t0);
    float yy3 = smoothstep(23.0, 24.0, t0);	
   
    float3 sv = lerp(sv1, sv0, yy1);
           sv = lerp(sv, sv0, yy2);
           sv = lerp(sv, sv1, yy3);  

    float3 hv = normalize(sv.xyz+wed.xyz);	
    float gf = saturate((0.01*0.001) + dot(hv, wn.xyz));
		  gf = pow(gf, 1500.0 * 5.0);
    float gm = (gf*gf);	
    float gl = (gm/0.01)* 0.2;	
//-------------------------------------
//-------------------------------------	
	float distfade = saturate(IN.vposition.z * 0.014);
	      rfl1.xyz = lerp(rfl1, skymix, smix);
	      rfl1.xyz*= 1.0-distfade*distfade;

    float4 dc = diffColor;
    float4 tc = tex;	
    float4 tc2 = tex;		
    float4 gl0 = 1.0;
    float4 gl1 = 1.0;	
	float3 st4 = normalize(dc.xyz);
	float3 ct4=dc.xyz/st4.xyz;
	       ct4=pow(ct4, CarContrast);
		   st4.xyz = pow(st4.xyz, CarSaturate);
	       dc.xyz = ct4*st4.xyz;
		   
   float Headlights = 0.25;
   float3 lf = lerp(Headlights*7.0, Headlights*7.0, x1);
          lf = lerp(lf, Headlights*6.0, x2);
          lf = lerp(lf, Headlights, x3);
          lf = lerp(lf, Headlights, x4);
          lf = lerp(lf, Headlights, xE);
          lf = lerp(lf, Headlights, x5);
          lf = lerp(lf, Headlights, x6);	 
          lf = lerp(lf, Headlights, x7);
		  lf = lerp(lf, Headlights, xG);
		  lf = lerp(lf, Headlights, xZ);
          lf = lerp(lf, Headlights*5.0, x8);		
          lf = lerp(lf, Headlights*6.5, x9);
		
	float3 lighting = 1.53*(dc+(ambColor*lf))*tex;
		
	float3 st5 = normalize(tc.xyz);
	float3 ct5 = tc.xyz/st5.xyz;
	       ct5 = pow(ct5, 9.0);
		   st5.xyz = pow(st5.xyz, -5.9);
	       tc.xyz = ct5*st5.xyz;
		
	float3 st6 = normalize(tc2.xyz);
	float3 ct6 = tc2.xyz/st6.xyz;
	       ct6 = pow(ct6, 8.0);
		   st6.xyz = pow(st6.xyz, 4.0);
	       tc2.xyz = ct6*st6.xyz;
		
	float3 nt = gl0.xyz*0.6;
	       nt.xyz = min(nt, gl0);
	float3 tr = gl0.xyz*(-3.0)/saturate(opacity + 0.02 + VehicleParameters1.y);
	       gl0.xyz = lerp(tr, nt, saturate(opacity+VehicleParameters1.y));
	float3 nt0 = gl1.xyz*0.0;
	       nt0.xyz = min(nt0, gl1);
	float3 tr0 = gl1.xyz*10.0/saturate(opacity + 0.02 + VehicleParameters1.y);
	       gl1.xyz = lerp(tr0, nt0, saturate(opacity+VehicleParameters1.y));

    float4 sc1 = specColor*1000.0;
    float4 sc2 = tex;
	float3 g2 = normalize(sc1.xyz);
	float3 g4 = normalize(sc2.xyz);
	float3 s2 = sc1.xyz/g2.xyz;
	float3 s4 = sc2.xyz/g4.xyz;
	       s2 = pow(s2, -15.0);
	       s4 = pow(s4, -34.0);
		   g4.xyz = pow(g4.xyz, 3.0);
	       sc2.xyz = s4*g4.xyz;
		   
	float mask0 = max(specColor.x, max(specColor.y, specColor.z));
	      mask0 = saturate(mask0*1000.0);				
	float mask1 = max(s2.x, max(s2.y, s2.z));
	      mask1 = saturate(mask1);
	float mask2 = max(sc2.x, max(sc2.y, sc2.z));
	      mask2 = saturate(mask2);					
	float mask3 = max(tc.x, max(tc.y, tc.z));
	      mask3 = saturate(mask3);				
	float mask4 = max(tc2.x, max(tc2.y, tc2.z));
	      mask4 = saturate(mask4);						
	float mask5 = max(gl0.x, max(gl0.y, gl0.z));
	      mask5 = saturate(mask5);	
    float mask6 = max(gl1.x, max(gl1.y, gl1.z));	
		  mask6 = saturate(mask6*0.3)*1.0;				

    float3 sv3 = normalize(float3(0.0, 0.0, 1.0)+wed.xyz);   
    float3 sv2 = normalize(sv.xyz+wed.xyz);
    float3 sv4 = normalize(sv.xyz);
    float3 ScCurrent;
    float3 ScNext;
if (wx.x==0,1) ScCurrent = sv2;
if (wx.y==0,1) ScNext = sv2;
if (wx.x==4) ScCurrent = sv3;
if (wx.x==7) ScCurrent = sv3;
if (wx.x==8) ScCurrent = sv3;
if (wx.x==9) ScCurrent = sv3;
if (wx.x==12) ScCurrent = sv3;
if (wx.x==15) ScCurrent = sv3;
if (wx.x==16) ScCurrent = sv3;
if (wx.y==4) ScNext = sv3;
if (wx.y==7) ScNext = sv3;
if (wx.y==8) ScNext = sv3;
if (wx.y==9) ScNext = sv3;
if (wx.y==12) ScNext = sv3;
if (wx.y==15) ScNext = sv3;
if (wx.y==16) ScNext = sv3;	  

    float3 svmix = lerp(ScCurrent, ScNext, wx.z);	 
    float3 ScCurrent2;
    float3 ScNext2;
if (wx.x==0,1) ScCurrent2 = sv4;
if (wx.y==0,1) ScNext2 = sv4;
if (wx.x==4) ScCurrent2 = sv3;
if (wx.x==7) ScCurrent2 = sv3;
if (wx.x==8) ScCurrent2 = sv3;
if (wx.x==9) ScCurrent2 = sv3;
if (wx.x==12) ScCurrent2 = sv3;
if (wx.x==15) ScCurrent2 = sv3;
if (wx.x==16) ScCurrent2 = sv3;
if (wx.y==4) ScNext2 = sv3;
if (wx.y==7) ScNext2 = sv3;
if (wx.y==8) ScNext2 = sv3;
if (wx.y==9) ScNext2 = sv3;
if (wx.y==12) ScNext2 = sv3;
if (wx.y==15) ScNext2 = sv3;
if (wx.y==16) ScNext2 = sv3;	  

float re4 = L3;
float gc5;
float gn5;
if (wx.x==0,1) gc5 = L2;
if (wx.y==0,1) gn5 = L2;
if (wx.x==4) gc5 = re4;
if (wx.x==7) gc5 = re4;
if (wx.x==8) gc5 = re4;
if (wx.x==9) gc5 = re4;
if (wx.x==12) gc5 = re4;
if (wx.x==15) gc5 = re4;
if (wx.x==16) gc5 = re4;
if (wx.y==4) gn5 = re4;
if (wx.y==7) gn5 = re4;
if (wx.y==8) gn5 = re4;
if (wx.y==9) gn5 = re4;
if (wx.y==12) gn5 = re4;
if (wx.y==15) gn5 = re4;
if (wx.y==16) gn5 = re4;	

   float svmix5 = lerp(gc5, gn5, wx.z);	
   float3 svmix2 = lerp(ScCurrent2, ScNext2, wx.z);	
   float3 np = normalize(wn.xyz);
   float factor = 0.18*svmix5 - dot(-svmix, np);
         factor = pow(factor, 12.0*L1);
   float factor2 = 0.85*svmix5 - dot(-svmix2, np);
         factor2 = pow(factor2, 2.0*L1);
   float factor3 = 0.08 - dot(-svmix, np);
         factor3 = pow(factor3, 25.0);		 
   float factor0 = EF2 - dot(-wed, np);
         factor0 = pow(factor0, EF1);	
   float fr0 = (factor0*factor0); 
         fr0/= 2.5;		 		 
   float fr = (factor*factor); 
         fr/= 2.5;		 
   float fr2 = (factor2*factor2); 
         fr2/= 2.5;	
   float fr3 = (factor3*factor3); 
         fr3/= 2.5;	

float re2 = ftColor3*0.4;
float gc2;
float gn2;
if (wx.x==0,1) gc2 = ftColor3;
if (wx.y==0,1) gn2 = ftColor3;
if (wx.x==4) gc2 = re2;
if (wx.x==7) gc2 = re2;
if (wx.x==8) gc2 = re2;
if (wx.x==9) gc2 = re2;
if (wx.x==12) gc2 = re2;
if (wx.x==15) gc2 = re2;
if (wx.x==16) gc2 = re2;
if (wx.y==4) gn2 = re2;
if (wx.y==7) gn2 = re2;
if (wx.y==8) gn2 = re2;
if (wx.y==9) gn2 = re2;
if (wx.y==12) gn2 = re2;
if (wx.y==15) gn2 = re2;
if (wx.y==16) gn2 = re2;		
		   
    float wmix3 = lerp(gc2, gn2, wx.z);			 
	 
   float3 tl = lerp(0.06, 0.10, x1);
          tl = lerp(tl, 0.10, x2);
          tl = lerp(tl, wmix3*0.8, x3);
          tl = lerp(tl, wmix3*0.8, x4);
          tl = lerp(tl, wmix3*0.8, xE);
          tl = lerp(tl, wmix3*0.8, x5);
          tl = lerp(tl, wmix3*0.8, x6); 
          tl = lerp(tl, 0.20, x7);
		  tl = lerp(tl, 0.10, xG);
		  tl = lerp(tl, 0.08, xZ);
          tl = lerp(tl, 0.06, x8);	
          tl = lerp(tl, 0.06, x9);	   
   	
   r0.xyz = (saturate(tl*lighting)*2.0);	 
	 
   float3 tl0 = lerp(2.8, 2.0, x1);
          tl0 = lerp(tl0, 2.0, x2);
          tl0 = lerp(tl0, 1.2, x3);
          tl0 = lerp(tl0, 1.2, x4);
          tl0 = lerp(tl0, 1.2, xE);
          tl0 = lerp(tl0, 1.2, x5);
          tl0 = lerp(tl0, 1.2, x6);	 
          tl0 = lerp(tl0, 1.2, x7);
		  tl0 = lerp(tl0, 1.2, xG);
		  tl0 = lerp(tl0, 1.2, xZ);
          tl0 = lerp(tl0, 2.0, x8);		
          tl0 = lerp(tl0, 2.8, x9);	
		  
float3 ge0 = 0.0;
float3 gc6;
float3 gn6;
if (wx.x==0,1) gc6 = float3(0.235, 0.157, 0.0);
if (wx.y==0,1) gn6 = float3(0.235, 0.157, 0.0);
if (wx.x==4) gc6 = ge0;
if (wx.x==7) gc6 = ge0;
if (wx.x==8) gc6 = ge0;
if (wx.x==9) gc6 = ge0;
if (wx.x==12) gc6 = ge0;
if (wx.x==15) gc6 = ge0;
if (wx.x==16) gc6 = ge0;
if (wx.y==4) gn6 = ge0;
if (wx.y==7) gn6 = ge0;
if (wx.y==8) gn6 = ge0;
if (wx.y==9) gn6 = ge0;
if (wx.y==12) gn6 = ge0;
if (wx.y==15) gn6 = ge0;
if (wx.y==16) gn6 = ge0;		
		   
   float3 wmix6 = lerp(gc6, gn6, wx.z);			  
   float srs = saturate(dc);	 
   float3 cdl = lerp(float3(0.392, 0.392, 0.392), wmix6, pow(srs, 0.05));

float re3 = flighting*1.4;
float gc3;
float gn3;
if (wx.x==0,1) gc3 = flighting;
if (wx.y==0,1) gn3 = flighting;
if (wx.x==4) gc3 = re3;
if (wx.x==7) gc3 = re3;
if (wx.x==8) gc3 = re3;
if (wx.x==9) gc3 = re3;
if (wx.x==12) gc3 = re3;
if (wx.x==15) gc3 = re3;
if (wx.x==16) gc3 = re3;
if (wx.y==4) gn3 = re3;
if (wx.y==7) gn3 = re3;
if (wx.y==8) gn3 = re3;
if (wx.y==9) gn3 = re3;
if (wx.y==12) gn3 = re3;
if (wx.y==15) gn3 = re3;
if (wx.y==16) gn3 = re3;		
		   
    float wmix4 = lerp(gc3, gn3, wx.z);		   

   float ecd3 = 0.12*wmix4;   
   float3 cc = lerp(0.0, 0.00, x1);
          cc = lerp(cc, cdl*ecd3*0.2, x2);
          cc = lerp(cc, cdl*ecd3*0.6, x3);
          cc = lerp(cc, cdl*ecd3*0.9, x4);
          cc = lerp(cc, cdl*ecd3, xE);
          cc = lerp(cc, cdl*ecd3, x5);
          cc = lerp(cc, cdl*ecd3*0.5, x6);	 
          cc = lerp(cc, 0.00, x7);
		  cc = lerp(cc, 0.00, xG);
		  cc = lerp(cc, 0.00, xZ);
          cc = lerp(cc, 0.00, x8);
          cc*= 1.0;	
		  
	  
   
  float3 fday = 6.5*refd;
  float3 fnight = 6.0*refn;
  float3 sl = lerp(fnight, fday, yy1);
         sl = lerp(sl, fday, yy2);
         sl = lerp(sl, fnight, yy3); 

   float ecd2 = 0.32*wmix4;	 
   float3 cc2 = lerp(0.0, 0.00, x1);
          cc2 = lerp(cc2, float3(0.627, 0.627, 0.627)*ecd2*0.4, x2);
          cc2 = lerp(cc2, float3(0.627, 0.627, 0.627)*ecd2*0.4, x3);
          cc2 = lerp(cc2, float3(0.667, 0.667, 0.667)*ecd2*0.8, x4);
          cc2 = lerp(cc2, ecd2, xE);
          cc2 = lerp(cc2, ecd2, x5);
          cc2 = lerp(cc2, float3(0.667, 0.667, 0.667)*ecd2*0.5, x6);		 
          cc2 = lerp(cc2, 0.0, x7);
		  cc2 = lerp(cc2, 0.0, xG);
		  cc2 = lerp(cc2, 0.0, xZ);
          cc2 = lerp(cc2, 0.0, x8);
          cc2*= 1.0;

   float ecd = 0.05*wmix4;		  
   float3 cc1 = lerp(0.0, 0.00, x1);
          cc1 = lerp(cc1, float3(0.627, 0.627, 0.627)*ecd*0.4, x2);
          cc1 = lerp(cc1, float3(0.627, 0.627, 0.627)*ecd*5.4, x3);
          cc1 = lerp(cc1, float3(0.667, 0.667, 0.667)*ecd*3.4, x4);
          cc1 = lerp(cc1, ecd*3.4, xE);
          cc1 = lerp(cc1, ecd*3.4, x5);
          cc1 = lerp(cc1, float3(0.667, 0.667, 0.667)*ecd*0.5, x6);	 
          cc1 = lerp(cc1, 0.00, x7);
		  cc1 = lerp(cc1, 0.00, xG);
		  cc1 = lerp(cc1, 0.00, xZ);
          cc1 = lerp(cc1, 0.00, x8);
          cc1*= 1.0;	

   float cc0 = lerp(2.0, 2.0, x1);
         cc0 = lerp(cc0, 3.0, x2);
         cc0 = lerp(cc0, 3.0, x3);
         cc0 = lerp(cc0, 3.0, x4);
         cc0 = lerp(cc0, 3.0, xE);
         cc0 = lerp(cc0, 3.0, x5);
         cc0 = lerp(cc0, 3.0, x6);
         cc0 = lerp(cc0, 3.0, x7);
		 cc0 = lerp(cc0, 3.0, xG);
		 cc0 = lerp(cc0, 3.0, xZ);
         cc0 = lerp(cc0, 2.0, x8);  
  	 
   float xt = saturate(fr0*tex*cc0);
	 
float3 re0 = 15.0*refw;
float3 gc;
float3 gn;
if (wx.x==0,1) gc = sl;
if (wx.y==0,1) gn = sl;
if (wx.x==4) gc = re0;
if (wx.x==7) gc = re0;
if (wx.x==8) gc = re0;
if (wx.x==9) gc = re0;
if (wx.x==12) gc = re0;
if (wx.x==15) gc = re0;
if (wx.x==16) gc = re0;
if (wx.y==4) gn = re0;
if (wx.y==7) gn = re0;
if (wx.y==8) gn = re0;
if (wx.y==9) gn = re0;
if (wx.y==12) gn = re0;
if (wx.y==15) gn = re0;
if (wx.y==16) gn = re0;		
		   
    float3 wmix1 = lerp(gc, gn, wx.z);		

float re1 = refglass*1.3;
float gc1;
float gn1;
if (wx.x==0,1) gc1 = refglass;
if (wx.y==0,1) gn1 = refglass;
if (wx.x==4) gc1 = re1;
if (wx.x==7) gc1 = re1;
if (wx.x==8) gc1 = re1;
if (wx.x==9) gc1 = re1;
if (wx.x==12) gc1 = re1;
if (wx.x==15) gc1 = re1;
if (wx.x==16) gc1 = re1;
if (wx.y==4) gn1 = re1;
if (wx.y==7) gn1 = re1;
if (wx.y==8) gn1 = re1;
if (wx.y==9) gn1 = re1;
if (wx.y==12) gn1 = re1;
if (wx.y==15) gn1 = re1;
if (wx.y==16) gn1 = re1;		
		   
    float wmix2 = lerp(gc1, gn1, wx.z);	
	

	 
	float3	specular=0.0;	
	for (int li=0; li<8; li++)
	{	 	 		 
		float3 sv5 = normalize(lightDirection[li].xyz);
		float specfact = saturate(dot(sv5, wn.xyz));
		      specfact = pow(specfact, 10.0);
		      specular+= saturate(lightDiffuse[li]-lightSpecular[li]) * specfact;		 
	}	

   float cc3 = lerp(0.0, 0.0, x1);
         cc3 = lerp(cc3, 0.3, x2);
         cc3 = lerp(cc3, 1.0, x3);
         cc3 = lerp(cc3, 1.0, x4);
         cc3 = lerp(cc3, 1.0, xE);
         cc3 = lerp(cc3, 1.0, x5);
         cc3 = lerp(cc3, 1.0, x6);
         cc3 = lerp(cc3, 0.0, x7);
		 cc3 = lerp(cc3, 0.0, xG);
		 cc3 = lerp(cc3, 0.0, xZ);
         cc3 = lerp(cc3, 0.0, x8);  	
   
   
    float4 dc1 = diffColor;
	float3 st8 = normalize(dc1.xyz);
	float3 ct8=dc1.xyz/st8.xyz;
	       ct8=pow(ct8, LightingContrast);
		   st8.xyz = pow(st8.xyz, LightingSaturate);
	       dc1.xyz = ct8*st8.xyz;
		   
           r0.xyz+= dc1*specular*tex*6.0;  // Ночное Освещение 1
           r0.xyz+= specular*tex*0.1;  // Ночное Освещение 2	 
           r0.xyz+= (float3(0.902, 1.0, 1.0)*fr3*tex*cc);	 	 
           r0.xyz+= saturate(fr2*dc1*tex*cc2)*2.2;
           r0.xyz+= saturate(fr*dc1*tex*cc1)*2.2;	
		   	   
float ch1 = 2.5*carchrome;
float ch2 = 1.5*carchrome;
float gc4; 
float gn4; 
if (wx.x==0,1) gc4 = ch2;
if (wx.y==0,1) gn4 = ch2;
if (wx.x==4) gc4 = ch1;
if (wx.x==7) gc4 = ch1;
if (wx.x==8) gc4 = ch1;
if (wx.x==9) gc4 = ch1;
if (wx.x==12) gc4 = ch1;
if (wx.x==15) gc4 = ch1;
if (wx.x==16) gc4 = ch1;
if (wx.y==4) gn4 = ch1;
if (wx.y==7) gn4 = ch1;
if (wx.y==8) gn4 = ch1;
if (wx.y==9) gn4 = ch1;
if (wx.y==12) gn4 = ch1;
if (wx.y==15) gn4 = ch1;
if (wx.y==16) gn4 = ch1;		
		   
    float wmix5 = lerp(gc4, gn4, wx.z);		
    float gg = saturate(1.1*mask0*mask3*mask2*mask5)*1.0;

	if (cre==true) r0.xyz = lerp(r0, saturate(0.06*rfl1)*11.9*wmix5, gg);
	//r0.xyz = lerp(r0, lerp(r0+wmix1*tl0*saturate(0.06*rfl1)*1.7, r0, xt*EFresnelFactor), 1.0*mask0*mask4*mask5); // Отражения на тачке
    r0.xyz = lerp(r0, lerp(wmix1*tl0*saturate(0.06*rfl1)*1.7, r0, xt*EFresnelFactor), 1.0*mask0*mask4*mask5); // Отражения на тачке	
    r0.xyz = lerp(r0, r0+(saturate(0.08*rfl1)*3.5)*wmix2, 1.0*mask6); // Отражения на стеклах	
		   
	float3 d1 = r0;
           r0.xyz+= saturate(mask0*gl*0.95*cc3)*4.0;
	float3 d0 = r0;
float3 dcurrent;
float3 dnext;
if (wx.x==0,1) dcurrent = d0;
if (wx.y==0,1) dnext = d0;
if (wx.x==4) dcurrent = d1;
if (wx.x==7) dcurrent = d1;
if (wx.x==8) dcurrent = d1;
if (wx.x==9) dcurrent = d1;
if (wx.x==12) dcurrent = d1;
if (wx.x==15) dcurrent = d1;
if (wx.x==16) dcurrent = d1;
if (wx.y==4) dnext = d1;
if (wx.y==7) dnext = d1;
if (wx.y==8) dnext = d1;
if (wx.y==9) dnext = d1;
if (wx.y==12) dnext = d1;
if (wx.y==15) dnext = d1;
if (wx.y==16) dnext = d1;	  

	r0.xyz = lerp(dcurrent, dnext, wx.z);
	float4 r3 = r0;
	float3 st7 = normalize(r3.xyz);
	float3 ct7=r3.xyz/st7.xyz;
	       ct7=pow(ct7, 1.0);
		   st7.xyz = pow(st7.xyz, 1.0);
	       r3.xyz = ct7*st7.xyz;
	
           r0.xyz = lerp(r0, r0+wLhights*r3*(-0.4), 1.0*mask2*mask5);
           r0.xyz = lerp(r0, r0+mask2*mask5*tcube, 0.02); 
           r0.w = lerp(1.0, opacity, carglass);
	float fadefact = (FogParam.w - IN.vposition.z) / (FogParam.w - FogParam.z);
           r0.xyz = lerp(FogFarColor.xyz, r0.xyz, saturate(fadefact));
	
	float nonlineardepth = (IN.vposition.z/IN.vposition.w);	
	float3 ssnormal = normalize(IN.viewnormal);
	       ssnormal.yz = -ssnormal.yz;
	float4 r1;
	       r1.xyz = ssnormal*0.5+0.5;
	       r1.w = 251.0/255.0;
		   
	float4 r2;	
	       r2 = nonlineardepth;
	       r2.w = 1.0;
	PS_OUTPUT3	OUT;
	OUT.Color[0]=r0;
	OUT.Color[1]=r1;
	OUT.Color[2]=r2;
	return OUT;
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

technique DrawTransparent
{
    pass p0
    {
	VertexShader = compile vs_3_0 VS_Draw();
	PixelShader  = compile ps_3_0 PS_Draw();
	}
}

