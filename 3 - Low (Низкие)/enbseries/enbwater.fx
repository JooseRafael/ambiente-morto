
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
float FieldOfView; float GameTime; float4 SunDirection; float4 CustomShaderConstants1[8];
float4 MatrixVP[4]; float4 MatrixInverseVP[4]; float4 MatrixVPRotation[4]; float4 MatrixInverseVPRotation[4];
float4 MatrixView[4]; float4 MatrixInverseView[4]; float4 CameraPosition;
float4x4 MatrixWVP; float4x4 MatrixWVPInverse; float4x4 MatrixWorld; float4x4 MatrixProj;
float4 FogParam; float4 FogFarColor; float4 WaterParameters1; float4 WaterParameters2;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//Textures
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

texture2D texOriginal;
texture2D texRefl;
texture2D texReflPrev;
texture2D texEnv;
texture2D texDepth;
texture2D texNormal;
texture2D texWave < string ResourceName="WaterFoam.png"; >;
texture2D texWater < string ResourceName="WaterRelief.png"; >;
texture2D texAlpha < string ResourceName="WaterAlpha.png"; >;
sampler2D SamplerOriginal = sampler_state { Texture   = <texOriginal>; };

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
//Sampler Inputs
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

sampler2D SamplerAlpha = sampler_state
{
	Texture   = <texAlpha>;
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

sampler2D SamplerWave = sampler_state
{
	Texture   = <texWave>;
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

sampler2D SamplerExternal = sampler_state
{
	Texture   = <texWater>;
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


sampler2D SamplerRefl = sampler_state
{
	Texture   = <texRefl>;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerRefl2 = sampler_state
{
	Texture   = <texReflPrev>;
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
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

struct PS_OUTPUT3
{
	float4 Color[3] : COLOR0;
};

struct VS_INPUT
{
	float3	pos : POSITION;
	float4	diff : COLOR0;
	float4	spec : COLOR1;
	float2	txcoord0 : TEXCOORD0;
};

struct VS_OUTPUT
{
	float4	pos : POSITION;
	float4	diff : COLOR0;
	float2	txcoord0 : TEXCOORD0;
	float4	txcoord1 : TEXCOORD1;
	float3	normal : TEXCOORD2;
	float3	vnormal : TEXCOORD3;
	float4	vposition : TEXCOORD4;
	float3	wposition : TEXCOORD5;
	float3	eyedir : TEXCOORD6;
};

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

bool ECausticsEnable
<
	string UIName="Water: EnableCaustics";
> = {true};

float ECausticsIntensity
<
	string UIName="Water: CausticsIntensity";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10.0;
> = {2.0};

float EWaveBumpiness
<
	string UIName="Water: WaveBumpiness";
	string UIWidget="Spinner";
	float UIMin=0.2;
	float UIMax=3.0;
> = {1.0};

float reflct
<
	string UIName="Reflection: Contrast";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=3.0;
> = {1.0};

float reflst
<
	string UIName="Reflection: Saturate";
	string UIWidget="Spinner";
	float UIMin=0.2;
	float UIMax=3.0;
> = {1.0};

float reflbr
<
	string UIName="Reflection: Brightness";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=3.0;
> = {1.0};

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

VS_OUTPUT VS_Draw(VS_INPUT IN)
{
	VS_OUTPUT	OUT;
	float4 pos = float4(IN.pos.x,IN.pos.y,IN.pos.z,1.0);
	float4 tpos = mul(pos, MatrixWVP);
	float3 cnormal = float3(0.0,0.0,1.0);
	OUT.normal = mul(cnormal, MatrixWorld);
	OUT.vnormal = normalize(mul(cnormal.xyz, MatrixWVP));
	OUT.pos = tpos;
	OUT.vposition = tpos;
	OUT.wposition = mul(pos, MatrixWorld);
	float4 uv;
	uv.xy = IN.txcoord0.xy;
	uv.zw = 0.0;
	OUT.txcoord0 = IN.txcoord0;
	uv.xy = IN.txcoord0.xy+WaterParameters2.w;
	uv.zw = IN.txcoord0.xy-WaterParameters2.w;
	uv.xy = 0.07*pos.xy+WaterParameters2.w;
	OUT.txcoord1 = uv;
	OUT.diff = IN.diff;
	float3 campos = CameraPosition;
	OUT.eyedir=(mul(pos, MatrixWorld) - campos);
 	return OUT;
}

float3 hash1(float3 p3)
{
	p3 = frac(p3 * float3(0.1031, 0.11369, 0.13787));
	p3 += dot(p3, p3.yxz+19.19);
	p3=float3((p3.x+p3.y), (p3.x+p3.z), (p3.y+p3.z)) * p3.zyx;
	return frac(p3);
}

float WorleyNoise(float3 n, float scale)
{
	n *= scale;
	float2	dis=1.0;
	for (float x=-1.0; x<1.01; x+=1.0)
	{
		for (float y=-1.0; y<1.01; y+=1.0)
		{
			for (float z=-1.0; z<1.01; z+=1.0)
			{
				float3	xyz=float3(x,y,z);
				float3	frcn=frac(n);
				float3	p=n-frcn + xyz;
				float3	pp = p;
				if (pp.z < 0.0) pp.z = scale + pp.z;
				if (pp.z > 1.0) pp.z = pp.z - scale;
				float2	sqrdist;
				frcn=xyz - frcn;
				p = hash1(pp) + frcn;
				sqrdist.x=dot(p, p);
				pp.z+=0.012;
				p = hash1(pp) + frcn;
				sqrdist.y=dot(p, p);
				dis.xy = min(dis.xy, sqrdist.xy);
			}
		}
	}
	dis.x=dis.x*dis.y;
	return saturate(sqrt(dis.x));
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
       tv.y = -tv.y;   
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

float4 reflection(float3 n, float2 coord)
{
   float3 wpos = wpd1(coord.xy);
   float3 v = {1.0, 1.0, 4.3};
          n = (normalize(n.xyz*v)*1.0); 
   float3 n0 = reflect(wpos.xyz, n.xyz);		 	 
   float3 n1 = ((1000.0/0.1)*n0)/1000.0; 
   float3 Cam = CameraPosition;			
   float3 r0 = (wpos+n1);
   float2 r1 = wpd2(r0.xyz);	
   float4 r2 = tex2Dlod(SamplerRefl2, float4(r1.xy, 0.0, 0.0));
	      r2.xyz+=0.000001;
   float3 xncol = normalize(r2.xyz);
   float3 scl=r2.xyz/xncol.xyz;
	      scl=pow(scl, 1.0);	
	      xncol.xyz = pow(xncol.xyz, 1.0);
	      r2.xyz = scl*xncol.xyz;
	      r2.xyz*= 1.0;	
	      r2.w = r1.y<0.0||r1.y>1.0 ? 0.0:1.0;		 
   return r2;
}

////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////SA_DirectX/////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

PS_OUTPUT3 PS_Draw(VS_OUTPUT IN, float2 vPos : VPOS, uniform sampler2D ColorTexture)
{
	float4 r0 = 0.0;
	float2 of = 0.0;	
	float4 color;
	float2 coord;
	float2 uv1;
	float2 uv2;
	float2 uv3;
	float2 uv4;
	float2 uv5;	
	float3 wed;
	float3 un0;
	float3 un1;
	float4 origcolor;
	float2 txcoord;
	float3 normal;
	float4 tvec;
	
	txcoord.xy = (IN.vposition.xy /IN.vposition.w)*float2(0.5, -0.5) + 0.5;
	txcoord.xy+= 0.5*float2(ScreenSize.y, ScreenSize.y*ScreenSize.z);
	wed = normalize(-IN.eyedir.xyz);
    float tx = Timer.x*100.0*5.0;
	float ts = 0.06;
	uv1.xy = 0.38*IN.txcoord1+float2(ts* 10.00, -ts*5.70)*tx;
	uv2.xy = 1.80*IN.txcoord1+float2( ts*11.00, -ts*7.75)*tx;
	uv3.xy = 10.0*IN.txcoord1+float2(-ts*5.78,  ts*7.72)*tx;
	uv4.xy = 0.13*IN.txcoord1+float2(ts*8.40, -ts*10.10)*tx;
	uv5.xy = 1.10*IN.txcoord1+float2(-ts*5.40, ts*15.10)*tx;		
	color = tex2D(ColorTexture, uv1);
	color+= tex2D(ColorTexture, uv2);
	color+= tex2D(ColorTexture, uv3);	
	color+= tex2D(ColorTexture, uv4);
	color+= tex2D(ColorTexture, uv5);	
	origcolor = color*0.25;
	
    float4 fx = tex2D(SamplerAlpha, 0.07*0.0010 + uv4);
	float2 uvshift = -IN.eyedir.xy * (origcolor.x*2.0-1.0) * 0.001*0.30;
	uv1.xy+= uvshift;
	uv2.xy+= uvshift;
	uv3.xy+= uvshift;	
	uv4.xy+= uvshift;
	uv5.xy+= uvshift;	
	color = tex2D(ColorTexture, uv1);
	color+= tex2D(ColorTexture, uv2);
	color+= tex2D(ColorTexture, uv3);	
	color+= tex2D(ColorTexture, uv4);
	color+= tex2D(ColorTexture, uv5);	
	origcolor = color*0.25;	
	float2 kernel[4]=
	{
		float2( 1.0, 0.0),
		float2(-1.0, 0.0),
		float2( 0.0, 1.0),
		float2( 0.0,-1.0)
	};
	for (int i=0; i<4; i++)
	{
		coord = kernel[i];
		float4 tc;
		float4 tc2;
		float4 tc3;
		float4 tc4;
		float4 tc5;			
		tc = tex2D(ColorTexture, 0.07*coord*0.0010 + uv1.xy);
		tc2 = tex2D(ColorTexture, 0.07*coord*0.0010 + uv2.xy);
		tc3 = tex2D(ColorTexture, 0.07*coord*0.0010 + uv3.xy);
		tc4 = tex2D(ColorTexture, 0.07*coord*0.0010 + uv4.xy);
		tc5 = tex2D(ColorTexture, 0.07*coord*0.0010 + uv5.xy);		
		tc.x+= tc2.x;
		tc.x+= tc3.x;
		tc.x+= tc4.x;	
		tc.x+= tc5.x;			
		of+= coord*(color.x-tc.x);
	}
	un0.xy = of.xy*10.0*2.20*EWaveBumpiness;
	un0.z = sqrt(1.0-dot(of.xy, of.xy));
	un0 = normalize(mul(un0, MatrixWorld));	
	un1.xy = of.xy * 10.0*5.50*EWaveBumpiness;
	un1.z = sqrt(1.0-dot(of.xy, of.xy));
	un1 = normalize(mul(un1, MatrixWorld));
	tvec.xyz = un0;
	normal.x = dot(tvec.xyz, MatrixView[0]);
	normal.y = dot(tvec.xyz, MatrixView[1]);
	normal.z = dot(tvec.xyz, MatrixView[2]);
	normal = normalize(normal);
	of.xy = normal.xy;
	of.y = -of.y;
	of.xy*= 0.05*float2(0.5, 0.0);
	float objdepth = IN.vposition.z/IN.vposition.w;
	float nonlinearobjdepth = objdepth;
	      objdepth = 1.0/max(1.0-objdepth, 0.000000001);
	float planardepth = tex2D(SamplerDepth, txcoord.xy).r;
	float depth = tex2D(SamplerDepth, txcoord.xy).r;

	if (nonlinearobjdepth>depth)
	{
		depth = planardepth;
		of.xy = 0.0;
	}

	float scenedepth = depth;
	      planardepth = 1.0/max(1.0-planardepth,0.000000001);
	      depth = 1.0/max(1.0-depth,0.000000001);
	float depthfact = (depth-objdepth)*WaterParameters1.z;
	      depthfact = depthfact*WaterParameters1.x;
	      depthfact = depthfact/(depthfact*0.50+1.0);

	if (scenedepth>0.99999) depthfact = 1.0;

	float backside = saturate(WaterParameters1.w*10.0);
	      depthfact*= 1.0-backside;
	float shorefade = (planardepth-objdepth);
	      shorefade = shorefade*WaterParameters1.x;
	      shorefade = saturate(13.40*shorefade*1000.0 - 0.05);
	      depthfact*= shorefade;
	
    float2 of2 = of;	
	       of2.xy*= shorefade*0.0;
	       of.xy*= shorefade;

	tvec.xyz = un0;
	tvec.xy*= shorefade;
	normal.x = dot(tvec.xyz, MatrixView[0]);
	normal.y = dot(tvec.xyz, MatrixView[1]);
	normal.z = dot(tvec.xyz, MatrixView[2]);
	normal = normalize(normal);
	un0 = normalize(tvec.xyz);
	tvec.xyz = un1;
	un1 = normalize(tvec.xyz);	
	
	float3 wpos = wpd1(txcoord.xy);
    float3 v = {1.0, 1.0, 4.30};
	float3 e4 = (normalize(un0.xyz*v)*1.0);
	float3 n = reflect(wpos.xyz, e4.xyz);		 
	float3 n0 = ((1000.0/0.1)*n)/1000.0;		
	float3 rd0 = (wpos+n0);
	float2 rd1 = wpd2(rd0.xyz);	
	float nf0 = saturate(16.0*(rd1.y));
	      nf0 = pow(nf0, 1.4);		

    float4 tv;
	       tv.xy = (txcoord.xy+of2.xy*0.2)*2.0-1.0;
	       tv.y = -tv.y;
	       tv.z = scenedepth;
	       tv.w = 1.0;
	float4 wp;		   
	       wp.x = dot(tv, MatrixInverseVPRotation[0]);
	       wp.y = dot(tv, MatrixInverseVPRotation[1]);
	       wp.z = dot(tv, MatrixInverseVPRotation[2]);
	       wp.w = dot(tv, MatrixInverseVPRotation[3]);
	       wp.xyz/= wp.w;
	       wp.xyz+= CameraPosition;

	float deepfade = max(IN.wposition.z-wp.z, 0.0);
		  deepfade*= 0.1;
		  deepfade = deepfade/(deepfade+1.0);

	      tv.xyz = wp-CameraPosition;
	float distfact = dot(tv.xyz, tv.xyz);
	      distfact = saturate(1.0-distfact*0.000025);

	float3 causticspos = wp;
	       causticspos*= 0.25;
	       causticspos.z*= 0.25;
	       causticspos.z+= Timer.x*1677.7216 * 2.0;
	       causticspos.z = frac(causticspos.z);

	float caustic = 0.0;
	float usecaustic = 0.0;
	if (ECausticsEnable==true) usecaustic = 1.0;
	      usecaustic*= distfact;
	      usecaustic*= 0.9-deepfade;
	if (usecaustic>=0.1)
	{
		caustic = WorleyNoise(causticspos.xyz, 8.0);
	}

	caustic = caustic*caustic;
	caustic*= lerp(1.0, caustic, deepfade);
	caustic = caustic-0.1;
	caustic*= shorefade;
	caustic*= 1.0-deepfade;
	caustic*= saturate((WaterParameters1.w+CameraPosition.z - wp.z));
	caustic*= distfact;
	caustic*= ECausticsIntensity;

	float reflamount = saturate(1.0-dot(un0, wed));
	      reflamount*= reflamount;
	      reflamount*= reflamount;	
	float dfx = lerp(1.7, 0.95, pow(ENightDayFactor, 15.0));	
	      reflamount = pow(reflamount, dfx);	
	float reflamount2= saturate(1.0-dot(un1, wed));
	      reflamount2*= reflamount2;
	      reflamount2*= reflamount2;	
	      reflamount2*= shorefade;
	      reflamount2 = reflamount2*1.0;	
	      reflamount*= shorefade;
	      reflamount = reflamount*1.0;
	
	float t = GameTime;
	float x1 = smoothstep(0.0, 4.0, t);
	float x2 = smoothstep(4.0, 5.0, t);
	float x3 = smoothstep(5.0, 6.0, t);
	float x4 = smoothstep(6.0, 7.0, t);
	float xE = smoothstep(8.0, 11.0, t);
	float x5 = smoothstep(16.0, 17.0, t);  
	float x50 = smoothstep(13.0, 14.0, t);  
	float x51 = smoothstep(14.0, 17.95, t);  
	float x60 = smoothstep(17.98, 18.0, t);     
	float x6 = smoothstep(18.0, 19.0, t);
	float x7 = smoothstep(19.0, 20.0, t);
	float xG = smoothstep(20.0, 21.0, t);  
	float xZ = smoothstep(21.0, 22.0, t);
	float x8 = smoothstep(22.0, 23.0, t);
	float x9 = smoothstep(23.0, 24.0, t);
   
   float3 t1 = lerp(0.28, 0.28, x1);
          t1 = lerp(t1, 0.7, x2);
          t1 = lerp(t1, 1.0, x3);
          t1 = lerp(t1, 1.0, x4);
          t1 = lerp(t1, 1.0, xE);
          t1 = lerp(t1, 1.0, x5);
          t1 = lerp(t1, 0.8, x6);		 
          t1 = lerp(t1, 0.4, x7);
		  t1 = lerp(t1, 0.3, xG);
		  t1 = lerp(t1, 0.2, xZ);
          t1 = lerp(t1, 0.2, x8);  
          t1 = lerp(t1, 0.28, x9);   
		  
   float3 t2 = lerp(0.04, 0.05, x1);
          t2 = lerp(t2, 0.1, x2);
          t2 = lerp(t2, 1.3, x3);
          t2 = lerp(t2, 1.6, x4);
          t2 = lerp(t2, 1.6, xE);
          t2 = lerp(t2, 1.6, x5);
          t2 = lerp(t2, 1.6, x6);		 
          t2 = lerp(t2, 0.6, x7);
		  t2 = lerp(t2, 0.5, xG);
		  t2 = lerp(t2, 0.2, xZ);
          t2 = lerp(t2, 0.1, x8);  
          t2 = lerp(t2, 0.04, x9); 

	float4 env = tex2D(SamplerEnv, txcoord.xy);	
	       env.xyz+= 0.000001;	
	float3 st0 = normalize(env.xyz);
	float3 ct0 = env.xyz/st0.xyz;
	       ct0 = pow(ct0, 1.3);
	       st0.xyz = pow(st0.xyz, 1.0);
	       env.xyz = ct0*st0.xyz;	
	       env.xyz*= 3.0;		  
		
float4 wx = WeatherAndTime;		
float3 wc0;
float3 wn0;
float3 cday = env*0.6;
float3 cweather = t2;

if (wx.x==0,1) wc0 = cday;
if (wx.y==0,1) wn0 = cday;
if (wx.x==4) wc0 = cweather;
if (wx.x==7) wc0 = cweather;
if (wx.x==8) wc0 = cweather;
if (wx.x==9) wc0 = cweather;
if (wx.x==12) wc0 = cweather;
if (wx.x==15) wc0 = cweather;
if (wx.x==16) wc0 = cweather;
if (wx.y==4) wn0 = cweather;
if (wx.y==7) wn0 = cweather;
if (wx.y==8) wn0 = cweather;
if (wx.y==9) wn0 = cweather;
if (wx.y==12) wn0 = cweather;
if (wx.y==15) wn0 = cweather;
if (wx.y==16) wn0 = cweather;  

    float3 wmix0 = lerp(wc0, wn0, wx.z);
    float3 refl0 = 4.5*saturate(wmix0*0.2);
	float4 refl = reflection(normalize(un0.xyz), txcoord.xy+of.xy);
	float skyreflmix;		 
	      skyreflmix = saturate(wed.z)*saturate(wed.z);   	   
		   refl.xyz = lerp(refl0, 8.0*saturate(refl*0.25),  pow(0.01*0.5, skyreflmix));
		   refl0 = lerp(refl0, refl.xyz, nf0*refl.w);
		   
	float3 st1 = normalize(refl0.xyz);
	float3 ct2 = refl0.xyz/st1.xyz;
	       ct2 = pow(ct2, reflct);
	       st1.xyz = pow(st1.xyz, reflst);
	       refl0.xyz = ct2*st1.xyz;	
	       refl0.xyz*= reflbr;
		   
	float stopfact = saturate(-dot(un0, wed));
	float anglefact = 1.0-saturate(dot(IN.normal.xyz, -wed) * 3.0);
	anglefact*= anglefact;
	anglefact*= anglefact;
	stopfact = saturate(stopfact-anglefact);
	stopfact = saturate((0.1-stopfact)*10.0);
	stopfact*= backside;

	float alphafix = saturate((origcolor.a-0.95)*50.0);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	r0 = tex2D(SamplerRefl, txcoord.xy);		
    r0.xyz*= (1.0+caustic);
    float3 t3 = lerp(0.5, 0.5, x1);
           t3 = lerp(t3, 0.2, x2);
           t3 = lerp(t3, 0.7, x3);
           t3 = lerp(t3, 0.7, x4);
           t3 = lerp(t3, 1.0, xE);
           t3 = lerp(t3, 1.0, x5);
           t3 = lerp(t3, 0.7, x6);		 
           t3 = lerp(t3, 0.3, x7);
		   t3 = lerp(t3, 0.2, xG);
		   t3 = lerp(t3, 0.2, xZ);
           t3 = lerp(t3, 0.8, x8); 
           t3 = lerp(t3, 0.4, x9);	
    float3 t6 = lerp(0.1, 0.1, x1);
           t6 = lerp(t6, 0.2, x2);
           t6 = lerp(t6, 0.7, x3);
           t6 = lerp(t6, 0.7, x4);
           t6 = lerp(t6, 1.0, xE);
           t6 = lerp(t6, 1.0, x5);
           t6 = lerp(t6, 0.7, x6);		 
           t6 = lerp(t6, 0.3, x7);
		   t6 = lerp(t6, 0.2, xG);
		   t6 = lerp(t6, 0.2, xZ);
           t6 = lerp(t6, 0.2, x8); 
           t6 = lerp(t6, 0.1, x9);	   		   
    r0.xyz*= t3;  
	r0.xyz=lerp(r0.xyz, r0.xyz*float3(0.0392, 0.0784, 0.098), saturate(depthfact));
    float4 res3 = tex2D(SamplerWave, uv2.xy*2.0);
    float4 res3x = tex2D(SamplerWave, uv5.xy*1.0);		 		
    float3 t0 = lerp(0.1, 0.1, x1);
           t0 = lerp(t0, 0.3, x2);
           t0 = lerp(t0, 0.7, x3);
           t0 = lerp(t0, 1.0, x4);
           t0 = lerp(t0, 1.0, xE);
           t0 = lerp(t0, 1.0, x5);
           t0 = lerp(t0, 0.4, x6);		 
           t0 = lerp(t0, 0.3, x7);
           t0 = lerp(t0, 0.2, xG);
           t0 = lerp(t0, 0.2, xZ);
           t0 = lerp(t0, 0.1, x8);  
           t0 = lerp(t0, 0.1, x9);				  		  
    float3 t5 = lerp(0.1, 0.1, x1);
           t5 = lerp(t5, 0.5, x2);
           t5 = lerp(t5, 1.3, x3);
           t5 = lerp(t5, 1.0, x4);
           t5 = lerp(t5, 1.0, xE);
           t5 = lerp(t5, 1.0, x5);
           t5 = lerp(t5, 1.0, x6);		 
           t5 = lerp(t5, 0.7, x7);
           t5 = lerp(t5, 0.5, xG);
           t5 = lerp(t5, 0.4, xZ);
           t5 = lerp(t5, 0.2, x8);  
           t5 = lerp(t5, 0.1, x9);		  
float3 wc1;
float3 wn1;		  
float3 pday = float3(0.471, 0.902, 0.902);
float3 pweather = t5;
if (wx.x==0,1) wc1 = pday;
if (wx.y==0,1) wn1 = pday;
if (wx.x==4) wc1 = pweather;
if (wx.x==7) wc1 = pweather;
if (wx.x==8) wc1 = pweather;
if (wx.x==9) wc1 = pweather;
if (wx.x==12) wc1 = pweather;
if (wx.x==15) wc1 = pweather;
if (wx.x==16) wc1 = pweather;
if (wx.y==4) wn1 = pweather;
if (wx.y==7) wn1 = pweather;
if (wx.y==8) wn1 = pweather;
if (wx.y==9) wn1 = pweather;
if (wx.y==12) wn1 = pweather;
if (wx.y==15) wn1 = pweather;
if (wx.y==16) wn1 = pweather;	  

    float3 wmix1 = lerp(wc1, wn1, wx.z);			  
	float4 resOrig = tex2D(SamplerRefl, txcoord.xy);		
	r0.xyz = lerp(r0*wmix1*t0, t0*float3(0.0392, 0.0784, 0.098)*0.2, saturate(1.50*depthfact*WaterParameters1.y));		
	r0.xyz = lerp(lerp(0.5*resOrig*t3, r0, saturate(3.5*depthfact*WaterParameters1.y)), r0, 0.5);	
	r0.xyz = lerp(lerp(lerp(r0, float3(1.0, 1.0, 1.0)*t6*4.0, res3x*res3), r0, saturate(16.00*depthfact*WaterParameters1.y)), r0,  0.5);	
	r0.xyz = lerp(r0.xyz, origcolor.xyz*float3(0.0392, 0.0784, 0.098), stopfact);		
	r0.w = lerp(IN.diff.w*origcolor.a, 1.0, alphafix);
	r0.xyz = lerp(origcolor.xyz*float3(0.0392, 0.0784, 0.098), r0.xyz, alphafix);	
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	

    float3 sv0 = SunDirection.xyz;
    float3 sv1 = normalize(float3(0.0, 0.0, 0.0));
    float yy1 = smoothstep(0.0, 1.0, t);	
    float yy2 = smoothstep(1.0, 23.0, t);
    float yy3 = smoothstep(23.0, 24.0, t);
   
    float3 sv3 = lerp(sv1, sv0, yy1);
           sv3 = lerp(sv3, sv0, yy2); 
           sv3 = lerp(sv3, sv1, yy3); 

    float3 tcolor = lerp(0.0, 0.0, x1);
           tcolor = lerp(tcolor, 0.0, x2);
           tcolor = lerp(tcolor, 0.0, x3);
           tcolor = lerp(tcolor, float3(0.0392, 0.0784, 0.098)*0.3, x4);
           tcolor = lerp(tcolor, float3(0.0392, 0.0784, 0.098), xE);
           tcolor = lerp(tcolor, float3(0.0392, 0.0784, 0.098), x5);
           tcolor = lerp(tcolor, 0.0, x6);		 
           tcolor = lerp(tcolor, 0.0, x7);
           tcolor = lerp(tcolor, 0.0, xG);
		   tcolor = lerp(tcolor, 0.0, xZ);
           tcolor = lerp(tcolor, 0.0, x8);  
           tcolor = lerp(tcolor, 0.0, x9);
    float3 tcolor2 = lerp(0.0, 0.0, x1);
           tcolor2 = lerp(tcolor2, 0.0, x2);
           tcolor2 = lerp(tcolor2, 0.05, x3);
           tcolor2 = lerp(tcolor2, 0.05, x4);
           tcolor2 = lerp(tcolor2, 0.05, xE);
           tcolor2 = lerp(tcolor2, 0.05, x5);
           tcolor2 = lerp(tcolor2, 0.05, x6);		 
           tcolor2 = lerp(tcolor2, 0.05, x7);
		   tcolor2 = lerp(tcolor2, 0.05, xG);
		   tcolor2 = lerp(tcolor2, 0.05, xZ);
           tcolor2 = lerp(tcolor2, 0.05, x8);  
           tcolor2 = lerp(tcolor2, 0.0, x9);
 
float3 wc2;
float3 wn2; 
if (wx.x==0,1) wc2 = tcolor;
if (wx.y==0,1) wn2 = tcolor;
if (wx.x==4) wc2 = tcolor2;
if (wx.x==7) wc2 = tcolor2;
if (wx.x==8) wc2 = tcolor2;
if (wx.x==9) wc2 = tcolor2;
if (wx.x==12) wc2 = tcolor2;
if (wx.x==15) wc2 = tcolor2;
if (wx.x==16) wc2 = tcolor2;
if (wx.y==4) wn2 = tcolor2;
if (wx.y==7) wn2 = tcolor2;
if (wx.y==8) wn2 = tcolor2;
if (wx.y==9) wn2 = tcolor2;
if (wx.y==12) wn2 = tcolor2;
if (wx.y==15) wn2 = tcolor2;
if (wx.y==16) wn2 = tcolor2;	  

    float3 wmix2 = lerp(wc2, wn2, wx.z);	
           r0.xyz+= lerp(wmix2, 0.5*refl0, (reflamount));		   
    float3 sv4 = normalize(sv3.xyz+wed.xyz);
    float3 n2 = normalize(un1.xyz);
    float factor = 0.01 - dot(-sv4, n2);
          factor = pow(factor, 100.0*1.58);
		  factor*= pow(reflamount2*reflamount2*reflamount2*reflamount2*1.0, 0.05);
    float fr = (factor*factor); 
          fr/= 1.0;
    float3 t4 = lerp(0.0, 0.0, x1);
           t4 = lerp(t4, 0.0, x2);
           t4 = lerp(t4, float3(0.431, 0.275, 0.196)*0.7, x3);
           t4 = lerp(t4, float3(0.431, 0.275, 0.196)*0.6, x4);
           t4 = lerp(t4, float3(0.216, 0.216, 0.216), xE);
           t4 = lerp(t4, float3(0.216, 0.216, 0.216), x50);
           t4 = lerp(t4, float3(0.216, 0.216, 0.216)*0.5,x51);	
           t4 = lerp(t4, float3(0.431, 0.275, 0.196)*0.6, x60);		  
           t4 = lerp(t4, float3(0.431, 0.275, 0.196)*0.0, x6);		 
           t4 = lerp(t4, 0.0, x7);
           t4 = lerp(t4, 0.0, x8);  
           t4 = lerp(t4, 0.0, x9);
float3 wc3;
float3 wn3;		  
float3 ScolorDay = t4*4.0;
float3 ScolorPasm = 0.0;
if (wx.x==0,1) wc3 = ScolorDay;
if (wx.y==0,1) wn3 = ScolorDay;
if (wx.x==4) wc3 = ScolorPasm;
if (wx.x==7) wc3 = ScolorPasm;
if (wx.x==8) wc3 = ScolorPasm;
if (wx.x==9) wc3 = ScolorPasm;
if (wx.x==12) wc3 = ScolorPasm;
if (wx.x==15) wc3 = ScolorPasm;
if (wx.x==16) wc3 = ScolorPasm;
if (wx.y==4) wn3 = ScolorPasm;
if (wx.y==7) wn3 = ScolorPasm;
if (wx.y==8) wn3 = ScolorPasm;
if (wx.y==9) wn3 = ScolorPasm;
if (wx.y==12) wn3 = ScolorPasm;
if (wx.y==15) wn3 = ScolorPasm;
if (wx.y==16) wn3 = ScolorPasm;	  
    float3 wmix3 = lerp(wc3, wn3, wx.z);			  
    r0.xyz+= saturate(fr*wmix3)*1.8;	
	fx.xyz+= 0.000001;	
	float3 st4 = normalize(fx.xyz);
	float3 ct4 = fx.xyz/st4.xyz;
	       ct4 = pow(ct4, 5.50);
	       st4.xyz = pow(st4.xyz, 3.93);
	       fx.xyz = ct4*st4.xyz;
           fx.xyz*= 0.23;
    float af0 = saturate(fx);
	r0.xyz = lerp(r0, lerp(r0, float3(1.0, 1.0, 1.0)*t6*3.0, res3x*res3), af0); 
	float nonlineardepth=(IN.vposition.z/IN.vposition.w);
	float3 ssnormal=normal;
	       ssnormal.yz=-ssnormal.yz;
	float4 r1;
	       r1.xyz = ssnormal*0.5+0.5;	
	       r1.w = 85.0/255.0;
	float4 r2;
	       r2 = nonlineardepth*1.0;
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
	PixelShader  = compile ps_3_0 PS_Draw(SamplerExternal);

	}
}

technique DrawUnderwater
{
	pass p0
    {
	VertexShader = compile vs_3_0 VS_Draw();
	PixelShader  = compile ps_3_0 PS_Draw(SamplerExternal);
	}
}

