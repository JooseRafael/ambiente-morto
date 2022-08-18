//++++++++++++++++++++++++++++++++++++++++++++
// ENBSeries effect file
// visit http://enbdev.com for updates
// Copyright (c) 2007-2013 Boris Vorontsov
//++++++++++++++++++++++++++++++++++++++++++++



//post processing mode. Change value (could be 1, 2, 3, 4). Every mode have own internal parameters, look below
#ifndef POSTPROCESS
 #define POSTPROCESS	2
#endif


float4 WeatherAndTime;

//+++++++++++++++++++++++++++++
//internal parameters, can be modified
//+++++++++++++++++++++++++++++
//modify these values to tweak various color processing
//POSTPROCESS 1
float	EAdaptationMinV1=0.01;
float	EAdaptationMaxV1=0.07;
float	EContrastV1=0.95;
float	EColorSaturationV1=1.0;
float	EToneMappingCurveV1=6.0;

//POSTPROCESS 2
//float	EBrightnessV2=2.5;



float	EAdaptationMinV2
<
	string UIName="EAdaptationMinV2";
	string UIWidget="Spinner";
	float UIMin=-100.0;//not zero!!!
	float UIMax=100.0;
> = {0.05};


float	EAdaptationMaxV2
<
	string UIName="EAdaptationMaxV2";
	string UIWidget="Spinner";
	float UIMin=-100.0;//not zero!!!
	float UIMax=100.0;
> = {0.05};  //0.125

float	EToneMappingCurveV2
<
	string UIName="EToneMappingCurveV2";
	string UIWidget="Spinner";
	float UIMin=-100.0;//not zero!!!
	float UIMax=100.0;
> = {8.00};

float	EIntensityContrastV2
<
	string UIName="EIntensityContrastV2";
	string UIWidget="Spinner";
	float UIMin=-100.0;//not zero!!!
	float UIMax=100.0;
> = {1.00};


float	EColorSaturationV2
<
	string UIName="EColorSaturationV2";
	string UIWidget="Spinner";
	float UIMin=-100.0;//not zero!!!
	float UIMax=100.0;
> = {1.30};

float	EToneMappingOversaturationV2
<
	string UIName="EToneMappingOversaturationV2";
	string UIWidget="Spinner";
	float UIMin=-10000.0;//not zero!!!
	float UIMax=10000.0;
> = {180.00};


//POSTPROCESS 3
float	EAdaptationMinV3=0.05;
float	EAdaptationMaxV3=0.125;
float	EToneMappingCurveV3=4.0;
float	EToneMappingOversaturationV3=60.0;

//POSTPROCESS 4
float	EAdaptationMinV4=0.2;
float	EAdaptationMaxV4=0.125;
float	EBrightnessCurveV4=0.7;
float	EBrightnessMultiplierV4=0.45;
float	EBrightnessToneMappingCurveV4=0.5;



//parameters for ldr color correction, if enabled
float	ECCGamma
<
	string UIName="CC: Gamma";
	string UIWidget="Spinner";
	float UIMin=0.2;//not zero!!!
	float UIMax=6.0;
> = {1.00};

float	ECCInBlack
<
	string UIName="CC: In black";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=2.0;
> = {0.0};

float	ECCInWhite
<
	string UIName="CC: In white";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=2.0;
> = {0.54};

float	ECCOutBlack
<
	string UIName="CC: Out black";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=2.0;
> = {0.01};

float	ECCOutWhite
<
	string UIName="CC: Out white";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=2.0;
> = {1.0};

float	ECCBrightness
<
	string UIName="CC: Brightness";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10.0;
> = {0.93};

float	ECCContrastGrayLevel
<
	string UIName="CC: Contrast gray level";
	string UIWidget="Spinner";
	float UIMin=0.01;
	float UIMax=0.99;
> = {0.30};

float	ECCContrast
<
	string UIName="CC: Contrast";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10.0;
> = {0.94};

float	ECCSaturation
<
	string UIName="CC: Saturation";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=10.0;
> = {1.0};

float	ECCDesaturateShadows
<
	string UIName="CC: Desaturate shadows";
	string UIWidget="Spinner";
	float UIMin=0.0;
	float UIMax=1.0;
> = {0.0};

float3	ECCColorBalanceShadows <
	string UIName="CC: Color balance shadows";
	string UIWidget="Color";
> = {0.33, 0.33, 0.33};

float3	ECCColorBalanceHighlights <
	string UIName="CC: Color balance highlights";
	string UIWidget="Color";
> = {1.0, 1.0, 1.0};

float3	ECCChannelMixerR <
	string UIName="CC: Channel mixer R";
	string UIWidget="Color";
> = {1.0, 0.0, 0.0};

float3	ECCChannelMixerG <
	string UIName="CC: Channel mixer G";
	string UIWidget="Color";
> = {0.0, 1.0, 0.0};

float3	ECCChannelMixerB <
	string UIName="CC: Channel mixer B";
	string UIWidget="Color";
> = {0.0, 0.0, 1.0};



//+++++++++++++++++++++++++++++
//external parameters, do not modify
//+++++++++++++++++++++++++++++
//keyboard controlled temporary variables (in some versions exists in the config file). Press and hold key 1,2,3...8 together with PageUp or PageDown to modify. By default all set to 1.0
float4	tempF1; //0,1,2,3
float4	tempF2; //5,6,7,8
float4	tempF3; //9,0
//x=generic timer in range 0..1, period of 16777216 ms (4.6 hours), w=frame time elapsed (in seconds)
float4	Timer;
//x=Width, y=1/Width, z=ScreenScaleY, w=1/ScreenScaleY
float4	ScreenSize;
//changes in range 0..1, 0 means that night time, 1 - day time
float	ENightDayFactor;
//changes 0 or 1. 0 means that exterior, 1 - interior
float	EInteriorFactor;
//changes in range 0..1, 0 means full quality, 1 lowest dynamic quality (0.33, 0.66 are limits for quality levels)
float	EAdaptiveQualityFactor;
//enb version of bloom applied, ignored if original post processing used
float	EBloomAmount;
//fov in degrees
float	FieldOfView;


texture2D texs0;//color
//texture2D texs1;//ds2 BloomTexSampler
//texture2D texs2;//ds2 LuminanceTexSampler
texture2D texs3;//bloom enb
texture2D texs4;//adaptation enb
//texture2D texs5;//ds2 ToneMapTableSampler
//texture2D texs6;//ds2 LUTTexSampler
texture2D texs7;//palette enb

texture2D texNormal;
texture2D texDepth;

sampler2D _s0 = sampler_state
{
	Texture   = <texs0>;
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

sampler2D _s3 = sampler_state
{
	Texture   = <texs3>;
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

sampler2D _s4 = sampler_state
{
	Texture   = <texs4>;
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

sampler2D _s7 = sampler_state
{
	Texture   = <texs7>;
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

sampler2D SamplerNormal = sampler_state
{
	Texture   = <texNormal>;
	MinFilter = POINT;
	MagFilter = POINT;
	MipFilter = NONE;//NONE;
	AddressU  = Clamp;
	AddressV  = Clamp;
	SRGBTexture=FALSE;
	MaxMipLevel=0;
	MipMapLodBias=0;
};

sampler2D SamplerDepth = sampler_state
{
	Texture   = <texDepth>;
	MinFilter = POINT;
	MagFilter = POINT;
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
	float2 txcoord0 : TEXCOORD0;
};
struct VS_INPUT_POST
{
	float3 pos  : POSITION;
	float2 txcoord0 : TEXCOORD0;
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
VS_OUTPUT_POST VS_Quad(VS_INPUT_POST IN)
{
	VS_OUTPUT_POST OUT;

	OUT.vpos=float4(IN.pos.x,IN.pos.y,IN.pos.z,1.0);

	OUT.txcoord0.xy=IN.txcoord0.xy;

	return OUT;
}

float linearlizeDepth(float nonlinearDepth)
{
float2 dofProj=float2(0.0509804, 3098.0392);
float2 dofDist=float2(0.0, 0.0509804);

float4 depth=nonlinearDepth;


depth.y=-dofProj.x + dofProj.y;
depth.y=1.0/depth.y;
depth.z=depth.y * dofProj.y;
depth.z=depth.z * -dofProj.x;
depth.x=dofProj.y * -depth.y + depth.x;
depth.x=1.0/depth.x;

depth.y=depth.z * depth.x;

depth.x=depth.z * depth.x - dofDist.y;
depth.x+=dofDist.x * -0.5;

depth.x=max(depth.x, 0.0);

return depth.x;
}


float4 PS_GTASA(VS_OUTPUT_POST IN, float2 vPos : VPOS) : COLOR
{
	float4 _oC0=0.0; //output
	float4 _oC1=0.0; //output
	
	float4 _v0=0.0;
	
	float  lindepth = linearlizeDepth(tex2D(SamplerDepth, IN.txcoord0.xy).x)*1.0;
	
        float4 _c6=float4(0, 0, 0, 0);
        float4 _c7=float4(2.212500006, 2.715399981, 2.0720999986, 9.0);

	float4 r0;
	float4 r1;
	float4 r2;
	float4 r3;
	float4 r4;
	float4 r5;
	float4 r6;
	float4 r7;
	float4 r8;
	float4 r9;
	float4 r10;
	float4 r12;
        float4 r13;
        float4 r14;
        float4 r15;

	_v0.xy=IN.txcoord0.xy;
	r1=tex2D(_s0, _v0.xy); //color

	//apply bloom
	float4	xcolorbloom=tex2D(_s3, _v0.xy);

	xcolorbloom.xyz=xcolorbloom-r1;
	xcolorbloom.xyz=max(xcolorbloom, 0.0);
	r1.xyz+=xcolorbloom*EBloomAmount;
        
        r15=r1;
        _oC0.xyz=r1.xyz;
		_oC1.xyz=r1.xyz;
        float4 color=r1;



	//adaptation in time
	float4	Adaptation=tex2D(_s4, 0.5);
	float	grayadaptation=max(max(Adaptation.x, Adaptation.y), Adaptation.z);

#if (POSTPROCESS==1)

	grayadaptation=max(grayadaptation, 0.0);
	grayadaptation=min(grayadaptation, 50.0);
	color.xyz=color.xyz/(grayadaptation*EAdaptationMaxV1+EAdaptationMinV1);//*tempF1.x

	float cgray=dot(color.xyz, float3(0.27, 0.67, 0.06));
	cgray=pow(cgray, EContrastV1);
	float3 poweredcolor=pow(color.xyz, EColorSaturationV1);
	float newgray=dot(poweredcolor.xyz, float3(0.27, 0.67, 0.06));
	color.xyz=poweredcolor.xyz*cgray/(newgray+0.0001);

	float3	luma=color.xyz;
	float	lumamax=300.0;
	color.xyz=(color.xyz * (1.0 + color.xyz/lumamax))/(color.xyz + EToneMappingCurveV1);

#endif



#if (POSTPROCESS==2)

	grayadaptation=max(grayadaptation, 0.0);
	grayadaptation=min(grayadaptation, 50.0);
	color.xyz=color.xyz/(grayadaptation*EAdaptationMaxV2+EAdaptationMinV2);//*tempF1.x

	//color.xyz*=EBrightnessV2;
	color.xyz+=0.000001;
	float3 xncol=normalize(color.xyz);
	float3 scl=color.xyz/xncol.xyz;
	scl=pow(scl, EIntensityContrastV2);
	xncol.xyz=pow(xncol.xyz, EColorSaturationV2);
	color.xyz=scl*xncol.xyz;

	float	lumamax=EToneMappingOversaturationV2;
	color.xyz=(color.xyz * (1.0 + color.xyz/lumamax))/(color.xyz + EToneMappingCurveV2);

#endif


#if (POSTPROCESS==3)

	grayadaptation=max(grayadaptation, 0.0);
	grayadaptation=min(grayadaptation, 50.0);
	color.xyz=color.xyz/(grayadaptation*EAdaptationMaxV3+EAdaptationMinV3);//*tempF1.x

	float	lumamax=EToneMappingOversaturationV3;
	color.xyz=(color.xyz * (1.0 + color.xyz/lumamax))/(color.xyz + EToneMappingCurveV3);

#endif



#if (POSTPROCESS==4)

	grayadaptation=max(grayadaptation, 0.0);
	grayadaptation=min(grayadaptation, 50.0);
	color.xyz=color.xyz/(grayadaptation*EAdaptationMaxV4+EAdaptationMinV4);

	float Y = dot(color.xyz, float3(0.299, 0.587, 0.114)); //0.299 * R + 0.587 * G + 0.114 * B;
	float U = dot(color.xyz, float3(-0.14713, -0.28886, 0.436)); //-0.14713 * R - 0.28886 * G + 0.436 * B;
	float V = dot(color.xyz, float3(0.615, -0.51499, -0.10001)); //0.615 * R - 0.51499 * G - 0.10001 * B;
	Y=pow(Y, EBrightnessCurveV4);
	Y=Y*EBrightnessMultiplierV4;
//	Y=Y/(Y+EBrightnessToneMappingCurveV4);
//	float	desaturatefact=saturate(Y*Y*Y*1.7);
//	U=lerp(U, 0.0, desaturatefact);
//	V=lerp(V, 0.0, desaturatefact);
	color.xyz=V * float3(1.13983, -0.58060, 0.0) + U * float3(0.0, -0.39465, 2.03211) + Y;

	color.xyz=max(color.xyz, 0.0);
	color.xyz=color.xyz/(color.xyz+EBrightnessToneMappingCurveV4);

#endif



	//pallete texture (0.082+ version feature)
#ifdef E_CC_PALETTE
	color.rgb=saturate(color.rgb);
	float3	brightness=Adaptation.xyz;//tex2D(_s4, 0.5);//adaptation luminance
//	brightness=saturate(brightness);//old version from ldr games
	brightness=(brightness/(brightness+1.0));//new version
	brightness=max(brightness.x, max(brightness.y, brightness.z));//new version
	float3	palette;
	float4	uvsrc=0.0;
	uvsrc.y=brightness.r;
	uvsrc.x=color.r;
	palette.r=tex2Dlod(_s7, uvsrc).r;
	uvsrc.x=color.g;
	uvsrc.y=brightness.g;
	palette.g=tex2Dlod(_s7, uvsrc).g;
	uvsrc.x=color.b;
	uvsrc.y=brightness.b;
	palette.b=tex2Dlod(_s7, uvsrc).b;
	color.rgb=palette.rgb;
#endif //E_CC_PALETTE



#ifdef E_CC_PROCEDURAL
	float	tempgray;
	float4	tempvar;
	float3	tempcolor;
/*
	//these replaced by "levels"
	//+++ gamma
	if (ECCGamma!=1.0)
	color=pow(color, 1.0/ECCGamma);

	//+++ brightness like in photoshop
	color=color+ECCAditiveBrightness;

	//+++ lightness
	tempvar.x=saturate(ELightness);
	tempvar.y=saturate(1.0+ECCLightness);
	color=tempvar.x*(1.0-color) + (tempvar.y*color);
*/
	//+++ levels like in photoshop, including gamma, lightness, additive brightness
	color=max(color-ECCInBlack, 0.0) / max(ECCInWhite-ECCInBlack, 0.0001);
	if (ECCGamma!=1.0) color=pow(color, ECCGamma);
	color=color*(ECCOutWhite-ECCOutBlack) + ECCOutBlack;

	//+++ brightness
	color=color*ECCBrightness;

	//+++ contrast
	color=(color-ECCContrastGrayLevel) * ECCContrast + ECCContrastGrayLevel;

	//+++ saturation
	tempgray=dot(color, 0.3333);
	color=lerp(tempgray, color, ECCSaturation);

	//+++ desaturate shadows
	tempgray=dot(color, 0.3333);
	tempvar.x=saturate(1.0-tempgray);
	tempvar.x*=tempvar.x;
	tempvar.x*=tempvar.x;
	color=lerp(color, tempgray, ECCDesaturateShadows*tempvar.x);

	//+++ color balance
	color=saturate(color);
	tempgray=dot(color, 0.3333);
	float2	shadow_highlight=float2(1.0-tempgray, tempgray);
	shadow_highlight*=shadow_highlight;
	color.rgb+=(ECCColorBalanceHighlights*2.0-1.0)*color * shadow_highlight.x;
	color.rgb+=(ECCColorBalanceShadows*2.0-1.0)*(1.0-color) * shadow_highlight.y;

	//+++ channel mixer
	tempcolor=color;
	color.r=dot(tempcolor, ECCChannelMixerR);
	color.g=dot(tempcolor, ECCChannelMixerG);
	color.b=dot(tempcolor, ECCChannelMixerB);
#endif //E_CC_PROCEDURAL



	_oC0.w=1.0;
	_oC0.xyz=color.xyz;

	float4 no1 = tex2D(SamplerNormal, IN.txcoord0.xy);

	return _oC0;
}



//switch between vanilla and mine post processing
technique Shader_GTASA <string UIName="ENBSeries";>
{
	pass p0
	{
		VertexShader = compile vs_3_0 VS_Quad();
		PixelShader  = compile ps_3_0 PS_GTASA();

		ColorWriteEnable=ALPHA|RED|GREEN|BLUE;
		ZEnable=FALSE;
		ZWriteEnable=FALSE;
		CullMode=NONE;
		AlphaTestEnable=FALSE;
		AlphaBlendEnable=FALSE;
		SRGBWRITEENABLE=FALSE;
	}
}



//original shader of post processing
float4	PS_ORIG(VS_OUTPUT_POST IN) : COLOR
{
	float4 _oC0=tex2D(_s0, IN.txcoord0.xy);
	_oC0.w=1.0;
	return _oC0;
}

technique Shader_ORIGINALPOSTPROCESS <string UIName="Vanilla";>
{
	pass p0
	{
		VertexShader = compile vs_3_0 VS_Quad();
		PixelShader  = compile ps_3_0 PS_ORIG();
		ColorWriteEnable=ALPHA|RED|GREEN|BLUE;
		ZEnable=FALSE;
		ZWriteEnable=FALSE;
		CullMode=NONE;
		AlphaTestEnable=FALSE;
		AlphaBlendEnable=FALSE;
		SRGBWRITEENABLE=FALSE;
    }
}

