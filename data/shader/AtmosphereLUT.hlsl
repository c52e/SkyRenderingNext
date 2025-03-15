#include "AtmosphereCommon.hlsl"

/*
This code contains portions of source code from
https://github.com/ebruneton/precomputed_atmospheric_scattering/blob/master/atmosphere/functions.glsl
Copyright (c) 2017 Eric Bruneton
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
*/

DECLARE_ROOT(AtmosphereRoot)

void GetRMuFromTransmittanceTextureIndex(int2 index, int2 size, out float r, out float mu)
{
    float2 uv = float2(index) / float2(size - 1);
    float x_mu = uv.x;
    float x_r = uv.y;
    // Distance to top atmosphere boundary for a horizontal ray at ground level.
    float H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon, from which we can compute r:
    float rho = H * x_r;
    r = sqrt(rho * rho + bottom_radius * bottom_radius);
    // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
    // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon) -
    // from which we can recover mu:
    float d_min = top_radius - r;
    float d_max = rho + H;
    float d = d_min + x_mu * (d_max - d_min);
    mu = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * r * d);
    mu = ClampCosine(mu);
}

float3 ComputeTransmittanceToTopAtmosphereBoundary(float r, float mu)
{
    const float SAMPLE_COUNT = transmittance_steps;

    float dx = DistanceToTopAtmosphereBoundary(r, mu) / SAMPLE_COUNT;

    float3 optical_length = (float3) (0.0);
    for (float i = 0.5; i < SAMPLE_COUNT; ++i)
    {
        float d_i = i * dx;
        // Distance between the current sample point and the planet center.
        float r_i = sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);

        float altitude_i = r_i - bottom_radius;

        optical_length += GetExtinctionCoefficient(altitude_i) * dx;
    }
#if STORE_OPTICAL_LENGTH_IN_TRANSMITTANCE_TEXTURE
    return optical_length;
#else
    return exp(-optical_length);
#endif
}

[numthreads(8, 8, 1)]
void CSTransmittance(int3 ThreadID : SV_DispatchThreadID)
{
    DECLARE_RESOURCE(RWTexture2D<float4>, transmittanceRW)
    
    int2 dim;
    transmittanceRW.GetDimensions(dim.x, dim.y);
    
    float r, mu;
    GetRMuFromTransmittanceTextureIndex(ThreadID.xy, dim, r, mu);
    float3 transmittance = ComputeTransmittanceToTopAtmosphereBoundary(r, mu);
    transmittanceRW[ThreadID.xy] = float4(transmittance, 1.0);
}



void GetScatteringCoefficient(float altitude, out float3 rayleigh, out float3 mie)
{
    rayleigh = rayleigh_scattering * clamp(exp(-altitude * inv_rayleigh_exponential_distribution), 0, 1);
    mie = mie_scattering * clamp(exp(-altitude * inv_mie_exponential_distribution), 0, 1);
}

#define MULTISCATTERING_COMPUTE_PROGRAM

float3 ComputeScatteredLuminance(sampler2D transmittance_texture
#ifndef MULTISCATTERING_COMPUTE_PROGRAM
    , sampler2D multiscattering_texture
    , float start_i
#endif
    , float3 earth_center, float3 start_position, float3 view_direction, float3 sun_direction,
    float marching_distance, float steps, out float3 transmittance
#ifdef MULTISCATTERING_COMPUTE_PROGRAM
    , out float3 L_f
#endif
)
{
    float r = length(start_position - earth_center);
    float rmu = dot(view_direction, start_position - earth_center);
    float mu_s = dot(view_direction, sun_direction);

    const float SAMPLE_COUNT = steps;
    float dx = marching_distance / SAMPLE_COUNT;

    transmittance = (float3)(1.0);
    float3 luminance = (float3)(0.0);
#ifdef MULTISCATTERING_COMPUTE_PROGRAM
    L_f = (float3)(0.0);
    float start_i = 0.5;
    float rayleigh_phase = IsotropicPhaseFunction();
    float mie_phase = IsotropicPhaseFunction();
#else
    float rayleigh_phase = RayleighPhaseFunction(mu_s);
    float mie_phase = MiePhaseFunction(mu_s, mie_phase_g);
#endif
    for (float i = start_i; i < SAMPLE_COUNT; ++i)
    {
        float d_i = i * dx;
        // Distance between the current sample point and the planet center.
        float r_i = sqrt(d_i * d_i + 2.0 * rmu * d_i + r * r);
        float3 position_i = start_position + view_direction * d_i;
        float altitude_i = r_i - bottom_radius;

        float3 rayleigh_scattering_i;
        float3 mie_scattering_i;
        GetScatteringCoefficient(altitude_i, rayleigh_scattering_i, mie_scattering_i);
        float3 scattering_i = rayleigh_scattering_i + mie_scattering_i;
        float3 scattering_with_phase_i = rayleigh_scattering_i * rayleigh_phase + mie_scattering_i * mie_phase;

        float3 extinction_i = GetExtinctionCoefficient(altitude_i);
        float3 transmittance_i = exp(-extinction_i * dx);
        float3 up_direction_i = normalize(position_i - earth_center);
        float mu_s_i = dot(sun_direction, up_direction_i);
        float3 luminance_i = scattering_with_phase_i * GetSunVisibility(transmittance_texture, r_i, mu_s_i);
#ifndef MULTISCATTERING_COMPUTE_PROGRAM
        float3 multiscattering_contribution = GetMultiscatteringContribution(
            multiscattering_texture, r_i, mu_s_i);
        luminance_i += multiscattering_contribution * scattering_i;
        luminance_i *= solar_illuminance;
#endif

        luminance += transmittance * (luminance_i - luminance_i * transmittance_i) / extinction_i;
#ifdef MULTISCATTERING_COMPUTE_PROGRAM
        L_f += transmittance * (scattering_i - scattering_i * transmittance_i) / extinction_i;
#endif
        transmittance *= transmittance_i;
    }
    return luminance;
}

void GetAltitudeMuSFromMultiscatteringTextureIndex(int2 index, int2 size, out float altitude, out float mu_s)
{
    float2 uv = float2(index) / float2(size - 1);
    float x_mu_s = uv.x;
    float x_altitude = uv.y;
    altitude = x_altitude * (top_radius - bottom_radius);
    mu_s = x_mu_s * 2.0 - 1.0;
}

float3 GetDirectionFromLocalIndex(int index)
{
    float unit_theta = (0.5 + float(index / 8)) / 8.0;
    float unit_phi = (0.5 + float(index % 8)) / 8.0;
    // Uniformly sample on a sphere
    float cos_theta = 1.0 - 2.0 * unit_theta;
    float sin_theta = sqrt(clamp(1.0 - cos_theta * cos_theta, 0, 1));

    float phi = 2 * PI * unit_phi;
    float cos_phi = cos(phi);
    float sin_phi = sin(phi);
    return float3(cos_phi * sin_theta, cos_theta, sin_phi * sin_theta);
}

groupshared float3 L_2nd_order_shared[64];
groupshared float3 f_ms_shared[64];

[numthreads(1, 1, 64)]
void CSMultiscattering(int3 ThreadID : SV_DispatchThreadID)
{
    DECLARE_RESOURCE(RWTexture2D<float4>, multiscatteringRW);
    DECLARE_RESOURCE_INDIRECT(Texture2D<float4>, transmittanceTex, cGlobal.atmosphere.transmittanceTex);
    
    int2 dim;
    multiscatteringRW.GetDimensions(dim.x, dim.y);
    
    float altitude, mu_s;
    GetAltitudeMuSFromMultiscatteringTextureIndex(ThreadID.xy, dim, altitude, mu_s);
    float3 earth_center = float3(0, -bottom_radius, 0);
    float3 start_position = float3(0, altitude, 0);
    float3 sun_direction = float3(0, mu_s, sqrt(1 - mu_s * mu_s));
    int local_index = ThreadID.z;
    float3 view_direction = GetDirectionFromLocalIndex(local_index);
    
    float r = altitude + bottom_radius;
    float mu = view_direction.y;
    bool intersect_bottom = RayIntersectsGround(r, mu);
    float marching_distance = DistanceToNearestAtmosphereBoundary(r, mu, intersect_bottom);

    float3 transmittance;
    float3 L_f;
    float3 luminance = ComputeScatteredLuminance(transmittanceTex, earth_center, start_position,
        view_direction, sun_direction, marching_distance, multiscattering_steps, transmittance, L_f);

    if (intersect_bottom)
    {
        float3 ground_position = start_position + view_direction * marching_distance;
        luminance += transmittance * ComputeGroundLuminance(
            transmittanceTex, earth_center, ground_position, sun_direction);
    }
    
    L_2nd_order_shared[local_index] = luminance;
    f_ms_shared[local_index] = L_f;

    GroupMemoryBarrierWithGroupSync();
    if (local_index < 32)
    {
        L_2nd_order_shared[local_index] += L_2nd_order_shared[local_index + 32];
        f_ms_shared[local_index] += f_ms_shared[local_index + 32];
    }
    GroupMemoryBarrierWithGroupSync();
    if (local_index < 16)
    {
        L_2nd_order_shared[local_index] += L_2nd_order_shared[local_index + 16];
        f_ms_shared[local_index] += f_ms_shared[local_index + 16];
    }
    GroupMemoryBarrierWithGroupSync();
    if (local_index < 8)
    {
        L_2nd_order_shared[local_index] += L_2nd_order_shared[local_index + 8];
        f_ms_shared[local_index] += f_ms_shared[local_index + 8];
    }
    GroupMemoryBarrierWithGroupSync();
    if (local_index < 4)
    {
        L_2nd_order_shared[local_index] += L_2nd_order_shared[local_index + 4];
        f_ms_shared[local_index] += f_ms_shared[local_index + 4];
    }
    GroupMemoryBarrierWithGroupSync();
    if (local_index < 2)
    {
        L_2nd_order_shared[local_index] += L_2nd_order_shared[local_index + 2];
        f_ms_shared[local_index] += f_ms_shared[local_index + 2];
    }
    GroupMemoryBarrierWithGroupSync();
    if (local_index < 1)
    {
        L_2nd_order_shared[local_index] += L_2nd_order_shared[local_index + 1];
        f_ms_shared[local_index] += f_ms_shared[local_index + 1];
    }
    GroupMemoryBarrierWithGroupSync();
    if (local_index > 0)
        return;

    float3 L_2nd_order = L_2nd_order_shared[0] / 64;
    float3 f_ms = f_ms_shared[0] / 64;

    float3 F_ms = 1 / (1 - f_ms);
    float3 multiscattering_contribution = L_2nd_order * F_ms;

    multiscatteringRW[ThreadID.xy] = float4(multiscattering_contribution, 1.0);
}


float2 GetMultiscatteringTextureUvFromRMuS(int2 size, float r, float mu_s)
{
    float x_mu_s = mu_s * 0.5 + 0.5;
    float x_r = (r - bottom_radius) / (top_radius - bottom_radius);
    return GetTextureCoordFromUnitRange(float2(x_mu_s, x_r), size);
}

float3 GetMultiscatteringContribution(sampler2D multiscattering_texture, float r, float mu_s)
{
    float2 uv = GetMultiscatteringTextureUvFromRMuS(textureSize(multiscattering_texture), r, mu_s);
    return texture(multiscattering_texture, uv).rgb;
}

void GetRMuMuSNuFromScatteringTextureUvwz(float4 uvwz,
    out float r, out float mu, out float mu_s, out float nu, out bool ray_r_mu_intersects_ground)
{
    // Distance to top atmosphere boundary for a horizontal ray at ground level.
    float H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon.
    float rho =
        H * GetUnitRangeFromTextureCoord(uvwz.w, SCATTERING_TEXTURE_R_SIZE);
    r = sqrt(rho * rho + bottom_radius * bottom_radius);

    if (uvwz.z < 0.5)
    {
        // Distance to the ground for the ray (r,mu), and its minimum and maximum
        // values over all mu - obtained for (r,-1) and (r,mu_horizon) - from which
        // we can recover mu:
        float d_min = r - bottom_radius;
        float d_max = rho;
        float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
            1.0 - 2.0 * uvwz.z, SCATTERING_TEXTURE_MU_SIZE / 2);
        mu = d == 0.0 ? -1.0 : ClampCosine(-(rho * rho + d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = true;
    }
    else
    {
        // Distance to the top atmosphere boundary for the ray (r,mu), and its
        // minimum and maximum values over all mu - obtained for (r,1) and
        // (r,mu_horizon) - from which we can recover mu:
        float d_min = top_radius - r;
        float d_max = rho + H;
        float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(
            2.0 * uvwz.z - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2);
        mu = d == 0.0 ? 1.0 : ClampCosine((H * H - rho * rho - d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = false;
    }

    float x_mu_s =
        GetUnitRangeFromTextureCoord(uvwz.y, SCATTERING_TEXTURE_MU_S_SIZE);
    float d_min = top_radius - bottom_radius;
    float d_max = H;
    float D = DistanceToTopAtmosphereBoundary(bottom_radius, kMuSMin);
    float A = (D - d_min) / (d_max - d_min);
    float a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
    float d = d_min + min(a, A) * (d_max - d_min);
    mu_s = d == 0.0 ? 1.0 : ClampCosine((H * H - d * d) / (2.0 * bottom_radius * d));

    nu = ClampCosine(uvwz.x * 2.0 - 1.0);
}

float mod(float x, float y)
{
    return x - y * floor(x / y);
}

void GetRMuMuSNuFromScatteringTextureFragCoord(
    float3 frag_coord,
    out float r, out float mu, out float mu_s, out float nu,
    out bool ray_r_mu_intersects_ground)
{
    const float4 SCATTERING_TEXTURE_SIZE = float4(
        SCATTERING_TEXTURE_NU_SIZE - 1,
        SCATTERING_TEXTURE_MU_S_SIZE,
        SCATTERING_TEXTURE_MU_SIZE,
        SCATTERING_TEXTURE_R_SIZE);
    float frag_coord_nu =
        floor(frag_coord.x / float(SCATTERING_TEXTURE_MU_S_SIZE));
    float frag_coord_mu_s =
        mod(frag_coord.x, float(SCATTERING_TEXTURE_MU_S_SIZE));
    float4 uvwz =
        float4(frag_coord_nu, frag_coord_mu_s, frag_coord.y, frag_coord.z) /
        SCATTERING_TEXTURE_SIZE;
    GetRMuMuSNuFromScatteringTextureUvwz(
        uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    // Clamp nu to its valid range of values, given mu and mu_s.
    nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)),
        mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
}

#define USE_MULTISCATTERING_LUT 1

void ComputeSingleScattering(
        float r, float mu, float mu_s, float nu,
        bool ray_r_mu_intersects_ground,
        out float3 rayleigh, out float3 mie)
{
    DECLARE_RESOURCE_INDIRECT(Texture2D<float4>, transmittanceTex, cGlobal.atmosphere.transmittanceTex);
    DECLARE_RESOURCE_INDIRECT(Texture2D<float4>, multiscatteringTex, cGlobal.atmosphere.multiscatteringTex);

    // Number of intervals for the numerical integration.
    const int SAMPLE_COUNT = 50;
    // The integration step, i.e. the length of each integration interval.
    float dx = DistanceToNearestAtmosphereBoundary(r, mu,
            ray_r_mu_intersects_ground) / float(SAMPLE_COUNT);
    // Integration loop.
    float3 rayleigh_sum = (float3)(0.0);
    float3 mie_sum = (float3)(0.0);
    float3 transmittance = (float3)(1);
    for (int i = 0; i < SAMPLE_COUNT; ++i)
    {
        float d_i = (float(i) + 0.5) * dx;
        // The Rayleigh and Mie single scattering at the current sample point.
        float3 rayleigh_i;
        float3 mie_i;
    
        float d = d_i;
        float r_d = ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
        float mu_s_d = ClampCosine((r * mu_s + d * nu) / r_d);
        float3 transmittance_to_sun = GetSunVisibility(transmittanceTex, r_d, mu_s_d);
        
        float altitude_i = r_d - bottom_radius;
        float3 rayleigh_scattering_i;
        float3 mie_scattering_i;
        GetScatteringCoefficient(altitude_i, rayleigh_scattering_i, mie_scattering_i);
        
        float3 extinction_i = GetExtinctionCoefficient(altitude_i);
        float3 transmittance_i = exp(-extinction_i * dx);

#if USE_MULTISCATTERING_LUT
        float3 multiscattering_contribution = GetMultiscatteringContribution(multiscatteringTex, r_d, mu_s_d);
        multiscattering_contribution /= RayleighPhaseFunction(nu);
        rayleigh_i = rayleigh_scattering_i * (transmittance_to_sun + multiscattering_contribution)
            + mie_scattering_i * multiscattering_contribution;
#else
        rayleigh_i = rayleigh_scattering_i * transmittance_to_sun;
#endif
        mie_i = mie_scattering_i * transmittance_to_sun;

        // See slide 28 at https://www.frostbite.com/2015/08/physically-based-unified-volumetric-rendering-in-frostbite/
        rayleigh_sum += transmittance * (rayleigh_i - rayleigh_i * transmittance_i) / extinction_i;
        mie_sum += transmittance * (mie_i - mie_i * transmittance_i) / extinction_i;
        transmittance *= transmittance_i;
    }
    rayleigh = rayleigh_sum * solar_illuminance;
    mie = mie_sum * solar_illuminance;
}

void ComputeSingleScatteringTexture(float3 frag_coord, out float3 rayleigh, out float3 mie)
{
    float r;
    float mu;
    float mu_s;
    float nu;
    bool ray_r_mu_intersects_ground;
    GetRMuMuSNuFromScatteringTextureFragCoord(frag_coord,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ComputeSingleScattering(r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
}

[numthreads(8, 8, 1)]
void CSScattering(uint3 ThreadID : SV_DispatchThreadID)
{
    DECLARE_RESOURCE(RWTexture3D<float4>, scatteringRW)
    DECLARE_RESOURCE(RWTexture3D<float4>, singleMieScatteringRW)
    
    float3 frag_coord = float3(ThreadID) + 0.5;
    
    float3 delta_rayleigh, delta_mie;
    ComputeSingleScatteringTexture(
        frag_coord, delta_rayleigh, delta_mie);
    float4 scattering = float4(delta_rayleigh.rgb, delta_mie.r);
    scatteringRW[ThreadID] = scattering;
#if !(COMBINED_SCATTERING_TEXTURES && USE_MULTISCATTERING_LUT)
    singleMieScatteringRW[ThreadID] = float4(delta_mie, 1.0);
#endif
}
