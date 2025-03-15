#ifndef _ATMOSPHERE_COMMON_H
#define _ATMOSPHERE_COMMON_H

#include "ShaderResources.hlsl"
#include "Common.hlsl"

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

#define STORE_OPTICAL_LENGTH_IN_TRANSMITTANCE_TEXTURE 1
#define COMBINED_SCATTERING_TEXTURES 0

static const float kMuSMin = cos(120.0 / 180.0 * PI);
static const float3 solar_illuminance = cGlobal.atmosphere.solarIlluminance;
static const float sun_angular_radius = cGlobal.atmosphere.sunAngularRadius;

static const float3 rayleigh_scattering = cGlobal.atmosphere.rayleighScattering;
static const float inv_rayleigh_exponential_distribution = cGlobal.atmosphere.invRayleighExponentialDistribution;

static const float3 mie_scattering = cGlobal.atmosphere.mieScattering;
static const float inv_mie_exponential_distribution = cGlobal.atmosphere.invMieExponentialDistribution;

static const float3 mie_absorption = cGlobal.atmosphere.mieAbsorption;
static const float ozone_center_altitude = cGlobal.atmosphere.ozoneCenterAltitude;

static const float3 ozone_absorption = cGlobal.atmosphere.ozoneAbsorption;
static const float inv_ozone_width = cGlobal.atmosphere.invOzoneWidth;

static const float3 ground_albedo = cGlobal.atmosphere.groundAlbedo;
static const float mie_phase_g = cGlobal.atmosphere.miePhaseG;

static const float bottom_radius = cGlobal.atmosphere.bottomRadius;
static const float top_radius = cGlobal.atmosphere.topRadius;
static const float transmittance_steps = cGlobal.atmosphere.transmittanceSteps;
static const float multiscattering_steps = cGlobal.atmosphere.multiscatteringSteps;

static const float3 earth_center = float3(0, -cGlobal.atmosphere.bottomRadius, 0);

#define sampler3D Texture3D<float4>
#define sampler2D Texture2D<float4>
#define texture(tex, uv) tex.Sample(bilinearSampler, uv)

uint2 textureSize(sampler2D tex)
{
    uint2 dim;
    tex.GetDimensions(dim.x, dim.y);
    return dim;
}

// https://ebruneton.github.io/precomputed_atmospheric_scattering/
// https://github.com/sebh/UnrealEngineSkyAtmosphere

float ClampCosine(float mu)
{
    return clamp(mu, -1.0, 1.0);
}

float ClampDistance(float d)
{
    return max(d, 0.0);
}

float ClampRadius(float r)
{
    return clamp(r, bottom_radius, top_radius);
}

float SafeSqrt(float a)
{
    return sqrt(max(a, 0.0));
}

float GetUnitRangeFromTextureCoord(float u, int texture_size)
{
    return (u - 0.5 / float(texture_size)) / (1.0 - 1.0 / float(texture_size));
}

float GetTextureCoordFromUnitRange(float x, int texture_size)
{
    return 0.5 / float(texture_size) + x * (1.0 - 1.0 / float(texture_size));
}

float2 GetTextureCoordFromUnitRange(float2 xy, int2 texture_size)
{
    return 0.5 / float2(texture_size) + xy * (1.0 - 1.0 / float2(texture_size));
}

bool RayIntersectsGround(float r, float mu)
{
    return mu < 0.0 && r * r * (mu * mu - 1.0) + bottom_radius * bottom_radius >= 0.0;
}

float DistanceToTopAtmosphereBoundary(float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1.0) + top_radius * top_radius;
    return ClampDistance(-r * mu + SafeSqrt(discriminant));
}

float DistanceToBottomAtmosphereBoundary(float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1.0) + bottom_radius * bottom_radius;
    return ClampDistance(-r * mu - SafeSqrt(discriminant));
}

float DistanceToNearestAtmosphereBoundary(
        float r, float mu, bool ray_r_mu_intersects_ground)
{
    return ray_r_mu_intersects_ground ?
        DistanceToBottomAtmosphereBoundary(r, mu) :
        DistanceToTopAtmosphereBoundary(r, mu);
}

float2 GetTransmittanceTextureUvFromRMu(int2 size, float r, float mu)
{
    // Distance to top atmosphere boundary for a horizontal ray at ground level.
    float H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon.
    float rho = SafeSqrt(r * r - bottom_radius * bottom_radius);
    // Distance to the top atmosphere boundary for the ray (r,mu), and its minimum
    // and maximum values over all mu - obtained for (r,1) and (r,mu_horizon).
    float d = DistanceToTopAtmosphereBoundary(r, mu);
    float d_min = top_radius - r;
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    return GetTextureCoordFromUnitRange(float2(x_mu, x_r), size);
}

float3 GetTransmittanceToTopAtmosphereBoundary(sampler2D transmittance_texture, float r, float mu)
{
    float2 uv = GetTransmittanceTextureUvFromRMu(textureSize(transmittance_texture), r, mu);
#if STORE_OPTICAL_LENGTH_IN_TRANSMITTANCE_TEXTURE
    return exp(-texture(transmittance_texture, uv).rgb);
#else
    return texture(transmittance_texture, uv).rgb;
#endif
}

#if STORE_OPTICAL_LENGTH_IN_TRANSMITTANCE_TEXTURE
float3 GetOpticalLengthToTopAtmosphereBoundary(sampler2D transmittance_texture, float r, float mu)
{
    float2 uv = GetTransmittanceTextureUvFromRMu(textureSize(transmittance_texture), r, mu);
    return texture(transmittance_texture, uv).rgb;
}
#endif

float3 GetTransmittance(sampler2D transmittance_texture,
        float r, float mu, float d, bool ray_r_mu_intersects_ground)
{

    float r_d = ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_d = ClampCosine((r * mu + d) / r_d);

#if STORE_OPTICAL_LENGTH_IN_TRANSMITTANCE_TEXTURE
    if (ray_r_mu_intersects_ground)
    {
        return min(
            exp(-GetOpticalLengthToTopAtmosphereBoundary(
                    transmittance_texture, r_d, -mu_d) +
                GetOpticalLengthToTopAtmosphereBoundary(
                    transmittance_texture, r, -mu)),
            (float3)(1.0));
    }
    else
    {
        return min(
            exp(-GetOpticalLengthToTopAtmosphereBoundary(
                transmittance_texture, r, mu) +
            GetOpticalLengthToTopAtmosphereBoundary(
                transmittance_texture, r_d, mu_d)),
            (float3)(1.0));
    }
#else
    if (ray_r_mu_intersects_ground) {
        return min(
            GetTransmittanceToTopAtmosphereBoundary(
                transmittance_texture, r_d, -mu_d) /
            GetTransmittanceToTopAtmosphereBoundary(
                transmittance_texture, r, -mu),
            (float3)(1.0));
    }
    else {
        return min(
            GetTransmittanceToTopAtmosphereBoundary(
                transmittance_texture, r, mu) /
            GetTransmittanceToTopAtmosphereBoundary(
                transmittance_texture, r_d, mu_d),
            (float3)(1.0));
    }
#endif
}

float3 GetSunVisibility(sampler2D transmittance_texture, float r, float mu_s)
{
    float sin_theta_h = bottom_radius / r;
    float cos_theta_h = -sqrt(max(1.0 - sin_theta_h * sin_theta_h, 0.0));
    return GetTransmittanceToTopAtmosphereBoundary(transmittance_texture, r, mu_s) *
        smoothstep(-sin_theta_h * sun_angular_radius,
            sin_theta_h * sun_angular_radius,
            mu_s - cos_theta_h);
}

float GetRayleighDensity(float altitude)
{
    return clamp(exp(-altitude * inv_rayleigh_exponential_distribution), 0, 1);
}

float GetMieDensity(float altitude)
{
    return clamp(exp(-altitude * inv_mie_exponential_distribution), 0, 1);
}

float GetOzoneDensity(float altitude)
{
    return max(0, altitude < ozone_center_altitude ?
            1 + (altitude - ozone_center_altitude) * inv_ozone_width :
            1 - (altitude - ozone_center_altitude) * inv_ozone_width);
}

float3 GetExtinctionCoefficient(float altitude)
{
    float3 rayleigh_extinction = rayleigh_scattering * GetRayleighDensity(altitude);

    float3 mie_extinction = (mie_scattering + mie_absorption) * GetMieDensity(altitude);

    float3 ozone_extinction = ozone_absorption * GetOzoneDensity(altitude);

    return rayleigh_extinction + mie_extinction + ozone_extinction;
}

float IsotropicPhaseFunction()
{
    return 1.0 / (4.0 * PI);
}

float RayleighPhaseFunction(float cos_theta)
{
    float k = 3.0 / (16.0 * PI);
    return k * (1.0 + cos_theta * cos_theta);
}

float CornetteShanks(float cos_theta, float g)
{
    float k = 3.0 / (8.0 * PI) * (1.0 - g * g) / (2.0 + g * g);
    return k * (1.0 + cos_theta * cos_theta) / pow(1.0 + g * g - 2.0 * g * cos_theta, 1.5);
}

float HenyeyGreenstein(float cos_theta, float g)
{
    float a = 1.0 - g * g;
    float b = 1.0 + g * g - 2.0 * g * cos_theta;
    b *= sqrt(b);
    return (0.25 * INV_PI) * a / b;
}

float MiePhaseFunction(float cos_theta, float g)
{
#if 0
    return CornetteShanks(cos_theta, g);
#else
    return HenyeyGreenstein(cos_theta, g);
#endif
}

float4 GetScatteringTextureUvwzFromRMuMuSNu(
        float r, float mu, float mu_s, float nu, bool ray_r_mu_intersects_ground)
{
    // Distance to top atmosphere boundary for a horizontal ray at ground level.
    float H = sqrt(top_radius * top_radius - bottom_radius * bottom_radius);
    // Distance to the horizon.
    float rho = SafeSqrt(r * r - bottom_radius * bottom_radius);
    float u_r = GetTextureCoordFromUnitRange(rho / H, SCATTERING_TEXTURE_R_SIZE);

    // Discriminant of the quadratic equation for the intersections of the ray
    // (r,mu) with the ground (see RayIntersectsGround).
    float r_mu = r * mu;
    float discriminant = r_mu * r_mu - r * r + bottom_radius * bottom_radius;
    float u_mu;
    if (ray_r_mu_intersects_ground)
    {
        // Distance to the ground for the ray (r,mu), and its minimum and maximum
        // values over all mu - obtained for (r,-1) and (r,mu_horizon).
        float d = -r_mu - SafeSqrt(discriminant);
        float d_min = r - bottom_radius;
        float d_max = rho;
        u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(d_max == d_min ? 0.0 :
            (d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
    }
    else
    {
        // Distance to the top atmosphere boundary for the ray (r,mu), and its
        // minimum and maximum values over all mu - obtained for (r,1) and
        // (r,mu_horizon).
        float d = -r_mu + SafeSqrt(discriminant + H * H);
        float d_min = top_radius - r;
        float d_max = rho + H;
        u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange(
            (d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
    }

    float d = DistanceToTopAtmosphereBoundary(bottom_radius, mu_s);
    float d_min = top_radius - bottom_radius;
    float d_max = H;
    float a = (d - d_min) / (d_max - d_min);
    float D = DistanceToTopAtmosphereBoundary(bottom_radius, kMuSMin);
    float A = (D - d_min) / (d_max - d_min);
    // An ad-hoc function equal to 0 for mu_s = mu_s_min (because then d = D and
    // thus a = A), equal to 1 for mu_s = 1 (because then d = d_min and thus
    // a = 0), and with a large slope around mu_s = 0, to get more texture 
    // samples near the horizon.
    float u_mu_s = GetTextureCoordFromUnitRange(
        max(1.0 - a / A, 0.0) / (1.0 + a), SCATTERING_TEXTURE_MU_S_SIZE);

    float u_nu = (nu + 1.0) / 2.0;
    return float4(u_nu, u_mu_s, u_mu, u_r);
}

float3 GetExtrapolatedSingleMieScattering(float4 scattering)
{
    // Algebraically this can never be negative, but rounding errors can produce
    // that effect for sufficiently short view rays.
    if (scattering.r <= 0.0)
    {
        return (float3)(0.0);
    }
    return scattering.rgb * scattering.a / scattering.r *
        (rayleigh_scattering.r / mie_scattering.r) *
        (mie_scattering / rayleigh_scattering);
}

float3 GetCombinedScattering(sampler3D scattering_texture,
        sampler3D single_mie_scattering_texture,
        float r, float mu, float mu_s, float nu,
        bool ray_r_mu_intersects_ground,
        out float3 single_mie_scattering)
{
    float4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(
        r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    float tex_coord_x = uvwz.x * float(SCATTERING_TEXTURE_NU_SIZE - 1);
    float tex_x = floor(tex_coord_x);
    float lerp = tex_coord_x - tex_x;
    float3 uvw0 = float3((tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
        uvwz.z, uvwz.w);
    float3 uvw1 = float3((tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE),
        uvwz.z, uvwz.w);
#if COMBINED_SCATTERING_TEXTURES
    float4 combined_scattering =
#if SCATTERING_TEXTURE_NU_SIZE == 1
        texture(scattering_texture, uvw0);
#else
        texture(scattering_texture, uvw0) * (1.0 - lerp) +
        texture(scattering_texture, uvw1) * lerp;
#endif
    float3 scattering = float3(combined_scattering);
    single_mie_scattering = GetExtrapolatedSingleMieScattering(combined_scattering);
#else
    float3 scattering = float3(
#if SCATTERING_TEXTURE_NU_SIZE == 1
        texture(scattering_texture, uvw0));
#else
        texture(scattering_texture, uvw0).rgb * (1.0 - lerp) +
        texture(scattering_texture, uvw1).rgb * lerp);
#endif
    single_mie_scattering = float3(
#if SCATTERING_TEXTURE_NU_SIZE == 1
        texture(single_mie_scattering_texture, uvw0));
#else
        texture(single_mie_scattering_texture, uvw0).rgb * (1.0 - lerp) +
        texture(single_mie_scattering_texture, uvw1).rgb * lerp);
#endif
#endif
    return scattering;
}

float3 GetSkyRadiance(sampler2D transmittance_texture,
        sampler3D scattering_texture, sampler3D single_mie_scattering_texture,
        float3 camera, float3 view_ray, float3 sun_direction, out float3 transmittance)
{
    // Compute the distance to the top atmosphere boundary along the view ray,
    // assuming the viewer is in space (or NaN if the view ray does not intersect
    // the atmosphere).
    float r = length(camera);
    float rmu = dot(camera, view_ray);

    // Compute the r, mu, mu_s and nu parameters needed for the texture lookups.
    float mu = rmu / r;
    float mu_s = dot(camera, sun_direction) / r;
    float nu = dot(view_ray, sun_direction);
    bool ray_r_mu_intersects_ground = RayIntersectsGround(r, mu);

    transmittance = ray_r_mu_intersects_ground ? (float3)(0.0) :
        GetTransmittanceToTopAtmosphereBoundary(transmittance_texture, r, mu);
    float3 single_mie_scattering;
    float3 scattering;
    scattering = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);
    return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
        MiePhaseFunction(nu, mie_phase_g);
}

float3 GetSkyRadianceToPoint(sampler2D transmittance_texture,
        sampler3D scattering_texture, sampler3D single_mie_scattering_texture,
        float3 camera, float3 view_ray, float d, float3 sun_direction, out float3 transmittance)
{
    // Compute the distance to the top atmosphere boundary along the view ray,
    // assuming the viewer is in space (or NaN if the view ray does not intersect
    // the atmosphere).
    float r = length(camera);
    float rmu = dot(camera, view_ray);
    float distance_to_top_atmosphere_boundary = -rmu -
        sqrt(rmu * rmu - r * r + top_radius * top_radius);

    // Compute the r, mu, mu_s and nu parameters for the first texture lookup.
    float mu = rmu / r;
    float mu_s = dot(camera, sun_direction) / r;
    float nu = dot(view_ray, sun_direction);
    bool ray_r_mu_intersects_ground = RayIntersectsGround(r, mu);

    transmittance = GetTransmittance(transmittance_texture,
        r, mu, d, ray_r_mu_intersects_ground);

    float3 single_mie_scattering;
    float3 scattering = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);

    // Compute the r, mu, mu_s and nu parameters for the second texture lookup.
    float r_p = ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_p = (r * mu + d) / r_p;
    float mu_s_p = (r * mu_s + d * nu) / r_p;

    float3 single_mie_scattering_p;
    float3 scattering_p = GetCombinedScattering(
        scattering_texture, single_mie_scattering_texture,
        r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
        single_mie_scattering_p);

    // Combine the lookup results to get the scattering between camera and point.
    scattering = scattering - transmittance * scattering_p;
    single_mie_scattering =
        single_mie_scattering - transmittance * single_mie_scattering_p;
#if COMBINED_SCATTERING_TEXTURES
    single_mie_scattering = GetExtrapolatedSingleMieScattering(
        float4(scattering, single_mie_scattering.r));
#endif

    // Hack to avoid rendering artifacts when the sun is below the horizon.
    single_mie_scattering = single_mie_scattering *
        smoothstep(float(0.0), float(0.01), mu_s);

    return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
        MiePhaseFunction(nu, mie_phase_g);
}


float3 ComputeGroundLuminance(sampler2D transmittance_texture, float3 earth_center, float3 position, float3 sun_direction)
{
    float3 up_direction = normalize(position - earth_center);
    float mu_s = dot(sun_direction, up_direction);
    float3 solar_illuminance_at_ground = GetSunVisibility(transmittance_texture, bottom_radius, mu_s);
#ifndef MULTISCATTERING_COMPUTE_PROGRAM
    solar_illuminance_at_ground *= solar_illuminance;
#endif
    float3 normal = normalize(position - earth_center);
    return INV_PI * clamp(dot(normal, sun_direction), 0, 1) * ground_albedo * solar_illuminance_at_ground;
}

#endif

