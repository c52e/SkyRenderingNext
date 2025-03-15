#pragma once

#ifndef GLM_ENABLE_EXPERIMENTAL
#define GLM_ENABLE_EXPERIMENTAL
#endif

#include <glm/geometric.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

namespace Math
{
using namespace glm;
using float2 = vec2;
using float3 = vec3;
using float4 = vec4;
using float3x3 = mat3;
using float4x4 = mat4;
using int2 = ivec2;
using int3 = ivec3;
using int4 = ivec4;
using uint2 = uvec2;
using uint3 = uvec3;
using uint4 = uvec4;

// https://www.pbr-book.org/3ed-2018/Utilities/Memory_Management#Arena-BasedAllocation
template <class T> inline T align(T x, T alignment)
{
    auto mask = alignment - 1;
    return (x + mask) & (~mask);
}

inline float VanDerCorput(size_t base, size_t index)
{
    float ret = 0.0f;
    float denominator = float(base);
    while (index > 0)
    {
        size_t multiplier = index % base;
        ret += float(multiplier) / denominator;
        index = index / base;
        denominator *= base;
    }
    return ret;
}

inline float2 Halton23(uint32_t frameIndex, uint32_t mod)
{
    uint32_t index = (frameIndex % mod) + 1;
    return float2{VanDerCorput(2, index), VanDerCorput(3, index)} - 0.5f;
}

} // namespace Math