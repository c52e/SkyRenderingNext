#pragma once

#include <json.hpp>

#include "Common.h"

namespace nlohmann
{

template <glm::length_t L, typename T, glm::qualifier Q>
struct adl_serializer<glm::vec<L, T, Q>>
{
    static void to_json(json &json_j, const glm::vec<L, T, Q> &json_t)
    {
        for (glm::length_t i = 0; i < L; ++i)
        {
            json_j.push_back(json_t[i]);
        }
    }

    static void from_json(const json &json_j, glm::vec<L, T, Q> &json_t)
    {
        ASSERT(json_j.is_array() && json_j.size() == L);
        for (glm::length_t i = 0; i < L; ++i)
        {
            json_t[i] = json_j[i];
        }
    }
};

template <glm::length_t C, glm::length_t R, typename T, glm::qualifier Q> 
struct adl_serializer<glm::mat<C, R, T, Q>>
{
    static void to_json(nlohmann::json &json_j, const glm::mat<C, R, T, Q> &json_t)
    {
        for (glm::length_t i = 0; i < C; ++i)
        {
            json_j.push_back(json_t[i]);
        }
    }

    static void from_json(const nlohmann::json &json_j, glm::mat<C, R, T, Q> &json_t)
    {
        ASSERT(json_j.is_array() && json_j.size() == C);
        for (glm::length_t i = 0; i < C; ++i)
        {
            json_t[i] = json_j[i];
        }
    }
};

} // namespace nlohmann

namespace ser = nlohmann;
