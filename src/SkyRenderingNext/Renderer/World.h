#pragma once

#include <memory>

#include "entt/entt.hpp"

#include "Graphics/GltfModel.h"

using Entity = entt::entity;

struct Transform
{
    Math::float4x4 m;
};

struct EntityName
{
    std::string value;
};

struct ResourcePath
{
    std::filesystem::path value;
};

class World : public entt::registry
{
  public:
    friend World& GetWorld();
  private:
    World() = default;
};

using ModelResourceHandle = std::shared_ptr<GltfModel>;

struct GltfModelLoader
{
    using result_type = std::shared_ptr<GltfModel>;

    result_type operator()(GraphicsCommandList *cmd, const std::filesystem::path &path) const
    {
        return std::make_shared<GltfModel>(cmd, path);
    }
};

World &GetWorld();

entt::resource_cache<GltfModel, GltfModelLoader>& GetModelCache();

void AddModel(GraphicsCommandList *cmd, const std::filesystem::path &path, const Math::float4x4 &matrix);

bool EditWorld();
