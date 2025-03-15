#pragma once

#include "GraphicsDevice.h"
#include "Texture2D.h"
#include "RenderBuffer.h"
#include "ShaderResources.h"

class GltfModel
{
  public:
    struct Material
    {
        Math::float4 baseColorFactor = Math::float4(1.0f);
        Math::float3 emissiveFactor = Math::float3(0.0f);
        float metallicFactor = 1.0f;
        float roughnessFactor = 1.0f;
        float aoFactor = 1.0f;

        int baseColorIndex = -1;
        int normalIndex = -1;
        int metallicRoughnessIndex = -1;
        int emissiveIndex = -1;
        int aoIndex = -1;

        bool isOpaque = false;
        float alphaCutoff = 0.5f;
    };

    struct Node
    {
        int mesh{};
        std::vector<int> children;

        Math::float3 translation{};
        Math::float3 scale{1.0f};
        Math::quat rotation{};
        Math::float4x4 matrix = Math::float4x4(1.0f);

        Math::float4x4 LocalMatrix() const
        {
            return Math::translate(Math::float4x4(1.0f), translation) * Math::float4x4(rotation) *
                   Math::scale(Math::float4x4(1.0f), scale) * matrix;
        }
    };

    struct Geometry
    {
        int materialIndex;
        UINT vertexOffset;
        UINT vertexCount;
        UINT indexOffset;
        UINT indexCount;
    };

    struct Mesh
    {
        std::vector<Geometry> geometries;
    };

    struct GpuData
    {
        UINT rtas;
    };

    GltfModel(GraphicsCommandList *cmd, const std::filesystem::path &path);

    template <class Callback> void ForEachInstanceDesc(const Callback &callback)
    {
        for (size_t i = 0; i < m_Nodes.size(); ++i)
        {
            auto meshIndex = m_Nodes[i].mesh;
            if (meshIndex < 0)
            {
                continue;
            }
            callback(GetBlasAddress(meshIndex), m_Matrices[i]);
        }
    }

    template <class InstanceCallback, class GeometryCallback>
    void ForEachInstanceAndGeometry(const InstanceCallback &instanceCallback,
                                    const GeometryCallback &geometryCallback)
    {
        for (size_t i = 0; i < m_Nodes.size(); ++i)
        {
            auto meshIndex = m_Nodes[i].mesh;
            if (meshIndex < 0)
            {
                continue;
            }

            instanceCallback();
            for (const auto &g : m_Meshes[meshIndex].geometries)
            {
                geometryCallback(g);
            }
        }
    }

    template <class Callback> void ForEachMaterial(const Callback &callback)
    {
        auto &device = GetGraphicsDevice();
        auto getTexture = [this, &device](int index, bool isNormal = false) {
            const auto &buildinTextures = device.GetBuildinTextures();
            if (index < 0)
            {
                return isNormal ? buildinTextures.GetNormal().GetSrvIndex() : buildinTextures.GetWhite().GetSrvIndex();
            }
            return m_Images[m_Textures[index]].GetSrvIndex();
        };

        for (const auto &m : m_Materials)
        {
            callback(ShaderResource::Material{
                .baseColorFactor = m.baseColorFactor,
                .emissiveFactor = m.emissiveFactor,
                .metallicFactor = m.metallicFactor,
                .roughnessFactor = m.roughnessFactor,
                .aoFactor = m.aoFactor,
                .baseColorTex = getTexture(m.baseColorIndex),
                .normalTex = getTexture(m.normalIndex, true),
                .metallicRoughnessTex = getTexture(m.metallicRoughnessIndex),
                .emissiveTex = getTexture(m.emissiveIndex),
                .aoTex = getTexture(m.aoIndex),
                .alphaCutoff = m.alphaCutoff,
            });
        }
    }

    const StructuredBuffer& GetVertexBuffer() const
    {
        return m_VertexBuffer;
    }

    const StructuredBuffer &GetIndexBuffer() const
    {
        return m_IndexBuffer;
    }

  private:
    void UpdateMatrices(const Math::float4x4 &m);

    void UpdateMatrices(const Math::float4x4 &m, int i);

    D3D12_GPU_VIRTUAL_ADDRESS GetBlasAddress(size_t index)
    {
        return m_Blases.GetGPUVirtualAddress(m_BlasOffsets[index]);
    }

    std::vector<Math::float4x4> m_Matrices;

    std::vector<Material> m_Materials;
    std::vector<Texture2D> m_Images;
    std::vector<int> m_Textures;
    std::vector<Node> m_Nodes;
    std::vector<Mesh> m_Meshes;
    std::vector<int> m_SceneNodes;

    StructuredBuffer m_IndexBuffer;
    StructuredBuffer m_VertexBuffer;
    GpuBuffer m_Blases;
    std::vector<size_t> m_BlasOffsets;
};
