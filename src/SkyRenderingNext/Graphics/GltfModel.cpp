#include "GltfModel.h"

#include <unordered_set>

#include <tiny_gltf.h>

namespace fs = std::filesystem;

using Node = GltfModel::Node;
using Geometry = GltfModel::Geometry;

std::vector<GltfModel::Material> LoadMaterials(const tinygltf::Model &input)
{
    std::vector<GltfModel::Material> materials;
    materials.reserve(input.materials.size());
    for (const auto &material : input.materials)
    {
        GltfModel::Material m;
        if (auto itr = material.values.find("baseColorFactor"); itr != material.values.end())
        {
            m.baseColorFactor = Math::make_vec4(itr->second.ColorFactor().data());
        }
        if (auto itr = material.values.find("baseColorTexture"); itr != material.values.end())
        {
            m.baseColorIndex = itr->second.TextureIndex();
        }

        m.normalIndex = material.normalTexture.index;
        m.metallicRoughnessIndex = material.pbrMetallicRoughness.metallicRoughnessTexture.index;
        m.metallicFactor = static_cast<float>(material.pbrMetallicRoughness.metallicFactor);
        m.roughnessFactor = static_cast<float>(material.pbrMetallicRoughness.roughnessFactor);
        m.emissiveIndex = material.emissiveTexture.index;
        for (int j = 0; j < 3; ++j)
            m.emissiveFactor[j] = static_cast<float>(material.emissiveFactor[j]);
        m.aoIndex = material.occlusionTexture.index;
        m.aoFactor = static_cast<float>(material.occlusionTexture.strength);
        m.isOpaque = material.alphaMode == "OPAQUE";
        m.alphaCutoff = static_cast<float>(material.alphaCutoff);

        materials.push_back(m);
    }

    return materials;
}

std::vector<Texture2D> LoadImages(GraphicsCommandList *cmd, tinygltf::Model &input,
                                  const std::unordered_set<int> &srgbIndices,
                                  const std::unordered_set<int> &normalIndices)
{
    std::vector<Texture2D> images;
    images.reserve(input.images.size());
    for (auto &image : input.images)
    {
        if (image.bits < 0)
        {
            LOG_WARN("Invalid image");
            images.emplace_back();
            continue;
        }

        Texture2D texture;
        if (image.bits != 8)
        {
            LOG_WARN("Image with bits = {} is not supported yet. Use white instead. File name: {}", image.bits, image.name);
            
            UINT8 data[]{0xff, 0xff, 0xff, 0xff};
            texture.CreateFromData(cmd, 1, 1, data, true);
        }
        else
        {
            auto width = image.width;
            auto height = image.height;
            auto component = image.component;
            auto pbuffer = &image.image;
            std::vector<UINT8> rgbaBuffer;
            if (component == 3)
            {
                DataAddAlphaChannel(rgbaBuffer, pbuffer->data(), width, height, component);
                component = 4;
                pbuffer = &rgbaBuffer;
            }
            if (component != 4)
            {
                throw std::runtime_error(std::format("component = {} , except 4", component));
            }
            auto imageIndex = static_cast<int>(images.size());
            if (normalIndices.contains(imageIndex))
            {
                size_t numPixels = static_cast<size_t>(image.width) * static_cast<size_t>(image.height);
                for (size_t i = 0; i < numPixels; ++i)
                {
                    (*pbuffer)[i * 4 + 1] = 255 - (*pbuffer)[i * 4 + 1];
                }
            }

            auto sRGB = srgbIndices.contains(imageIndex);
            texture.CreateFromData(cmd, width, height, pbuffer->data(), sRGB);
        }

        images.push_back(std::move(texture));
    }

    return images;
}

std::vector<int> LoadTextures(const tinygltf::Model &input)
{
    // TODO: Sampler

    std::vector<int> textures;
    textures.reserve(input.textures.size());
    for (const auto &texture : input.textures)
    {
        textures.push_back(texture.source);
    }
    return textures;
}

std::vector<Node> LoadNodes(const tinygltf::Model &input)
{
    std::vector<Node> nodes;
    nodes.reserve(input.nodes.size());
    for (const auto &node : input.nodes)
    {
        Node n;
        n.children = node.children;
        n.mesh = node.mesh;

        if (node.translation.size() == 3)
        {
            n.translation = glm::make_vec3(node.translation.data());
        }
        if (node.rotation.size() == 4)
        {
            n.rotation = glm::make_quat(node.rotation.data());
#ifndef GLM_FORCE_QUAT_DATA_XYZW
            // https://github.com/g-truc/glm/commit/820a2c0e625f26000c688d841836bb10483be34d
            n.rotation = {n.rotation[3], n.rotation[0], n.rotation[1], n.rotation[2]};
#endif
        }
        if (node.scale.size() == 3)
        {
            n.scale = Math::make_vec3(node.scale.data());
        }
        if (node.matrix.size() == 16)
        {
            n.matrix = Math::make_mat4x4(node.matrix.data());
        }
        nodes.push_back(n);
    }
    return nodes;
}

struct MeshesData
{
    std::vector<GltfModel::Mesh> meshes;
    std::vector<ShaderResource::Vertex> vertexBuffer;
    std::vector<int> indexBuffer;
} LoadMeshes(const tinygltf::Model &input)
{
    MeshesData res;
    res.meshes.reserve(input.meshes.size());

    for (const auto &mesh : input.meshes)
    {
        GltfModel::Mesh m;
        m.geometries.reserve(mesh.primitives.size());
        for (const auto &primitive : mesh.primitives)
        {
            if (primitive.mode != TINYGLTF_MODE_TRIANGLES)
            {
                throw std::runtime_error("primitive.mode != TINYGLTF_MODE_TRIANGLES");
            }

            auto baseIndex = static_cast<uint32_t>(res.indexBuffer.size());
            auto vertexStart = static_cast<uint32_t>(res.vertexBuffer.size());

            // Vertices
            const float *positionBuffer = nullptr;
            const float *normalBuffer = nullptr;
            const float *uvBuffer = nullptr;
            const float *tangentBuffer = nullptr;
            uint32_t vertexCount = 0;

            for (auto [ppBuffer, key] : {
                     std::make_pair(&positionBuffer, "POSITION"),
                     std::make_pair(&normalBuffer, "NORMAL"),
                     std::make_pair(&uvBuffer, "TEXCOORD_0"),
                     std::make_pair(&tangentBuffer, "TANGENT"),
                 })
            {
                if (auto itr = primitive.attributes.find(key); itr != primitive.attributes.end())
                {
                    const auto &accessor = input.accessors[itr->second];
                    if (accessor.componentType != TINYGLTF_PARAMETER_TYPE_FLOAT)
                    {
                        throw std::runtime_error("accessor.componentType != TINYGLTF_PARAMETER_TYPE_FLOAT");
                    }
                    const auto &view = input.bufferViews[accessor.bufferView];
                    *ppBuffer = reinterpret_cast<const float *>(
                        &(input.buffers[view.buffer].data[accessor.byteOffset + view.byteOffset]));
                    vertexCount = static_cast<uint32_t>(accessor.count);
                }
            }

            if (!positionBuffer)
            {
                throw std::runtime_error("POSITION not found");
            }
            if (!normalBuffer)
            {
                throw std::runtime_error("NORMAL not found");
            }

            for (uint32_t v = 0; v < vertexCount; v++)
            {
                ShaderResource::Vertex vert{};
                vert.position = Math::make_vec3(&positionBuffer[v * 3]);
                vert.normal = Math::make_vec3(&normalBuffer[v * 3]);
                vert.uv = uvBuffer ? Math::make_vec2(&uvBuffer[v * 2]) : Math::vec2(0.0f);
                if (tangentBuffer)
                {
                    vert.tangent = Math::make_vec4(&tangentBuffer[v * 4]);
                }
                else
                {
                    const Math::float3 up(0.0f, 1.0f, 0.0f);
                    const Math::float3 right(1.0f, 0.0f, 0.0f);
                    // TODO: calculate tangents
                    auto tangent =
                        Math::normalize(Math::cross(vert.normal, abs(Math::dot(vert.normal, up)) > 0.01f ? right : up));;
                    vert.tangent = Math::float4(tangent, 1.0f);
                        
                }
                res.vertexBuffer.push_back(vert);
            }

            uint32_t indexCount = 0;
            if (primitive.indices >= 0)
            {
                // Indices
                const auto &accessor = input.accessors[primitive.indices];
                const auto &view = input.bufferViews[accessor.bufferView];
                const auto &buffer = input.buffers[view.buffer];
                indexCount = static_cast<uint32_t>(accessor.count);

                auto copy_data = [indexCount, &res](auto buf) {
                    for (uint32_t index = 0; index < indexCount; index++)
                    {
                        res.indexBuffer.push_back(buf[index]);
                    }
                };
                auto p_data = &buffer.data[accessor.byteOffset + view.byteOffset];
                switch (accessor.componentType)
                {
                case TINYGLTF_PARAMETER_TYPE_UNSIGNED_INT:
                    copy_data(reinterpret_cast<const uint32_t *>(p_data));
                    break;
                case TINYGLTF_PARAMETER_TYPE_UNSIGNED_SHORT:
                    copy_data(reinterpret_cast<const uint16_t *>(p_data));
                    break;
                case TINYGLTF_PARAMETER_TYPE_UNSIGNED_BYTE:
                    copy_data(reinterpret_cast<const uint8_t *>(p_data));
                    break;
                default:
                    throw std::runtime_error("Unknown index data type");
                }
            }

            Geometry p;
            p.materialIndex = primitive.material;
            p.vertexOffset = vertexStart;
            p.vertexCount = vertexCount;
            p.indexOffset = baseIndex;
            p.indexCount = indexCount;
            m.geometries.push_back(p);
        }
        res.meshes.push_back(std::move(m));
    }

    return res;
}

GltfModel::GltfModel(GraphicsCommandList *cmd, const std::filesystem::path &path)
{
    tinygltf::TinyGLTF gltfCtx;
    tinygltf::Model model;
    std::string err;
    std::string warn;
    std::string ext = path.extension().string();

    LOG_INFO("Loading {}", path.string());

    bool ret = false;
    if (ext.compare(".glb") == 0)
    {
        ret = gltfCtx.LoadBinaryFromFile(&model, &err, &warn, path.string());
    }
    else
    {
        ret = gltfCtx.LoadASCIIFromFile(&model, &err, &warn, path.string());
    }

    if (!warn.empty())
    {
        LOG_WARN("{}", warn);
    }

    if (!err.empty())
    {
        LOG_ERROR("{}", err);
    }

    if (!ret)
    {
        throw std::runtime_error(std::format("Loading {} failed", path.string()));
    }

    m_Materials = LoadMaterials(model);
    m_Textures = LoadTextures(model);

    std::unordered_set<int> srgbIndices;
    std::unordered_set<int> normalIndices;
    for (const auto &m : m_Materials)
    {
        for (auto tex_i : {m.baseColorIndex, m.emissiveIndex})
        {
            if (tex_i < 0)
                continue;
            auto img_i = m_Textures[tex_i];
            if (img_i < 0)
                continue;
            srgbIndices.insert(img_i);
        }
        if (m.normalIndex >= 0)
        {
            auto img_i = m_Textures[m.normalIndex];
            if (img_i >= 0)
            {
                normalIndices.insert(img_i);
            }
        }
    }

    m_Images = LoadImages(cmd, model, srgbIndices, normalIndices);
    m_Nodes = LoadNodes(model);
    if (model.scenes.size() == 0)
    {
        throw std::runtime_error("model.scenes.size() == 0");
    }
    m_SceneNodes = model.scenes[0].nodes;
    auto meshesData = LoadMeshes(model);
    const auto &vertexBuffer = meshesData.vertexBuffer;
    const auto &indexBuffer = meshesData.indexBuffer;

    auto &device = GetGraphicsDevice();
    auto d3d12Device = device.GetD3D12Device();

    m_IndexBuffer = StructuredBuffer::Create(cmd, indexBuffer);
    m_VertexBuffer = StructuredBuffer::Create(cmd, vertexBuffer);

    struct MeshBlasInfo
    {
        std::vector<D3D12_RAYTRACING_GEOMETRY_DESC> descs;
        D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_INPUTS inputs;
        D3D12_RAYTRACING_ACCELERATION_STRUCTURE_PREBUILD_INFO prebuildInfo;
    };
    std::vector<MeshBlasInfo> meshInfos;
    meshInfos.reserve(meshesData.meshes.size());
    // UINT64 scratchSize = 0;
    for (const auto &mesh : meshesData.meshes)
    {
        MeshBlasInfo meshInfo{};
        meshInfo.descs.reserve(mesh.geometries.size());
        for (const auto &geometry : mesh.geometries)
        {
            D3D12_RAYTRACING_GEOMETRY_DESC geometrydesc{
                .Type = D3D12_RAYTRACING_GEOMETRY_TYPE_TRIANGLES,
                .Flags = D3D12_RAYTRACING_GEOMETRY_FLAG_NONE,
                .Triangles{
                    .IndexFormat = DXGI_FORMAT_R32_UINT,
                    .VertexFormat = DXGI_FORMAT_R32G32B32_FLOAT,
                    .IndexCount = geometry.indexCount,
                    .VertexCount = geometry.vertexCount,
                    .IndexBuffer = m_IndexBuffer.GetGPUVirtualAddress(geometry.indexOffset * sizeof(indexBuffer[0])),
                    .VertexBuffer{
                        .StartAddress =
                            m_VertexBuffer.GetGPUVirtualAddress(geometry.vertexOffset * sizeof(vertexBuffer[0])),
                        .StrideInBytes = sizeof(vertexBuffer[0]),
                    },
                },
            };
            if (m_Materials[geometry.materialIndex].isOpaque)
            {
                geometrydesc.Flags |= D3D12_RAYTRACING_GEOMETRY_FLAG_OPAQUE;
            }
            meshInfo.descs.push_back(geometrydesc);
        }
        meshInfo.inputs = D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_INPUTS{
            .Type = D3D12_RAYTRACING_ACCELERATION_STRUCTURE_TYPE_BOTTOM_LEVEL,
            .NumDescs = static_cast<UINT>(meshInfo.descs.size()),
            .DescsLayout = D3D12_ELEMENTS_LAYOUT_ARRAY,
            .pGeometryDescs = meshInfo.descs.data(),
        };
        d3d12Device->GetRaytracingAccelerationStructurePrebuildInfo(&meshInfo.inputs, &meshInfo.prebuildInfo);
        // scratchSize = std::max(scratchSize, meshInfo.prebuildInfo.ScratchDataSizeInBytes);
        meshInfos.push_back(std::move(meshInfo));
    }

    /*auto scratchBuffer =
        device.CreateBuffer(scratchSize, D3D12_RESOURCE_STATE_COMMON, D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);*/

    m_BlasOffsets.reserve(meshInfos.size());
    std::vector<size_t> scratchOffsets;
    scratchOffsets.reserve(meshInfos.size());
    size_t blasesSize = 0;
    size_t scratchSize = 0;
    for (const auto &meshInfo : meshInfos)
    {
        m_BlasOffsets.push_back(blasesSize);
        scratchOffsets.push_back(scratchSize);
        blasesSize += Math::align(meshInfo.prebuildInfo.ResultDataMaxSizeInBytes, 256ull);
        scratchSize += Math::align(meshInfo.prebuildInfo.ScratchDataSizeInBytes, 256ull);
    }
    auto scratchBuffer =
        GpuBuffer(scratchSize, D3D12_RESOURCE_STATE_COMMON, D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);

    m_Blases = GpuBuffer(blasesSize, D3D12_RESOURCE_STATE_RAYTRACING_ACCELERATION_STRUCTURE,
                         D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);
    for (size_t i = 0; i < meshInfos.size(); ++i)
    {
        D3D12_BUILD_RAYTRACING_ACCELERATION_STRUCTURE_DESC blasBuildDesc{
            .DestAccelerationStructureData = GetBlasAddress(i),
            .Inputs = meshInfos[i].inputs,
            .ScratchAccelerationStructureData = scratchBuffer.GetGPUVirtualAddress(scratchOffsets[i]),
        };

        cmd->BuildRaytracingAccelerationStructure(&blasBuildDesc, 0, nullptr);
    }
    m_Blases.BarrierUAV(cmd);

    m_Meshes = std::move(meshesData.meshes);

    UpdateMatrices(Math::float4x4(1.0f));
}

void GltfModel::UpdateMatrices(const Math::float4x4 &m)
{
    m_Matrices.resize(m_Nodes.size());
    for (auto node : m_SceneNodes)
    {
        UpdateMatrices(m, node);
    }
}

void GltfModel::UpdateMatrices(const Math::float4x4 &m, int i)
{
    m_Matrices[i] = m * m_Nodes[i].LocalMatrix();
    for (auto child : m_Nodes[i].children)
    {
        UpdateMatrices(m_Matrices[i], child);
    }
}
