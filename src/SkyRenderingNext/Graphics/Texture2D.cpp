#include "Texture2D.h"

#include "stb_image.h"

#include "RenderBuffer.h"

void Texture2D::CreateFromData(GraphicsCommandList *cmd, UINT width, UINT height, const void *data, bool sRGB)
{
    const auto format = sRGB ? DXGI_FORMAT_R8G8B8A8_UNORM_SRGB : DXGI_FORMAT_R8G8B8A8_UNORM;
    m_ResourceDesc = D3D12_RESOURCE_DESC {
        .Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D,
        .Alignment = 0,
        .Width = width,
        .Height = height,
        .DepthOrArraySize = 1,
        .MipLevels = 1,
        .Format = format,
        .SampleDesc{.Count = 1, .Quality = 0},
        .Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN,
        .Flags = D3D12_RESOURCE_FLAG_NONE,
    };

    auto &device = GetGraphicsDevice();
    auto d3d12Device = device.GetD3D12Device();
    {
        CD3DX12_HEAP_PROPERTIES heapProps(D3D12_HEAP_TYPE_DEFAULT);
        m_ResourceState = D3D12_RESOURCE_STATE_COPY_DEST;
        d3d12Device->CreateCommittedResource(&heapProps, D3D12_HEAP_FLAG_NONE, &m_ResourceDesc, m_ResourceState,
                                             nullptr, IID_PPV_ARGS(&m_Resource));
    }

    const auto uploadBufferSize = GetRequiredIntermediateSize(m_Resource.Get(), 0, 1);
    UploadBuffer uploadBuffer(uploadBufferSize);
    D3D12_SUBRESOURCE_DATA textureData{
        .pData = data,
        .RowPitch = width * 4,
        .SlicePitch = width * height * 4,
    };
    UpdateSubresources(cmd->GetD3D12CmdList(), m_Resource.Get(), uploadBuffer.GetResource(), 0, 0, 1, &textureData);
    BarrierTransition(cmd, D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);

    m_Srv = device.AllocBindlessDescriptor();
    D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc{
        .Format = format,
        .ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D,
        .Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING,
        .Texture2D{
            .MostDetailedMip = 0,
            .MipLevels = 1,
            .PlaneSlice = 0,
            .ResourceMinLODClamp = 0,
        },
    };
    d3d12Device->CreateShaderResourceView(m_Resource.Get(), &srvDesc, device.GetBindlessCpuHandle(m_Srv));
}

void Texture2D::CreateFromPath(GraphicsCommandList *cmd, const std::filesystem::path &path)
{
    int w, h, comp;
    FILE *file{};
    _wfopen_s(&file, path.c_str(), L"rb");
    auto rawData = stbi_load_from_file(file, &w, &h, &comp, 0);
    fclose(file);
    ASSERT(rawData);
    std::vector<UINT8> imgData;
    void *pdata = rawData;
    if (comp == 3)
    {
        DataAddAlphaChannel(imgData, rawData, w, h, comp);
        comp = 4;
        pdata = imgData.data();
        stbi_image_free(rawData);
        rawData = nullptr;
    }
    ASSERT(comp == 4);
    CreateFromData(cmd, w, h, pdata, true);
    if (rawData)
        stbi_image_free(rawData);
}

BuildinTextures::BuildinTextures(GraphicsCommandList *cmd)
{
    {
        UINT8 data[]{0xff, 0xff, 0xff, 0xff};
        m_White.CreateFromData(cmd, 1, 1, data, true);
    }
    {
        uint8_t data[]{0x7f, 0x7f, 0xff, 0xff};
        m_Normal.CreateFromData(cmd, 1, 1, data, false);
    }
}

void DataAddAlphaChannel(std::vector<UINT8> &dstVec, const UINT8 *src, int w, int h, int c)
{
    dstVec.resize(static_cast<uint64_t>(w) * h * 4);
    auto dst = dstVec.data();
    for (size_t i = 0; i < static_cast<uint64_t>(w) * h; ++i)
    {
        *dst++ = *src++;
        *dst++ = *src++;
        *dst++ = *src++;
        *dst++ = 255;
    }
}
