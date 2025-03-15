#include "RenderTexture.h"

inline bool IsDepthFormat(DXGI_FORMAT format)
{
    switch (format)
    {
    case DXGI_FORMAT_D16_UNORM:
    case DXGI_FORMAT_D24_UNORM_S8_UINT:
    case DXGI_FORMAT_D32_FLOAT:
    case DXGI_FORMAT_D32_FLOAT_S8X24_UINT:
        return true;
    default:
        return false;
    }
}

RenderTexture::RenderTexture(DXGI_FORMAT format, bool allowUav, std::wstring name)
    : RenderTexture({4, 4}, format, allowUav, std::move(name))
{
}

RenderTexture::RenderTexture(Math::uint2 size, const RenderTexture &fromTex)
    : RenderTexture(size, fromTex.GetDesc().Format, !fromTex.GetUav().Empty(), fromTex.m_Name)
{
}

RenderTexture::RenderTexture(Math::uint2 size, DXGI_FORMAT format, bool allowUav, std::wstring name)
    : ResourceBarrierBase(D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE), m_Name(std::move(name))
{
    m_ResourceDesc = D3D12_RESOURCE_DESC{
        .Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D,
        .Alignment = 0,
        .Width = size.x,
        .Height = size.y,
        .DepthOrArraySize = 1,
        .MipLevels = 1,
        .Format = format,
        .SampleDesc{.Count = 1, .Quality = 0},
        .Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN,
        .Flags =
            IsDepthFormat(format) ? D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL : D3D12_RESOURCE_FLAG_ALLOW_RENDER_TARGET,
    };

    if (allowUav)
    {
        m_ResourceDesc.Flags |= D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS;
    }

    CD3DX12_HEAP_PROPERTIES heapProps(D3D12_HEAP_TYPE_DEFAULT);
    auto &device = GetGraphicsDevice();
    auto d3d12Device = device.GetD3D12Device();
    d3d12Device->CreateCommittedResource(&heapProps, D3D12_HEAP_FLAG_NONE, &m_ResourceDesc, m_ResourceState, nullptr,
                                         IID_PPV_ARGS(&m_Resource));

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

    if (allowUav)
    {
        m_Uav = device.AllocBindlessDescriptor();
        D3D12_UNORDERED_ACCESS_VIEW_DESC uavDesc{
            .Format = format,
            .ViewDimension = D3D12_UAV_DIMENSION_TEXTURE2D,
            .Texture2D{
                .MipSlice = 0,
                .PlaneSlice = 0,
            },
        };
        d3d12Device->CreateUnorderedAccessView(m_Resource.Get(), nullptr, &uavDesc, device.GetBindlessCpuHandle(m_Uav));
    }

    m_Resource->SetName(m_Name.c_str());
}

RenderTexture::RenderTexture(Math::uint3 size, DXGI_FORMAT format, bool allowUav, std::wstring name)
    : ResourceBarrierBase(D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE), m_Name(std::move(name))
{
    m_ResourceDesc = D3D12_RESOURCE_DESC{
        .Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE3D,
        .Alignment = 0,
        .Width = size.x,
        .Height = size.y,
        .DepthOrArraySize = static_cast<UINT16>(size.z),
        .MipLevels = 1,
        .Format = format,
        .SampleDesc{.Count = 1, .Quality = 0},
        .Layout = D3D12_TEXTURE_LAYOUT_UNKNOWN,
        .Flags =
            IsDepthFormat(format) ? D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL : D3D12_RESOURCE_FLAG_ALLOW_RENDER_TARGET,
    };

    if (allowUav)
    {
        m_ResourceDesc.Flags |= D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS;
    }

    CD3DX12_HEAP_PROPERTIES heapProps(D3D12_HEAP_TYPE_DEFAULT);
    auto &device = GetGraphicsDevice();
    auto d3d12Device = device.GetD3D12Device();
    d3d12Device->CreateCommittedResource(&heapProps, D3D12_HEAP_FLAG_NONE, &m_ResourceDesc, m_ResourceState, nullptr,
                                         IID_PPV_ARGS(&m_Resource));

    m_Srv = device.AllocBindlessDescriptor();
    D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc{
        .Format = format,
        .ViewDimension = D3D12_SRV_DIMENSION_TEXTURE3D,
        .Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING,
        .Texture3D{
            .MostDetailedMip = 0,
            .MipLevels = 1,
            .ResourceMinLODClamp = 0,
        },
    };
    d3d12Device->CreateShaderResourceView(m_Resource.Get(), &srvDesc, device.GetBindlessCpuHandle(m_Srv));

    if (allowUav)
    {
        m_Uav = device.AllocBindlessDescriptor();
        D3D12_UNORDERED_ACCESS_VIEW_DESC uavDesc{
            .Format = format,
            .ViewDimension = D3D12_UAV_DIMENSION_TEXTURE3D,
            .Texture3D{
                .MipSlice = 0,
                .FirstWSlice = 0,
                .WSize = static_cast<UINT>(-1),
            },
        };
        d3d12Device->CreateUnorderedAccessView(m_Resource.Get(), nullptr, &uavDesc, device.GetBindlessCpuHandle(m_Uav));
    }

    m_Resource->SetName(m_Name.c_str());
}
