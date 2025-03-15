#pragma once

#include <filesystem>

#include "GraphicsDevice.h"
#include "ResourceBarrier.h"

void DataAddAlphaChannel(std::vector<UINT8> &dstVec, const UINT8 *src, int w, int h, int c);

class Texture2D : private ResourceBarrierBase<Texture2D>
{
    friend ResourceBarrierBase<Texture2D>;
  public:
    Texture2D() = default;

    void CreateFromData(GraphicsCommandList *cmd, UINT width, UINT height, const void *data, bool sRGB);

    void CreateFromPath(GraphicsCommandList *cmd, const std::filesystem::path &path);

    bool Empty() const noexcept
    {
        return !m_Resource;
    }

    ID3D12Resource *GetResource()
    {
        return m_Resource.Get();
    }

    const D3D12_RESOURCE_DESC &GetDesc() const noexcept
    {
        return m_ResourceDesc;
    }

    const DescriptorHandle &GetSrv() const noexcept
    {
        return m_Srv;
    }

    UINT GetSrvIndex() const noexcept
    {
        return GetSrv().index();
    }

  private:
    ResPtr<ID3D12Resource> m_Resource{};
    D3D12_RESOURCE_DESC m_ResourceDesc{};
    DescriptorHandle m_Srv{};
};

class BuildinTextures
{
  public:
    BuildinTextures(GraphicsCommandList *cmd);

    const Texture2D &GetWhite() const
    {
        return m_White;
    }

    const Texture2D &GetNormal() const
    {
        return m_Normal;
    }

    const Texture2D &GetTextureOrWhite(const Texture2D &tex) const
    {
        return tex.Empty() ? m_White : tex;
    }

    const Texture2D &GetTextureOrNormal(const Texture2D &tex) const
    {
        return tex.Empty() ? m_Normal : tex;
    }

  private:
    Texture2D m_White;
    Texture2D m_Normal;
};
